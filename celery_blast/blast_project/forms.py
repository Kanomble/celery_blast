from django import forms
from django.core.exceptions import ValidationError
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User
from .py_biopython import get_species_taxid_by_name, check_given_taxonomic_node, get_list_of_species_taxid_by_name, \
    get_list_of_species_taxids_by_list_of_scientific_names, fetch_protein_records, \
    check_if_protein_identifier_correspond_to_backward_taxid
from .py_django_db_services import get_all_succeeded_databases, get_database_by_id, check_if_taxid_is_in_database, \
    check_if_sequences_are_in_database, check_if_project_title_exists, check_if_database_title_exists, check_if_remote_project_title_exists

from string import punctuation, ascii_letters

''' CreateTaxonomicFileForm
post form for the create_taxonomic_file.html template

functions:
__init__ :
    directly sets the user_email field with the request.user.email attribute
    this is used as input for the get_species_taxid_by_name function

clean_species_name : 
    translates the scientific name into a taxonomic node that can be used in the view
    for triggering the get_species_taxids.sh script (blast_project/tasks.get_species_taxids_into_file)
    the resulting taxonomic node is returned after invoking the .is_valid() method inside the clean_data dict
    
fields:
    species_name: filled by user, valid input should be able to successfully 
    user_email: user e-mail provided by passing the request.user object into form object creation, can get altered by the user
'''


class CreateTaxonomicFileForm(forms.Form):
    species_name = forms.CharField(max_length=200, required=True, )
    user_email = forms.CharField(max_length=200)

    def __init__(self, user, *args, **kwargs):
        super(CreateTaxonomicFileForm, self).__init__(*args, **kwargs)
        self.fields['user_email'].charfield = user.email
        self.fields['user_email'].initial = user.email

    def clean_species_name(self):
        species_name = self.cleaned_data['species_name']
        user_email = self.fields['user_email'].charfield
        try:
            taxonomic_nodes = get_list_of_species_taxid_by_name(user_email, species_name)
            return species_name, taxonomic_nodes
        except Exception as e:
            raise ValidationError(
                "validation error in clean_species_name pls check your provided scientific name : {}".format(e))


'''CreateTaxonomicFileForMultipleScientificNames

    This form replaces the old CreateTaxonomicFileForm.  In addition to all attributes from the previous form, this form
    inherits a filename field, which will be used as the actual filename, thus several species names will end up in too
    long names. This form can replace the CreateTaxonomicFileForm. 

'''


class CreateTaxonomicFileForMultipleScientificNames(forms.Form):
    filename = forms.CharField(max_length=200, required=True)
    species_names = forms.CharField(max_length=600, required=True)
    user_email = forms.CharField(max_length=200)

    def __init__(self, user, *args, **kwargs):
        super(CreateTaxonomicFileForMultipleScientificNames, self).__init__(*args, **kwargs)
        self.fields['user_email'].charfield = user.email
        self.fields['user_email'].initial = user.email

    def clean_species_names(self):
        species_names = self.cleaned_data['species_names']
        user_email = self.fields['user_email'].charfield

        try:
            species_names = species_names.split(",")
            taxids, errors = get_list_of_species_taxids_by_list_of_scientific_names(user_email, species_names)
            if len(errors) > 0:
                error_string = ''
                for name in errors:
                    error_string = error_string + name + ' '
                raise ValidationError(
                    "{}".format(error_string))
            return taxids
        except Exception as e:
            raise ValidationError(
                "validation error in clean_species_name pls check your provided scientific name : {}".format(e))


# registration and login form
class CreateUserForm(UserCreationForm):
    class Meta:
        model = User
        fields = ['username', 'email', 'password1', 'password2']


''' ProjectCreationForm
    This form is used in the form_project_creation.html.
    All form fields are processed into the BlastProject model instance via the project_creation_view functions.
    Field validation is done in the clean function of this form.
    
    
    fields | types:
        - project_title | charfield
        - search_strategy | choicefield --> currently disabled
        
        # The user has to provide data for at least one of these fields.
        - query_sequence_file | filefield
        - query_sequence_text | charfield
        
        - project_forward_database | BlastDatabaseModelChoiceField
        - project_backward_database | BlastDatabaseModelChoiceField
        - species_name_for_backward_blast | charfield
        - user_email | charfield
'''


class ProjectCreationForm(forms.Form):
    class BlastDatabaseModelChoiceField(forms.ModelChoiceField):
        def label_from_instance(self, blast_database):
            return str(blast_database.database_name)

    project_title = forms.CharField(
        label="Project title",
        error_messages={
            'required': "A project title is required for saving project metadata into the database"})

    # for now just blastp is possible this field is not included in the html form
    search_strategy = forms.ChoiceField(
        required=False,
        choices=(('blastp', 'blastp'), ('blastn', 'blastn')),
        label="Search strategy",
        error_messages={
            'required': "Please specify a search strategy otherwise you won't be able to execute a BLAST run .."})

    species_name_for_backward_blast = forms.CharField(
        required=True,
        label='Scientific Names (conversion to Taxonomic Nodes) for Backward BLAST',
        error_messages={
            'required': "Specify ONE! Scientific Name for your backward BLAST - The species from which you obtained the query sequences!"})

    user_email = forms.CharField(
        max_length=200,
        required=False)

    query_sequence_file = forms.FileField(
        required=False,
        error_messages={
            'required': "Upload a query sequence file, this file will serve as the -query parameter for the forward BLAST analysis"})

    query_sequence_text = forms.CharField(
        label="Query Sequence IDs", max_length=800, required=False
    )

    project_forward_database = BlastDatabaseModelChoiceField(
        queryset=get_all_succeeded_databases(),
        empty_label=None
    )

    project_backward_database = BlastDatabaseModelChoiceField(
        queryset=get_all_succeeded_databases(),
        empty_label=None
    )

    def __init__(self, user, *args, **kwargs):
        super(ProjectCreationForm, self).__init__(*args, **kwargs)
        self.fields['user_email'].charfield = user.email
        self.fields['user_email'].initial = user.email

    ''' Validation of user input.
        
        project_title: the project title is checked for uniqueness.
            
            functions: check_if_project_title_exists
        
        species_name_for_backward_blast: The provided species name must be a valid scientific name, that can get trans-
        lated into taxonomic identifier, which will be used for database filtering. This taxonomic identifier (taxid) 
        has to reside in the backward database. This may cause problems as there is another taxid called species taxid.
        Organism names are validated with biopython.
        
            functions: get_species_taxid_by_name, check_if_taxid_is_in_database
        
        query_sequence_file and/or query_sequence_text: One of those fields have to be filled by the user, even if both
        fields have the required=False setting. The filename has to consist of ascii_letters or "_", "-" signs. It has 
        to end with .faa, .fasta or .fa. Currently there are no restrictions for the number of provided sequences.
        Query sequences have to reside in the backward database, this is validated by combaring query sequence identifiers.
        
            functions: fetch_protein_records, check_if_sequences_are_in_database
    '''

    def clean(self):
        try:
            cleaned_data = super().clean()

            species_name = cleaned_data['species_name_for_backward_blast']
            user_email = self.fields['user_email'].charfield

            if check_if_project_title_exists(cleaned_data['project_title']):
                self.add_error('project_title', 'This title is already in use, please specify another title.')

            try:
                species_name = species_name.strip()
                taxonomic_nodes = get_species_taxid_by_name(user_email, species_name)
            except Exception as e:
                self.add_error('species_name_for_backward_blast',
                               'ERROR during fetching the taxonomic node of {}\n. There might be a entrez server error,'
                               ' try again later. Exception: {}.'.format(species_name, e))

            backward_db = cleaned_data['project_backward_database']
            try:
                booleanbw = check_if_taxid_is_in_database(backward_db.id, taxonomic_nodes)
                if booleanbw == False:
                    self.add_error('species_name_for_backward_blast',
                                   'specified taxonomic node: {} does not reside in the selected BACKWARD database: {}'.format(
                                       taxonomic_nodes, backward_db.database_name))

                if booleanbw == True:
                    # taking the first taxonomic node in the provided list
                    cleaned_data['species_name_for_backward_blast'] = (species_name, taxonomic_nodes[0])

            except Exception:
                self.add_error('species_name_for_backward_blast',
                               'specified taxonomic node: {} does not reside in the selected BACKWARD database: {}'.format(
                                   taxonomic_nodes, backward_db.database_name))

            query_file = cleaned_data['query_sequence_file']
            query_sequences = cleaned_data['query_sequence_text']

            # upload a query file or specify valid protein identifiers
            if query_file is None and query_sequences == '':
                self.add_error('query_sequence_file',
                               "please upload a fasta file containing your sequences or specify valid protein identifier")

            # query file was uploaded
            if query_file != None:
                if query_file.name.endswith('.faa') != True and query_file.name.endswith('.fasta') != True and query_file.name.endswith('.fa') != True:
                    self.add_error('query_sequence_file', "please upload only fasta files, your filename should end with .faa, .fa or .fasta")

                if len(query_file.name.split(".")) != 2:
                    self.add_error('query_sequence_file',
                                         "there are no dots allowed except the filetype delimiter")
                else:
                    filename = query_file.name.split(".")[0]
                    for character in filename:
                        if character in punctuation:
                            if character != '_' and character != '-':
                                self.add_error('query_sequence_file',
                                                     "bad character: \"{}\" in query file name".format(character))
                        if character not in ascii_letters:
                            if character != '_' and character != '-':
                                self.add_error('query_sequence_file',
                                                     "bad character: \"{}\" in query file name".format(character))

                header = []
                # checks accession identifier of query sequences
                for chunk in query_file.chunks():
                    lines = chunk.decode().split("\n")
                    for line in lines:
                        if line.startswith('>'):
                            try:
                                acc = line.split(" ")[0].split('>')[-1].split(".")[0]
                                if "|" in acc or " " in acc:
                                    raise Exception(
                                        "{} is no valid protein identifier - Format: e.g. WP_8765432".format(acc))
                                header.append(acc)
                            except Exception as e:
                                self.add_error('query_sequence_file',
                                               'error during parsing of query_file : {}'.format(e))

                if len(header) > 300:
                    self.add_error('query_sequence_file',
                                   'You try to infer orthologs for more than 300 query sequences,'
                                   ' this is not allowed, consider to separate the query sequences.')
                else:
                    valid = check_if_sequences_are_in_database(backward_db.id, header)
                    if valid != True:
                        self.add_error('query_sequence_file',
                                       'following sequences do not reside in your backward database: {}'.format(valid))

                if len(header) != len(set(header)):
                    self.add_error('query_sequence_file',
                                   'there are duplicate proteins in your uploaded file, please remove the duplicate entries and upload the file again!')


                #returncocde, protein = check_if_protein_identifier_correspond_to_backward_taxid(header, taxonomic_nodes[0], user_email)
                returncode = 0
                if returncode != 0:
                    protein = "placeholder"
                    self.add_error("species_name_for_backward_blast","specified taxonomic node: {} does not have any of the specified protein identifier(s) ...".format(taxonomic_nodes[0]))
                    self.add_error("query_sequence_file", "following protein identifier does not match to your specified taxonomic node: {}".format(protein))
            # protein identifier have been uploaded
            elif query_sequences != '' and query_file is None:
                # check string for invalid characters
                query_sequences = query_sequences.replace(" ", "").split(",")
                query_sequences = [qseq.split(".")[0] for qseq in query_sequences]
                try:
                    valid = check_if_sequences_are_in_database(backward_db.id, query_sequences)
                    if valid != True:
                        self.add_error('query_sequence_file',
                                       'following sequences do not reside in your backward database: {}'.format(valid))

                    proteins, errors = fetch_protein_records(query_sequences, user_email)
                    if len(errors) > 0:
                        self.add_error('query_sequence_text', 'following sequences are unavailable: {}'.format(errors))
                    cleaned_data['query_sequence_text'] = proteins
                except Exception as e:
                    self.add_error("query_sequence_text", "please provide valid protein identifiers. "
                                                          "Check your backward database. Maybe your backward database is broken.")

                # self.add_error('query_sequence_text','not available yet')
                #returncocde, protein = check_if_protein_identifier_correspond_to_backward_taxid(query_sequences,
                #                                                                                taxonomic_nodes[0],
                #                                                                                user_email)

                returncode = 0

                if returncode != 0:
                    self.add_error("species_name_for_backward_blast",
                                   "specified taxonomic node: {} does not have any of the specified protein identifier(s) ...".format(
                                       taxonomic_nodes[0]))
                    protein = "placeholder"
                    self.add_error("query_sequence_text",
                                   "following protein identifier does not match to your specified taxonomic node: {}".format(
                                       protein))

            else:
                self.add_error('query_sequence_text',
                               "please upload a fasta file containing your sequences or specify valid protein identifier")



        except Exception as e:
            raise ValidationError(
                "validation error in project creation, due to this exception: {}".format(
                    e))

        print("[***]", cleaned_data)
        return cleaned_data


'''BlastSettingsFormForward
    
    Form for the BlastProject. During project creation, a snakemake configuration file is written into the project dir-
    ectory. This form contains the fields for the forward BLAST settings, that are written into this configuration file.
    
'''


class BlastSettingsFormForward(forms.Form):
    fw_e_value = forms.FloatField(
        label="FW E-Value", initial=0.001)
    fw_word_size = forms.IntegerField(
        label="FW Word Size", initial=3)
    fw_num_alignments = forms.IntegerField(
        label="FW Number of possible alignment outputs", initial=10000)
    fw_max_target_seqs = forms.IntegerField(
        label="FW max_target_seqs of possible alignment description outputs", initial=10000)
    fw_num_threads = forms.IntegerField(
        label="FW Threads", initial=1)
    fw_max_hsps = forms.IntegerField(
        label='FW max hsps', initial=500
    )


'''BlastSettingsFormBackward

    Form for the BlastProject. During project creation, a snakemake configuration file is written into the project dir-
    ectory. This form contains the fields for the backward BLAST settings, that are written into this configuration file.
    
'''


class BlastSettingsFormBackward(forms.Form):
    bw_e_value = forms.FloatField(
        label="BW E-Value", initial=0.001)
    bw_word_size = forms.IntegerField(
        label="BW Word Size", initial=3)
    bw_num_alignments = forms.IntegerField(
        label="BW Number of possible alignment outputs", initial=1)
    bw_max_target_seqs = forms.IntegerField(
        label="BW max_target_seqs of possible alignment description outputs", initial=1)
    bw_num_threads = forms.IntegerField(
        label="BW Threads", initial=1)
    bw_max_hsps = forms.IntegerField(
        label='BW max hsps', initial=500
    )

'''SymBLASTProjectSettingsForm
    
    This form is used within the write_snakemake_configuration_file function in the models.py file.
    It is used to set up additional options for the snakemake pipeline.

'''
class SymBLASTProjectSettingsForm(forms.Form):
    bitscore_filter = forms.IntegerField(
        label="Bitscore threshold", initial=50
    )
    max_amount_of_rbh_for_msa_and_phylogeny = forms.IntegerField(
        label="Maximum number of sequences to use for phylogenetic inference", initial=500
    )
    # documentation of trimal options: http://trimal.cgenomics.org/use_of_the_command_line_trimal_v1.2
    trimal_gt = forms.FloatField(
        label="1 - (fraction of sequences with a gap allowed)", initial=0.8,
        min_value=0.0, max_value=1.0
    )
    trimal_st = forms.FloatField(
        label="Minimum average similarity allowed.", initial=0.001,
        min_value=0.0, max_value=1.0
    )
    trimal_cons = forms.IntegerField(
        label="Minimum percentage of the positions in the original alignment to conserve.", initial=60,
        min_value=0, max_value=100
    )
    # documentation of mview options: https://desmid.github.io/mview/manual/manual.html
    mview_sort = forms.ChoiceField(
        choices=(("coverage", "cov"), ("percent identity", "pid"), ("coverage then percent identity", "cov:pid"),
                 ("percent identity then coverage", "pid:cov"))
    )
    mview_coloring = forms.ChoiceField(
        choices=(("any","any"), ("identity","identity"), ("mismatch","mismatch"), ("consensus","consensus"),
                 ("group","group"))
    )


''' RemoteProjectCreationForm
    
    Adaption to ProjectCreationForm for remote BLAST projects. Model fields are essentially the 
    same as in the ProjectCreationForm except for the r_project_forward_database. All fields
    are extended by the prefix "r_".
    Options for the forward database are either the nr, the refseq or the swissprot database.

    This form is used in the form_remote_project_creation.html.
    All form fields are processed into the RemoteBlastProject model instance via the project_creation_view functions.
    Field validation is done in the clean function of this form.

    fields | types:
        - r_project_title | charfield
        - r_search_strategy | choicefield --> currently disabled

        # The user has to provide data for at least one of these fields.
        - r_query_sequence_file | filefield
        - r_query_sequence_text | charfield

        - r_project_forward_database | nr, swissprot or refseq
        - r_project_backward_database | BlastDatabaseModelChoiceField
        - r_species_name_for_backward_blast | charfield
        - r_user_email | charfield
'''


class RemoteProjectCreationForm(forms.Form):
    class BlastDatabaseModelChoiceField(forms.ModelChoiceField):
        def label_from_instance(self, blast_database):
            return str(blast_database.database_name)

    BLAST_REMOTE_DATABASES = [
        ('nr', 'none redundant proteins'),
        ('refseq_protein', 'refseq protein database'),
        ('refseq_select_prot', 'refseq selected protein database'),
        ('swissprot', 'protein sequences from the swissprot database')]

    r_project_title = forms.CharField(
        label="Project title",
        error_messages={
            'required': "A project title is required for saving project metadata into the database"})

    # for now just blastp is possible this field is not included in the html form
    r_search_strategy = forms.ChoiceField(
        required=False,
        choices=(('blastp', 'blastp'), ('blastn', 'blastn')),
        label="Search strategy",
        error_messages={
            'required': "Please specify a search strategy otherwise you won't be able to execute a BLAST run .."})

    r_species_name_for_backward_blast = forms.CharField(
        required=True,
        label='Scientific Names (conversion to Taxonomic Nodes) for Backward BLAST',
        error_messages={
            'required': "Specify ONE! Scientific Name for your backward BLAST - The species from which you obtained the query sequences!"})

    r_user_email = forms.CharField(
        max_length=200,
        required=False)

    r_query_sequence_file = forms.FileField(
        required=False,
        error_messages={
            'required': "Upload a query sequence file, this file will serve as the -query parameter for the forward BLAST analysis"})

    r_query_sequence_text = forms.CharField(
        label="Query Sequence IDs", max_length=800, required=False
    )

    r_project_forward_database = forms.ChoiceField(
        choices=BLAST_REMOTE_DATABASES
    )

    r_project_backward_database = BlastDatabaseModelChoiceField(
        queryset=get_all_succeeded_databases(),
        empty_label=None
    )

    r_entrez_query = forms.CharField(
        max_length=500,
        required=True,
        error_messages={'required':"Add taxonomic groups for which you would like to infer orthologs e.g.: eubacteria[organism] OR mammalia[organism].\n"
                                   "Check the syntax by clicking the link above.\nIf you want to search within the whole database add the query: all[organism]."}
    )

    def __init__(self, user, *args, **kwargs):
        super(RemoteProjectCreationForm, self).__init__(*args, **kwargs)
        self.fields['r_user_email'].charfield = user.email
        self.fields['r_user_email'].initial = user.email

    ''' Validation of user input.

        species_name_for_backward_blast: The provided species name must be a valid scientific name, that can get trans-
        lated into taxonomic identifier, which will be used for database filtering. This taxonomic identifier (taxid) 
        has to reside in the backward database. This may cause problems as there is another taxid called species taxid.
        Organism names are validated with biopython.

            functions: get_species_taxid_by_name, check_if_taxid_is_in_database

        query_sequence_file and/or query_sequence_text: One of those fields have to be filled by the user, even if both
        fields have the required=False setting. The filename has to consist of ascii_letters or "_", "-" signs. It has 
        to end with .faa, .fasta or .fa. Currently there are no restrictions for the number of provided sequences.
        Query sequences have to reside in the backward database, this is validated by combaring query sequence identifiers.

            functions: fetch_protein_records, check_if_sequences_are_in_database
    '''

    def clean(self):
        try:
            cleaned_data = super().clean()

            species_name = cleaned_data['r_species_name_for_backward_blast']
            user_email = self.fields['r_user_email'].charfield

            if check_if_remote_project_title_exists(cleaned_data['r_project_title']):
                self.add_error('r_project_title', 'This title is already in use, please specify another title.')

            try:
                species_name = species_name.strip()
                taxonomic_nodes = get_species_taxid_by_name(user_email, species_name)
            except Exception as e:
                self.add_error('r_species_name_for_backward_blast',
                               'ERROR during fetching the taxonomic node of {}\n. There might be a entrez server error,'
                               ' try again later. Exception: {}.'.format(species_name, e))

            backward_db = cleaned_data['r_project_backward_database']
            try:
                booleanbw = check_if_taxid_is_in_database(backward_db.id, taxonomic_nodes)
                if booleanbw == False:
                    self.add_error('r_species_name_for_backward_blast',
                                   'specified taxonomic node: {} does not reside in the selected BACKWARD database: {}'.format(
                                       taxonomic_nodes, backward_db.database_name))

                if booleanbw == True:
                    # taking the first taxonomic node in the provided list
                    cleaned_data['r_species_name_for_backward_blast'] = (species_name, taxonomic_nodes[0])

            except Exception:
                self.add_error('r_species_name_for_backward_blast',
                               'specified taxonomic node: {} does not reside in the selected BACKWARD database: {}'.format(
                                   taxonomic_nodes, backward_db.database_name))

            query_file = cleaned_data['r_query_sequence_file']
            query_sequences = cleaned_data['r_query_sequence_text']

            # upload a query file or specify valid protein identifiers
            if query_file is None and query_sequences == '':
                self.add_error('r_query_sequence_file',
                               "please upload a fasta file containing your sequences or specify valid protein identifier")

            # query file was uploaded
            if query_file != None:
                if query_file.name.endswith('.faa') != True and query_file.name.endswith('.fasta') != True:
                    self.add_error('r_query_sequence_file', "please upload only fasta files!")

                if len(query_file.name.split(".")) != 2:
                    self.add_error('r_query_sequence_file',
                                         "there are no dots allowed except the filetype delimiter")
                else:
                    filename = query_file.name.split(".")[0]
                    for character in filename:
                        if character in punctuation:
                            if character != '_' and character != '-':
                                self.add_error('r_query_sequence_file',
                                                     "bad character: \"{}\" in query file name".format(character))
                        if character not in ascii_letters:
                            if character != '_' and character != '-':
                                self.add_error('r_query_sequence_file',
                                                     "bad character: \"{}\" in query file name".format(character))

                header = []
                # checks accession identifier of query sequences
                for chunk in query_file.chunks():
                    lines = chunk.decode().split("\n")
                    for line in lines:
                        if line.startswith('>'):
                            try:
                                acc = line.split(" ")[0].split('>')[-1].split(".")[0]
                                if "|" in acc or " " in acc:
                                    raise Exception(
                                        "{} is no valid protein identifier - Format: e.g. WP_8765432".format(acc))
                                header.append(acc)
                            except Exception as e:
                                self.add_error('r_query_sequence_file',
                                               'error during parsing of query_file : {}'.format(e))

                if len(header) > 20:
                    self.add_error('r_query_sequence_file',
                                   'You try to infer orthologs for more than 20 query sequences,'
                                   ' this is not allowed with remote BLAST searches, consider to separate the query sequences or use the local pipeline.')
                else:
                    valid = check_if_sequences_are_in_database(backward_db.id, header)
                    if valid != True:
                        self.add_error('r_query_sequence_file',
                                       'following sequences do not reside in your backward database: {}'.format(valid))

                if len(header) != len(set(header)):
                    self.add_error('r_query_sequence_file',
                                   'there are duplicate proteins in your uploaded file, please remove the duplicate entries and upload the file again!')

                #returncocde, protein = check_if_protein_identifier_correspond_to_backward_taxid(header,
                #                                                                                taxonomic_nodes[0],
                #                                                                                user_email)
                returncode = 0

                if returncode != 0:
                    self.add_error("r_species_name_for_backward_blast",
                                   "specified taxonomic node: {} does not have any of the specified protein identifier(s) ...".format(
                                       taxonomic_nodes[0]))

                    protein = "placeholder"
                    self.add_error("r_query_sequence_file",
                                   "following protein identifier does not match to your specified taxonomic node: {}".format(
                                       protein))
            # protein identifier have been uploaded
            elif query_sequences != '' and query_file is None:
                # check string for invalid characters
                query_sequences = query_sequences.replace(" ", "").split(",")
                query_sequences = [qseq.split(".")[0] for qseq in query_sequences]
                try:
                    valid = check_if_sequences_are_in_database(backward_db.id, query_sequences)
                    if valid != True:
                        self.add_error('r_query_sequence_file',
                                       'following sequences do not reside in your backward database: {}'.format(valid))

                    proteins, errors = fetch_protein_records(query_sequences, user_email)
                    if len(errors) > 0:
                        self.add_error('r_query_sequence_text', 'following sequences are unavailable: {}'.format(errors))
                    cleaned_data['r_query_sequence_text'] = proteins
                except Exception as e:
                    self.add_error("query_sequence_text", "please provide valid protein identifiers. "
                                                          "Maybe your backward database is broken, check your backward database.")

                # self.add_error('query_sequence_text','not available yet')
                returncocde, protein = check_if_protein_identifier_correspond_to_backward_taxid(query_sequences,
                                                                                                taxonomic_nodes[0],
                                                                                                user_email)
                if returncocde != 0:
                    self.add_error("r_species_name_for_backward_blast",
                                   "specified taxonomic node: {} does not have any of the specified protein identifier(s) ...".format(
                                       taxonomic_nodes[0]))
                    self.add_error("r_query_sequence_text",
                                   "following protein identifier does not match to your specified taxonomic node: {}".format(
                                       protein))

            else:
                self.add_error('r_query_sequence_text',
                               "please upload a fasta file containing your sequences or specify valid protein identifier")

        except Exception as e:
            raise ValidationError(
                "validation error in project creation, due to this exception: {}".format(
                    e))
        return cleaned_data

'''UploadGenomeForm

    Form for uploading a single (concatenated) fasta file. Several metadata information are necessary to format this file
    into a BlastDatabase. Those information are reflected by this form fields. Validation is done via a longer clean 
    function:
    
    taxonomic_node: The user has to provide a taxonomic_node or a file containing a mapping between taxonomic_nodes and
    sequence identifier (taxmap_file). Taxonomic nodes are validated with biopython.
        
        functions: check_given_taxonomic_node
    
    If the user uploads a concatenated fastafile, all metadata have toget uploaded via files. Those files should contain
    information separated by newlines, e.g. each organism name has to reside in their one line. The length of those files
    have to be the same. If the user provides 9 different organism names we can assume, that there must be 9 different
    assemblies in the uploaded fasta file. Organism names are validated with biopython.
    
        functions: get_species_taxid_by_name
    
           
'''


class UploadGenomeForm(forms.Form):
    database_title = forms.CharField(max_length=200, required=True)
    database_description = forms.CharField(max_length=200, required=True)
    assembly_entries = forms.IntegerField(min_value=1, required=True)
    # obsolete due to multiple genome form for just one genome file
    # TODO remove fields and fix associated view/form functions
    ###### OBSOLETE ######
    assembly_accession = forms.CharField(max_length=200, required=False)
    assembly_level = forms.ChoiceField(
        choices=(("Chromosome", "Chromosome"), ("Complete Genome", "Complete Genome"), ("Scaffold", "Scaffold"),
                 ("Contig", "Contig")),
        required=False)

    organism_name = forms.CharField(
        label="Organism Name", required=False
    )
    ###### OBSOLETE ######

    taxonomic_node = forms.IntegerField(
        min_value=2,
        label="Taxonomic ID", required=False)

    taxmap_file = forms.FileField(
        required=False
    )

    organism_name_file = forms.FileField(
        required=False
    )

    assembly_accessions_file = forms.FileField(
        required=False
    )

    assembly_level_file = forms.FileField(
        required=False
    )

    genome_fasta_file = forms.FileField(required=True)

    user_email = forms.CharField(max_length=200, required=False)

    def __init__(self, user, *args, **kwargs):
        super(UploadGenomeForm, self).__init__(*args, **kwargs)
        self.fields['user_email'].charfield = user.email
        self.fields['user_email'].initial = user.email

    def clean(self):
        cleaned_data = super(UploadGenomeForm, self).clean()
        genome_fasta_file = cleaned_data['genome_fasta_file']
        taxmap_file = cleaned_data['taxmap_file']
        taxonomic_node = cleaned_data['taxonomic_node']
        organism_file = cleaned_data['organism_name_file']
        assembly_accessions_file = cleaned_data['assembly_accessions_file']
        assembly_level_file = cleaned_data['assembly_level_file']
        user_email = cleaned_data['user_email']

        if 'database_title' not in list(cleaned_data.keys()):
            self.add_error('database_title', "This field is required, set an appropriate database title!")
        else:
            database_title = cleaned_data['database_title']
            if check_if_database_title_exists(database_title):
                self.add_error('database_title', 'This title is already in use, please specify another title.')

        if genome_fasta_file.name.endswith('.faa') is not True and genome_fasta_file.name.endswith(
                '.fasta') is not True:
            self.add_error('genome_fasta_file', 'specify a valid fasta file (with .faa or .fasta file name extension)')

        if taxonomic_node is not None:
            try:
                check_given_taxonomic_node(user_email, taxonomic_node)
            except Exception as e:
                self.add_error('taxonomic_node',
                               'there is no organism with your specified taxonomic node, exception : {}'.format(e))

        if taxmap_file is None and taxonomic_node is None:
            self.add_error('taxonomic_node', 'specify a taxonomic node for the formatting procedure')
            self.add_error('taxmap_file', 'specify a taxmap file for the formatting procedure')

        if taxmap_file is not None:

            if taxonomic_node is not None:
                self.add_error('taxonomic_node', 'just specify one, a taxonomic node or a taxmap file')
                self.add_error('taxmap_file', 'just specify one, a taxonomic node or a taxmap file')

            if organism_file is None:
                self.add_error('organism_name_file',
                               'if you upload a taxmap file you should also upload a file that contains the organism names of target genomes, separated by lines')

            protein_ids = 0
            for chunk in genome_fasta_file.chunks():
                protein_ids += chunk.decode().count(">")


            taxmap_ids = 0
            for chunk in taxmap_file.chunks():
                taxmap_ids += chunk.decode().count('\n')

            # check if taxonomic_nodes_exists
            organisms = 0
            for chunk in organism_file.chunks():

                lines = chunk.decode().split("\n")
                for line in lines:
                    if line != '':
                        try:
                            # check if there are taxids available (for provided organism names)
                            taxids = get_species_taxid_by_name(user_email, line)
                            # multiple taxids are valid
                            if len(taxids) == 0:
                                raise Exception
                        except:
                            self.add_error('organism_name_file',
                                           'there is no taxonomic node available for : {}'.format(line))
                organisms += chunk.decode().count('\n')

            if taxmap_ids != protein_ids:
                self.add_error('taxmap_file',
                               'the amount of provided acc_ids: {} does not match the provided amount of protein_ids: {}, '
                               'make sure there are no additional ">" signs within the fasta header'.format(
                                   taxmap_ids, protein_ids))

            if assembly_accessions_file is not None:
                amount_of_assemblies = 0
                for chunk in assembly_accessions_file.chunks():
                    lines = chunk.decode().split('\n')
                    for line in lines:
                        if line == '\r':
                            self.add_error('assembly_accessions_file',
                                           'there are lines without any informations in your assembly accessions file')
                    amount_of_assemblies += chunk.decode().count('\n')
                if amount_of_assemblies != organisms:
                    self.add_error('assembly_accessions_file',
                                   'the amount of assemblies: {} does not match the amount of provided organisms: {}'.format(
                                       amount_of_assemblies, organisms))

            if assembly_level_file is not None:
                amount_of_levels = 0
                for chunk in assembly_level_file.chunks():
                    amount_of_levels += chunk.decode().count('\n')
                if amount_of_levels != organisms:
                    self.add_error('assembly_level_file',
                                   'the amount of assembly levels: {} does not match the amount of provided organisms: {}'.format(
                                       amount_of_levels, organisms))
        return cleaned_data


'''UploadMultipleFilesGenomeForm

    Form for uploading multiple fasta files. The user has to provide a valid organism name, that can get translated into
    a taxonomic identifier. This form can replace the UploadGenomeForm, as it has a similar functionality. 

'''


class UploadMultipleFilesGenomeForm(forms.Form):
    database_title = forms.CharField(max_length=200, required=True)
    database_description = forms.CharField(max_length=200, required=True)

    genome_file_field_0 = forms.FileField(required=True)
    organism_name_0 = forms.CharField(max_length=200, required=True)
    extra_field_count = forms.CharField(initial="0", widget=forms.HiddenInput())
    user_email = forms.CharField(max_length=200, required=False)

    def __init__(self, user, *args, **kwargs):
        extra_fields = kwargs.pop('extra', 0)
        super(UploadMultipleFilesGenomeForm, self).__init__(*args, **kwargs)

        self.fields['user_email'].charfield = user.email
        self.fields['user_email'].initial = user.email
        self.fields['extra_field_count'].initial = extra_fields
        extra_fields = int(extra_fields)

        if extra_fields > 0:
            extra_fields += 1
        for index in range(extra_fields):
            self.fields['genome_file_field_{index}'.format(index=index)] = forms.FileField(required=False)
            self.fields['organism_name_{index}'.format(index=index)] = forms.CharField(required=False)

    # https://stackoverflow.com/questions/6142025/dynamically-add-field-to-a-form

    def clean(self):
        cleaned_data = super().clean()
        user_email = cleaned_data['user_email']
        if 'database_title' not in list(cleaned_data.keys()):
            self.add_error('database_title', "This field is required, set an appropriate database title!")
        else:
            database_title = cleaned_data['database_title']
            if check_if_database_title_exists(database_title):
                self.add_error('database_title', 'This title is already in use, please specify another title.')

        for field in self.fields:
            if "genome_file" in field:
                file = cleaned_data.get(field)
                if file == None:
                    self.add_error(field, 'You have to provide a genome file')

                elif file.name.split(".")[-1] in ["fasta", "faa", "fa"] == False:
                    self.add_error(field,
                                   'You have to upload a FASTA file, if you provide a valid FASTA file make sure to have a file ending with .faa, .fasta or .fa!')

            elif "organism" in field:
                if cleaned_data.get(field) == '':
                    self.add_error(field, 'You have to provide a valid scientific name')
                elif cleaned_data.get(field) == None:
                    self.add_error(field, "You have to provide a valid scientific name")
                else:
                    try:
                        taxids = get_species_taxid_by_name(user_email, cleaned_data.get(field))
                        # multiple taxids are valid
                        if len(taxids) == 0:
                            raise Exception
                    except Exception as e:
                        self.add_error(field, "{} : {} is no valid name!".format(e, cleaned_data.get(field)))

        return cleaned_data
