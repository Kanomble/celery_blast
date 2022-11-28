from django import forms
from django.core.exceptions import ValidationError
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User
from .py_biopython import get_species_taxid_by_name, check_given_taxonomic_node, get_list_of_species_taxid_by_name, \
    get_list_of_species_taxids_by_list_of_scientific_names, check_if_protein_identifier_correspond_to_backward_taxid
from .py_django_db_services import get_all_succeeded_databases, get_database_by_id, check_if_taxid_is_in_database, check_if_sequences_are_in_database

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

    species_name = forms.CharField(max_length=200, required=True,)
    user_email = forms.CharField(max_length=200)

    def __init__(self, user, *args, **kwargs):
        super(CreateTaxonomicFileForm, self).__init__(*args, **kwargs)
        self.fields['user_email'].charfield = user.email
        self.fields['user_email'].initial = user.email

    def clean_species_name(self):
        species_name = self.cleaned_data['species_name']
        user_email = self.fields['user_email'].charfield
        try:
            taxonomic_nodes = get_list_of_species_taxid_by_name(user_email,species_name)
            return species_name, taxonomic_nodes
        except Exception as e:
            raise ValidationError("validation error in clean_species_name pls check your provided scientific name : {}".format(e))

'''CreateTaxonomicFileForMultipleScientificNames
This form replaces the old CreateTaxonomicFileForm. 
In addition to all attributes from the previous form, this form
inherits a filename field, which will be used as the actual filename,
 thus several species names will end up in too long names.

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
            raise ValidationError("validation error in clean_species_name pls check your provided scientific name : {}".format(e))

#registration and login form
class CreateUserForm(UserCreationForm):
    class Meta:
        model = User
        fields = ['username','email','password1','password2']

''' ProjectCreationForm
    This form is used in the form_project_creation.html.
    All form fields are processed into the BlastProject model instance via the project_creation_view functions.
    Field validation is done in the clean function of this form.
    
    
    fields | types:
        - project_title | charfield
        - search_strategy | choicefield --> currently disabled
        - query_sequence_file | filefield
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

    #for now just blastp is possible this field is not included in the html form
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
            'required':"Specify a Scientific Name for your backward BLAST - use a comma separated list - names will be converted to taxids that will be written to a file which will serve as the -taxidlist parameter of your backward BLAST"})

    user_email = forms.CharField(
        max_length=200,
        required=False)

    query_sequence_file = forms.FileField(
        error_messages={
            'required':"Upload a query sequence file, this file will serve as the -query parameter for the forward BLAST analysis"})

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

    def clean(self):
        try:
            cleaned_data = super().clean()

            species_name = cleaned_data['species_name_for_backward_blast']
            user_email = self.fields['user_email'].charfield
            try:
                taxonomic_nodes = get_species_taxid_by_name(user_email, species_name)
            except Exception as e:
                self.add_error('species_name_for_backward_blast','your provided species has no taxonomic node - check on NCBI'.format(species_name))

            backward_db = cleaned_data['project_backward_database']
            try:
                booleanbw = check_if_taxid_is_in_database(backward_db.id, taxonomic_nodes)
                if booleanbw == False:
                    self.add_error('species_name_for_backward_blast',
                                   'specified taxonomic node: {} does not reside in the selected BACKWARD database: {}'.format(
                                       taxonomic_nodes, backward_db.database_name))

                if booleanbw == True:
                    #taking the first taxonomic node in the provided list
                    cleaned_data['species_name_for_backward_blast'] = (species_name, taxonomic_nodes[0])

            except Exception:
                self.add_error('species_name_for_backward_blast',
                               'specified taxonomic node: {} does not reside in the selected BACKWARD database: {}'.format(
                                   taxonomic_nodes, backward_db.database_name))


            query_file = cleaned_data['query_sequence_file']

            #check for fasta files
            if query_file.name.endswith('.faa') != True and query_file.name.endswith('.fasta') != True:
                self.add_error("query_sequence_file","please upload only fasta files!")
            else:
                header = []
                #checks accession identifier of query sequences
                for chunk in query_file.chunks():
                    lines = chunk.decode().split("\n")
                    for line in lines:
                        if line.startswith('>'):
                            try:
                                acc = line.split(" ")[0].split('>')[-1].split(".")[0]
                                if "|" in acc or " " in acc:
                                    raise Exception("{} is no valid protein identifier - Format: e.g. WP_8765432".format(acc))
                                header.append(acc)
                            except Exception as e:
                                self.add_error('query_sequence_file','error during parsing of query_file : {}'.format(e))
                #maximum number of query sequences = 300
                if len(header) > 300:
                    self.add_error('query_sequence_file','You try to infer orthologs for more than 300 query sequences,'
                                                         ' this is not allowed, consider to separate the query sequences.')

                #checks if query sequences reside in the backward database
                else:
                    valid = check_if_sequences_are_in_database(backward_db.id, header)
                    if valid != True:
                        self.add_error('query_sequence_file','following sequences do not reside in your backward database: {}'.format(valid))
                #check for duplicate entries in the query file
                if len(header) != len(set(header)):
                    self.add_error('query_sequence_file','there are duplicate proteins in your uploaded file, please remove the duplicate entries and upload the file again!')

                #check if provided taxid corresponds to query sequence origin
                #what if the user provides custom query sequence ids?
                #try:
                #    retcode=check_if_protein_identifier_correspond_to_backward_taxid(header,taxonomic_nodes[0],user_email)
                #    if retcode != 0:
                #        self.add_error('species_name_for_backward_blast',
                #                       'specified taxonomic node: {} does not match with query sequences, following nodes have been translated from the protein queries: {}'.format(
                #                           taxonomic_nodes[0],' '.join(retcode)))
                #except Exception:
                #    #TODO add taxid for query sequences
                #    self.add_error('species_name_for_backward_blast',
                #                   'specified taxonomic node: {} does not match with query sequences'.format(
                #                       taxonomic_nodes[0]))

        except Exception as e:
            raise ValidationError(
                "validation error in project creation, due to this exception: {}".format(
                     e))
        return cleaned_data

#TODO documentation
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

#TODO documentation
class BlastSettingsFormBackward(forms.Form):
    bw_e_value = forms.FloatField(
        label="BW E-Value", initial=0.001)
    bw_word_size = forms.IntegerField(
        label="BW Word Size", initial=3)
    bw_num_alignments = forms.IntegerField(
        label="BW Number of possible alignment outputs",initial=1)
    bw_max_target_seqs = forms.IntegerField(
        label="BW max_target_seqs of possible alignment description outputs", initial=1)
    bw_num_threads = forms.IntegerField(
        label="BW Threads",initial=1)
    bw_max_hsps = forms.IntegerField(
        label='BW max hsps', initial=500
    )

class UploadGenomeForm(forms.Form):
    genome_fasta_file = forms.FileField(
        error_messages={
                'required':"Upload a genome FASTA file with protein sequences, that can get formatted to a BLAST database."}
    )

    database_title = forms.CharField(max_length=200, required=True)
    database_description = forms.CharField(max_length=200, required=True)
    assembly_entries = forms.IntegerField(min_value=1,required=True)
    assembly_accession = forms.CharField(max_length=200,required=False)

    assembly_level = forms.ChoiceField(
        choices=(("Chromosome","Chromosome"),("Complete Genome","Complete Genome"),("Scaffold","Scaffold"),("Contig","Contig")),
        required=True)

    organism_name = forms.CharField(
        label="Organism Name", required=False
    )


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

    user_email = forms.CharField(max_length=200,required=False)

    def __init__(self,user,*args,**kwargs):
        super(UploadGenomeForm,self).__init__(*args,**kwargs)
        self.fields['user_email'].charfield = user.email
        self.fields['user_email'].initial = user.email

    def clean_genome_fasta_file(self):
        genome_fasta_file = self.cleaned_data['genome_fasta_file']
        if genome_fasta_file.name.endswith('.faa') != True and genome_fasta_file.name.endswith('.fasta') != True:
            raise ValidationError('specify a valid fasta file (with .faa or .fasta file name extension)')
        else:
            return genome_fasta_file

    #TODO add controll
    def clean(self):
        cleaned_data = super().clean()
        genome_fasta_file = cleaned_data['genome_fasta_file']
        taxmap_file = cleaned_data['taxmap_file']
        taxonomic_node = cleaned_data['taxonomic_node']
        organism_file = cleaned_data['organism_name_file']
        assembly_accessions_file = cleaned_data['assembly_accessions_file']
        assembly_level_file = cleaned_data['assembly_level_file']
        user_email = cleaned_data['user_email']

        if taxonomic_node != None:
            try:
                check_given_taxonomic_node(user_email,taxonomic_node)
            except Exception as e:
                self.add_error('taxonomic_node','there is no organism with your specified taxonomic node, exception : {}'.format(e))

        if taxmap_file == None and taxonomic_node == None:
           self.add_error('taxonomic_node','specify a taxonomic node for the formatting procedure')
           self.add_error('taxmap_file','specify a taxmap file for the formatting procedure')

        if taxmap_file != None:

            if taxonomic_node != None:
                self.add_error('taxonomic_node', 'just specify one, a taxonomic node or a taxmap file')
                self.add_error('taxmap_file', 'just specify one, a taxonomic node or a taxmap file')

            if organism_file == None:
                self.add_error('organism_name_file','if you upload a taxmap file you should also upload a file that contains the organism names of target genomes, separated by lines')

            protein_ids = 0
            for chunk in genome_fasta_file.chunks():
                protein_ids += chunk.decode().count(">")

            taxmap_ids = 0
            for chunk in taxmap_file.chunks():
                taxmap_ids += chunk.decode().count('\n')

            #check if taxonomic_nodes_exists
            organisms = 0
            for chunk in organism_file.chunks():
                lines = chunk.decode().split("\n")
                for line in lines:
                    if line != '':
                        try:
                            taxids = get_species_taxid_by_name(user_email,line)
                            if taxids != True:
                                raise Exception
                        except:
                            self.add_error('organism_name_file','there is no taxonomic node available for : {}'.format(line))
                organisms += chunk.decode().count('\n')

            if taxmap_ids != protein_ids:
                self.add_error('taxmap_file','the amount of provided acc_ids: {} does not match the provided amount of protein_ids: {}'.format(taxmap_ids,protein_ids))


            if assembly_accessions_file != None:
                amount_of_assemblies = 0
                for chunk in assembly_accessions_file.chunks():
                    lines = chunk.decode().split('\n')
                    #print(lines)
                    for line in lines:
                        if line == '\r':
                            self.add_error('assembly_accessions_file','there are lines without any informations in your assembly accessions file')
                    amount_of_assemblies += chunk.decode().count('\n')
                if amount_of_assemblies != organisms:
                    self.add_error('assembly_accessions_file','the amount of assemblies: {} does not match the amount of provided organisms: {}'.format(amount_of_assemblies,organisms))

            if assembly_level_file != None:
                amount_of_levels = 0
                for chunk in assembly_level_file.chunks():
                    amount_of_levels += chunk.decode().count('\n')
                if amount_of_levels != organisms:
                    self.add_error('assembly_level_file','the amount of assembly levels: {} does not match the amount of provided organisms: {}'.format(amount_of_levels,organisms))

class UploadMultipleFilesGenomeForm(forms.Form):
    database_title = forms.CharField(max_length=200, required=True)
    database_description = forms.CharField(max_length=200, required=True)

    genome_file_field_0 = forms.FileField(required=True)
    organism_name_0 = forms.CharField(max_length=200, required=True)
    extra_field_count = forms.CharField(initial="0",widget=forms.HiddenInput())
    user_email = forms.CharField(max_length=200,required=False)


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

    #https://stackoverflow.com/questions/6142025/dynamically-add-field-to-a-form

    def clean(self):
        cleaned_data = super().clean()
        user_email = cleaned_data['user_email']

        for field in self.fields:
            if "genome_file" in field:
                file = cleaned_data.get(field)
                if file == None:
                    self.add_error(field,'You have to provide a genome file')

                elif file.name.split(".")[-1] in ["fasta","faa","fa"] == False:                    self.add_error(field,'You have to upload a FASTA file, if you provide a valid FASTA file make sure to have a file ending with .faa, .fasta or .fa!')

            elif "organism" in field:
                if cleaned_data.get(field) == '':
                    self.add_error(field,'You have to provide a valid scientific name')
                elif cleaned_data.get(field) == None:
                    self.add_error(field,"You have to provide a valid scientific name")
                else:
                    try:
                        taxids = get_species_taxid_by_name(user_email,cleaned_data.get(field))
                        if taxids != True:
                            raise Exception
                    except Exception as e:
                        self.add_error(field,"{} : {} is no valid name!".format(e,cleaned_data.get(field)))




