from django import forms
from django.core.exceptions import ValidationError
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User
from django.forms import ModelChoiceField

from pandas import read_csv
from .py_biopython import get_species_taxid_by_name, check_given_taxonomic_node
from .py_django_db_services import get_all_succeeded_databases, get_database_by_id, check_if_taxid_is_in_database
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
            taxonomic_node = get_species_taxid_by_name(user_email,species_name)
            return species_name, taxonomic_node
        except Exception as e:
            raise ValidationError("validation error in clean_species_name pls check your provided scientific name : {}".format(e))

#registration and login form
class CreateUserForm(UserCreationForm):
    class Meta:
        model = User
        fields = ['username','email','password1','password2']

''' ProjectCreationForm
    This form is used in the form_project_creation.html.
    All form fields are processed into the BlastProject model instance via the project_creation_view functions 
    in the views.py file. 
    
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

    species_name_for_backward_blast = forms.CharField(
        required=True,
        label='Scientific Names (conversion to Taxonomic Nodes) for Backward BLAST',
        error_messages={
            'required':"Specify a Scientific Name for your backward BLAST - use a comma separated list - names will be converted to taxids that will be written to a file which will serve as the -taxidlist parameter of your backward BLAST"})

    user_email = forms.CharField(
        max_length=200,
        required=False)

    def __init__(self, user, *args, **kwargs):
        super(ProjectCreationForm, self).__init__(*args, **kwargs)
        self.fields['user_email'].charfield = user.email
        self.fields['user_email'].initial = user.email

    #TODO check if backward organism reside in the backward blast database
    def clean_species_name_for_backward_blast(self):
        species_name = self.cleaned_data['species_name_for_backward_blast']
        user_email = self.fields['user_email'].charfield
        try:
            taxonomic_node = get_species_taxid_by_name(user_email, species_name)
            backward_db = self.cleaned_data['project_backward_database']
            forward_db = self.cleaned_data['project_forward_database']
            booleanbw = check_if_taxid_is_in_database(backward_db.id,taxonomic_node)
            booleanfw = check_if_taxid_is_in_database(forward_db.id,taxonomic_node)
            if booleanbw == False:
                self.add_error('species_name_for_backward_blast','specified taxonomic node: {} does not reside in the selected BACKWARD database: {}'.format(taxonomic_node,backward_db.database_name))
            if booleanfw == False:
                self.add_error('species_name_for_backward_blast','specified taxonomic node: {} does not reside in the selected FORWARD database: {}'.format(taxonomic_node,forward_db.database_name))
            if booleanfw == True and booleanbw == True:
                return species_name, taxonomic_node
        except Exception as e:
            raise ValidationError(
                "validation error in clean_species_name pls check your provided scientific name : {}, due to this exception: ".format(species_name,e))

    def clean_query_sequence_file(self):
        query_file = self.cleaned_data['query_sequence_file']
        if query_file.name.endswith('.faa') != True and query_file.name.endswith('.fasta') != True:
            raise ValidationError("please upload only fasta files!")
        else:
            return query_file

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
                            get_species_taxid_by_name(user_email,line)
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
