from django import forms
from django.core.exceptions import ValidationError
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User
from django.forms import ModelChoiceField

from .py_biopython import get_species_taxid_by_name
from .py_django_db_services import get_all_succeeded_databases
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

#TODO documentation
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

    project_database = BlastDatabaseModelChoiceField(
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

    def clean_species_name_for_backward_blast(self):
        species_name = self.cleaned_data['species_name_for_backward_blast']
        user_email = self.fields['user_email'].charfield
        try:
            taxonomic_node = get_species_taxid_by_name(user_email, species_name)
            return species_name, taxonomic_node
        except Exception as e:
            raise ValidationError(
                "validation error in clean_species_name pls check your provided scientific name : {}".format(e))

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
