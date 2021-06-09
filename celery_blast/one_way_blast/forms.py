from django import forms
from django.core.exceptions import ValidationError
from blast_project.py_django_db_services import get_all_succeeded_databases
from refseq_transactions.forms import get_taxonomic_files_tuple
from blast_project.py_services import list_taxonomic_files

#TODO documentation
class OneWayProjectCreationForm(forms.Form):
    class BlastDatabaseModelChoiceField(forms.ModelChoiceField):
        def label_from_instance(self, blast_database):
            return str(blast_database.database_name)

    project_title = forms.CharField(
        label="Project title",
        error_messages={
            'required': "A project title is required for saving project metadata into the database"})

    query_sequence_file = forms.FileField(
        error_messages={
            'required':"Upload a query sequence file, this file will serve as the -query parameter for the forward BLAST analysis"})

    project_database = BlastDatabaseModelChoiceField(
        queryset=get_all_succeeded_databases(),
        empty_label=None
    )

    user_email = forms.CharField(
        max_length=200,
        required=False)

    def __init__(self, user, *args, **kwargs):
        super(OneWayProjectCreationForm, self).__init__(*args, **kwargs)
        self.fields['user_email'].charfield = user.email
        self.fields['user_email'].initial = user.email

    def clean_query_sequence_file(self):
        query_file = self.cleaned_data['query_sequence_file']
        if query_file.name.endswith('.faa') != True and query_file.name.endswith('.fasta') != True:
            raise ValidationError("please upload only fasta files!")
        else:
            return query_file

#TODO documentation
class BlastSettingsForm(forms.Form):
    e_value = forms.DecimalField(
        label="E-Value", initial=0.001)
    word_size = forms.IntegerField(
        label="Word Size", initial=3)
    num_alignments = forms.IntegerField(
        label="Number of possible alignment outputs", initial=10000)
    num_threads = forms.IntegerField(
        label="Threads", initial=1)
    max_hsps = forms.IntegerField(
        label='Max HSPS', initial=500
    )


#TODO documentation
class OneWayRemoteProjectCreationForm(forms.Form):

    BLAST_SEARCH_PROGRAMS = [('blastp', 'search against protein databases'), ('blastn', 'search against nucleotide databases')]
    BLAST_REMOTE_DATABASES = [( 'nr','none redundant proteins'), ('nt','none redundant dna')]

    r_project_title = forms.CharField(
        label="Project title",
        error_messages={
            'required': "A project title is required for saving project metadata into the database"})

    r_query_sequence_file = forms.FileField(
        required=True,
        error_messages={
            'required':"Upload a query sequence file, this file will serve as the -query parameter for the forward BLAST analysis"})


    r_project_database = forms.ChoiceField(
        choices=BLAST_REMOTE_DATABASES
    )

    r_search_strategy = forms.ChoiceField(
        choices=BLAST_SEARCH_PROGRAMS
    )

    r_user_email = forms.CharField(
        max_length=200,
        required=False)

    r_entrez_query = forms.CharField(
        max_length=200,
        required=False,
        error_messages = {
            'required': "Upload a query sequence file, this file will serve as the -query parameter for the forward BLAST analysis"}
    )

    def __init__(self, user, *args, **kwargs):
        super(OneWayRemoteProjectCreationForm, self).__init__(*args, **kwargs)
        self.fields['r_user_email'].charfield = user.email
        self.fields['r_user_email'].initial = user.email

    def clean_r_query_sequence_file(self):
        query_file = self.cleaned_data['r_query_sequence_file']
        if query_file.name.endswith('.faa') != True and query_file.name.endswith('.fasta') != True:
            raise ValidationError("please upload only fasta files!")
        else:
            return query_file

