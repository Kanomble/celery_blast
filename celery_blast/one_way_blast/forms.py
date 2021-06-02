from django import forms
from django.core.exceptions import ValidationError
from blast_project.py_django_db_services import get_all_succeeded_databases

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
        super(ProjectCreationForm, self).__init__(*args, **kwargs)
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