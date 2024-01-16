from string import punctuation, ascii_letters

from blast_project.py_biopython import fetch_protein_records
from blast_project.py_django_db_services import get_all_succeeded_databases
from django import forms


# TODO documentation
class OneWayProjectCreationForm(forms.Form):
    class BlastDatabaseModelChoiceField(forms.ModelChoiceField):
        def label_from_instance(self, blast_database):
            return str(blast_database.database_name)

    project_title = forms.CharField(
        label="Project title",
        error_messages={
            'required': "A project title is required for saving project metadata into the database"})

    query_sequence_file = forms.FileField(
        required=False,
        error_messages={
            'required': "Upload a query sequence file, this file will serve as the -query parameter for the forward BLAST analysis"})

    query_sequence_text = forms.CharField(
        label="Query Sequence IDs", max_length=800, required=False
    )

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

    def clean(self):
        cleaned_data = super().clean()
        query_file = cleaned_data['query_sequence_file']
        query_sequences = cleaned_data['query_sequence_text']
        user_email = cleaned_data['user_email']

        # upload a query file or specify valid protein identifiers
        if query_file is None and query_sequences == '':
            self.add_error('query_sequence_file',
                           "please upload a fasta file containing your sequences or specify valid protein identifier")

        # query file was uploaded
        if query_file is not None:
            if query_file.name.endswith('.faa') != True and query_file.name.endswith('.fasta') != True:
                raise self.add_error('query_sequence_file', "please upload only fasta files!")

            if len(query_file.name.split(".")) != 2:
                raise self.add_error('query_sequence_file', "there are no dots allowed except the filetype delimiter")
            else:
                filename = query_file.name.split(".")[0]
                for character in filename:
                    if character in punctuation:
                        if character != '_' and character != '-':
                            raise self.add_error('query_sequence_file',
                                                 "bad character: \"{}\" in query file name".format(character))
                    if character not in ascii_letters:
                        if character != '_' and character != '-':
                            raise self.add_error('query_sequence_file',
                                                 "bad character: \"{}\" in query file name".format(character))

            header = []
            for chunk in query_file.chunks():
                lines = chunk.decode().split("\n")
                for line in lines:
                    if line.startswith('>'):
                        try:
                            acc = line.split(" ")[0].split('>')[-1].split(".")[0]
                            header.append(acc)
                        except Exception as e:
                            self.add_error('query_sequence_file', 'error during parsing of query_file : {}'.format(e))

            if len(header) > 300:
                self.add_error('query_sequence_file', 'You try to infer orthologs for more than 300 query sequences,'
                                                      ' this is not allowed, consider to separate the query sequences.')

            if len(header) != len(set(header)):
                self.add_error('query_sequence_file',
                               'there are duplicate proteins in your uploaded file, please remove the duplicate entries and upload the file again!')
        # protein identifier have been uploaded
        elif query_sequences != '' and query_file is None:
            # check string for invalid characters
            query_sequences = query_sequences.replace(" ", "").split(',')
            query_sequences = [qseq.split(".")[0] for qseq in query_sequences]
            try:
                proteins, errors = fetch_protein_records(query_sequences, user_email)
                if len(errors) > 0:
                    self.add_error('query_sequence_text', 'following sequences are unavailable: {}'.format(errors))
                cleaned_data['query_sequence_text'] = proteins
            except Exception as e:
                self.add_error("query_sequence_text", "please provide valid protein identifiers")
        else:
            self.add_error('query_sequence_text',
                           "please upload a fasta file containing your sequences or specify valid protein identifier")

        return cleaned_data


# TODO documentation
class BlastSettingsForm(forms.Form):
    e_value = forms.FloatField(
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


# TODO documentation
class OneWayRemoteProjectCreationForm(forms.Form):
    BLAST_SEARCH_PROGRAMS = [('blastp', 'search against protein databases'),
                             ]
    BLAST_REMOTE_DATABASES = [
        ('nr', 'none redundant proteins'),
        ('refseq_protein', 'refseq protein database'),
        ('refseq_select_prot', 'refseq selected protein database'),
        ('swissprot', 'protein sequences from the swissprot database')]
    #        ('refseq_rna', 'refseq rna database'), ('blastn', 'search against nucleotide databases')
    #         ('refseq_select_rna', 'refseq selected rna database'),
    #        ('nt', 'none redundant dna'),
    #    ('env_nr', 'nr for not yet known organisms (env_nr)'),
    #    ('env_nt', 'nt for not yet known organisms (env_nt)'),

    r_project_title = forms.CharField(
        label="Project title",
        error_messages={
            'required': "A project title is required for saving project metadata into the database"})

    r_query_sequence_file = forms.FileField(
        required=False,
        error_messages={
            'required': "Upload a query sequence file, this file will serve as the -query parameter for the BLAST analysis"})

    r_query_sequence_text = forms.CharField(
        label="Query Sequence IDs",
        max_length=800,
        required=False
    )

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
    )

    def __init__(self, user, *args, **kwargs):
        super(OneWayRemoteProjectCreationForm, self).__init__(*args, **kwargs)
        self.fields['r_user_email'].charfield = user.email
        self.fields['r_user_email'].initial = user.email

    def clean(self):
        cleaned_data = super().clean()
        query_file = cleaned_data['r_query_sequence_file']
        query_sequences = cleaned_data['r_query_sequence_text']
        user_email = cleaned_data['r_user_email']

        # upload a query file or specify valid protein identifiers
        if query_file == None and query_sequences == '':
            self.add_error('r_query_sequence_file',
                           "please upload a fasta file containing your sequences or specify valid protein identifier")

        # query file was uploaded
        if query_file != None:
            if query_file.name.endswith('.faa') != True and query_file.name.endswith('.fasta') != True:
                raise self.add_error('r_query_sequence_file', "please upload only fasta files!")

            if len(query_file.name.split(".")) != 2:
                raise self.add_error('r_query_sequence_file', "there are no dots allowed except the filetype delimiter")
            else:
                filename = query_file.name.split(".")[0]
                for character in filename:
                    if character in punctuation:
                        if character != '_' and character != '-':
                            raise self.add_error('r_query_sequence_file',
                                                 "bad character: \"{}\" in query file name".format(character))
                    if character not in ascii_letters:
                        if character != '_' and character != '-':
                            raise self.add_error('r_query_sequence_file',
                                                 "bad character: \"{}\" in query file name".format(character))

            header = []
            for chunk in query_file.chunks():
                lines = chunk.decode().split("\n")
                for line in lines:
                    if line.startswith('>'):
                        try:
                            acc = line.split(" ")[0].split('>')[-1].split(".")[0]
                            header.append(acc)
                        except Exception as e:
                            self.add_error('r_query_sequence_file', 'error during parsing of query_file : {}'.format(e))

            if len(header) > 300:
                self.add_error('r_query_sequence_file', 'You try to infer orthologs for more than 300 query sequences,'
                                                        ' this is not allowed, consider to separate the query sequences.')

            if len(header) != len(set(header)):
                self.add_error('r_query_sequence_file',
                               'there are duplicate proteins in your uploaded file, please remove the duplicate entries and upload the file again!')
        # protein identifier have been uploaded
        elif query_sequences != '' and query_file is None:
            # check string for invalid characters
            query_sequences = query_sequences.replace(" ", "").split(',')
            query_sequences = [qseq.split(".")[0] for qseq in query_sequences]
            try:
                proteins, errors = fetch_protein_records(query_sequences, user_email)
                if len(errors) > 0:
                    self.add_error('r_query_sequence_text', 'following sequences are unavailable: {}'.format(errors))
                cleaned_data['r_query_sequence_text'] = proteins
            except Exception as e:
                self.add_error("r_query_sequence_text", "please provide valid protein identifiers: exception: {}".format(e))

            # self.add_error('r_query_sequence_text','not available yet')
        else:
            self.add_error('r_query_sequence_text',
                           "please upload a fasta file containing your sequences or specify valid protein identifier")

        return cleaned_data
