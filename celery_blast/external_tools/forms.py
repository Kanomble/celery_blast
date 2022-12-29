from django import forms
from .models import ExternalTools

'''
    EntrezSearchForm for the EntrezSearch model, which instances hold information 
    of conducted entrez searches. Supported databases can be included into the DATABASES list.
'''
class EntrezSearchForm(forms.Form):

    DATABASES = [('pubmed','pubmed'),
                 ('protein','protein')]

    entrez_query = forms.CharField(initial="ribosomes AND cyanobacteria AND review [PT]", widget=forms.widgets.Textarea)

    database = forms.ChoiceField(choices=DATABASES,
                                 initial=('protein','protein'))


'''RpsBLASTSettingsForm

    Form for RPSBLAST settings.
    
'''
class RpsBLASTSettingsForm(forms.Form):

    query_sequence = forms.ChoiceField(
    )

    rps_e_value = forms.FloatField(
        label="RPSBLAST E-Value",
        initial=0.001, min_value=0.0000001, max_value=10)

    rps_num_alignments = forms.IntegerField(
        label="RPSBLAST Number of possible alignment outputs",
        initial=10000, min_value=1, max_value=50000)

    rps_max_target_seqs = forms.IntegerField(
        label="RPSBLAST max_target_seqs of possible alignment description outputs",
        initial=10000,min_value=1, max_value=50000)

    rps_num_threads = forms.IntegerField(
        label="RPSBLAST Threads",
        initial=1, min_value=1, max_value=16)

    rps_max_hsps = forms.IntegerField(
        label='RPSBLAST max hsps',
        initial=500, min_value=1, max_value=1000)

    def __init__(self, query_sequences_rdy_for_cdd, *args, **kwargs):
        super(RpsBLASTSettingsForm, self).__init__(*args, **kwargs)
        self.fields['query_sequence'] = forms.ChoiceField(
            choices=tuple(
                [(query_name.query_accession_id,query_name.query_accession_id) for query_name in query_sequences_rdy_for_cdd]))