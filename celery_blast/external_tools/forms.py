from django import forms

'''
EntrezSearchForm for the EntrezSearch model, which instances hold information 
of conducted entrez searches. Supported databases can be included into the DATABASES list.
'''
class EntrezSearchForm(forms.Form):
    # for adding  more databases, the NCBI database need to be added here, in the models.py get_pandas_table function and in entrez_search_service.py to the execute_entrez_search function
    DATABASES = [('pubmed','Pubmed'),
                 ('protein','Protein'),
                 ('assembly','Assembly'),
                 ('cdd','CDD'),
                 ('protfam','Protfam')]

    entrez_query = forms.CharField(initial="ribosomes AND cyanobacteria AND review [PT]", widget=forms.widgets.Textarea)

    database = forms.ChoiceField(choices=DATABASES,
                                 initial=('protein','protein'))