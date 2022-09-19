from django import forms

'''
EntrezSearchForm for the EntrezSearch model, which instances hold information 
of conducted entrez searches. Supported databases can be included into the DATABASES list.
'''
class EntrezSearchForm(forms.Form):

    DATABASES = [('pubmed','pubmed'),
                 ('protein','protein'),
                 ('assembly','assembly'),
                 ('cdd','cdd'),
                 ('protfam','protfam')]

    entrez_query = forms.CharField(initial="ribosomes AND cyanobacteria AND review [PT]", widget=forms.widgets.Textarea)

    database = forms.ChoiceField(choices=DATABASES,
                                 initial=('protein','protein'))