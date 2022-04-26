from django import forms

class EntrezSearchForm(forms.Form):

    DATABASES = [('pubmed','pubmed'),
                 ('protein','protein')]

    entrez_query = forms.CharField(max_length=500,
                                   required=True)

    database = forms.ChoiceField(choices=DATABASES,
                                 initial=('protein','protein'))