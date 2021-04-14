from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User
from django.core.exceptions import ValidationError

from .py_biopython import get_species_taxid_by_name

''' CreateTaxonomicFileForm
post form for the create_taxonomic_file.html template

__init__ :
    directly sets the user_email field with the request.user.email attribute
    this is used as input for the get_species_taxid_by_name function

clean_species_name : 
    translates the scientific name into a taxonomic node that can be used in the view
    for triggering the get_species_taxids.sh script (blast_project/tasks.get_species_taxids_into_file)
'''
class CreateTaxonomicFileForm(forms.Form):
    species_name = forms.CharField(max_length=200, required=True,)
    user_email = forms.CharField(max_length=200)
    taxonomic_node = forms.IntegerField(required=False)

    def __init__(self, user, *args, **kwargs):
        super(CreateTaxonomicFileForm, self).__init__(*args, **kwargs)
        self.fields['user_email'].charfield = user.email
        self.fields['user_email'].initial = user.email
        #can be turned on, this would hide the relevant field div in the template
        #self.fields['user_email'].widget = forms.HiddenInput()

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