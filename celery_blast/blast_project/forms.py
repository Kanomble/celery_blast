from django import forms
from django.core.exceptions import ValidationError
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User

from .py_biopython import get_species_taxid_by_name

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