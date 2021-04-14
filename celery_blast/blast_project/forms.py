from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User

#create_taxonomic_file.html form
class CreateTaxonomicFileForm(forms.Form):
    species_name = forms.CharField(max_length=200, required=True,)
    user_email = forms.CharField(max_length=200)

    def __init__(self, user, *args, **kwargs):
        super(CreateTaxonomicFileForm, self).__init__(*args, **kwargs)
        self.fields['user_email'].charfield = user.email
        self.fields['user_email'].initial = user.email
        #self.fields['user_email'].widget = forms.HiddenInput()

#registration and login form
class CreateUserForm(UserCreationForm):
    class Meta:
        model = User
        fields = ['username','email','password1','password2']