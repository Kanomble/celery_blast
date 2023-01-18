from django import forms
from os import listdir
from django.core.exceptions import ValidationError
from blast_project.py_biopython import get_list_of_species_taxids_by_list_of_scientific_names

'''list_taxonomic_files

    utilization in create_taxonomic_file_view and refseqdatabaseform 
    returns a list of all files and their corresponding total line length in the media/taxonomic_node_files folder that end with .taxids

'''
def list_taxonomic_files():
    try:
        files_in_taxonomic_node_files = listdir('media/taxonomic_node_files/')
        files = []
        length = []
        for file in files_in_taxonomic_node_files:
            lines = 0
            with open('media/taxonomic_node_files/'+file) as f:
                for line in f:
                    lines = lines + 1
            if file.endswith('.taxids'):
                files.append(file)
                length.append(lines)
        #[file for file in files_in_taxonomic_node_files if file.endswith('.taxids')]
        return files, length
    except Exception as e:
        raise Exception('[-] exception ocurred in blast_project/py_services.list_taxonomic_files : {}'.format(e))

'''get_taxonomic_files_tuple

    This function is used in the RefseqDatabaseForm for displaying available 
    taxonomic files. It returns a tuple with taxonomic_file_names.

    :returns form_choice_field_input
        :type tuple[list:str, list:str]
'''
def get_taxonomic_files_tuple():
    taxid_files_list = list_taxonomic_files()[0]
    form_choice_field_input = tuple(zip(taxid_files_list, taxid_files_list))
    return form_choice_field_input

'''RefseqDatabaseForm
    
    Input form class for downloading refseq proteoms from the NCBI-FTP server.
    
'''
class RefseqDatabaseForm(forms.Form):
    ASSEMBLY_LEVELS = [
        ('Scaffold','Scaffold'),
        ('Chromosome','Chromosome'),
        ('Contig','Contig'),
        ('Complete Genome','Complete Genome')]

    assembly_levels = forms.MultipleChoiceField(
        label="Assembly Completeness Level",
        required=False,
        choices=ASSEMBLY_LEVELS,
        widget=forms.CheckboxSelectMultiple()
    )

    database_name = forms.CharField(
        label="Database title",
        error_messages={'required': 'A database title is required!'},
        required=True
    )

    database_description = forms.CharField(
        label="Database description",
        error_messages={'required' : 'Add a short description about the purpose of this database!'}
    )

    taxid_file = forms.FileField(
        label="Optional: File for limiting refseq databases by taxonomy",
        required=False
    )

    taxid_text_field = forms.CharField(
        required=False, max_length=800
    )

    taxid_uploaded_file = forms.ChoiceField(
        choices=get_taxonomic_files_tuple(),
        required=False
    )

    user_email = forms.CharField(max_length=200)


    def __init__(self, user, *args, **kwargs):
        super(RefseqDatabaseForm, self).__init__(*args, **kwargs)
        self.fields['taxid_uploaded_file'].choices = get_taxonomic_files_tuple()
        self.fields['user_email'].charfield = user.email
        self.fields['user_email'].initial = user.email

    def clean(self):
        try:
            cleaned_data = super().clean()
            taxid_file = cleaned_data['taxid_file']
            taxid_uploaded_file = cleaned_data['taxid_uploaded_file']
            taxid_text_field = cleaned_data['taxid_text_field']
            user_email = self.fields['user_email'].charfield

            #if taxid_file empty --> taxid_file == None
            #if taxid_uploaded_file empty --> taxid_uploaded_file == ''
            #if taxid_text_field empty --> taxid_text_field == ''
            if taxid_file != None and taxid_uploaded_file != '' and taxid_text_field != '':
               self.add_error('taxid_file','Choose just one taxonomic limitation: file upload, available taxonomic files or provide valid scientific names!')

            if taxid_file != None and taxid_uploaded_file != '':
               self.add_error('taxid_file','Choose just one taxonomic limitation: file upload, available taxonomic files or provide valid scientific names!')

            if taxid_file != None and taxid_text_field != '':
               self.add_error('taxid_file','Choose just one taxonomic limitation: file upload, available taxonomic files or provide valid scientific names!')

            if taxid_uploaded_file != '' and taxid_text_field != '':
               self.add_error('taxid_file','Choose just one taxonomic limitation: file upload, available taxonomic files or provide valid scientific names!')

            if taxid_file != None:
                for line in list_taxonomic_files()[0]:
                    if taxid_file.name == line:
                        self.add_error('taxid_file',"This taxid_file already exist!")

                if taxid_file.name.endswith('.taxid') != True and taxid_file.name.endswith('.taxids') != True:
                    self.add_error('taxid_file', "Are you sure that this file is a taxonomic nodes file?"
                                                 "\n Your file should only contain one taxonomic nodes for each line "
                                                 "AND should be named to {species}.taxid or {species}.taxids"
                                   )

            elif taxid_text_field != '':
                species_names = taxid_text_field.split(",")

                try:
                    taxids, errors = get_list_of_species_taxids_by_list_of_scientific_names(user_email, species_names)
                    if len(errors) > 0:
                        error_string = ''
                        for name in errors:
                            error_string = error_string + name + ' '
                            self.add_error("taxid_text_field", "{}".format(name))

                    cleaned_data['taxid_text_field'] = taxids
                except Exception as e:
                    self.add_error("taxid_text_field","{}".format(e))

        except Exception as e:
            raise ValidationError(
                "validation error in refseq database creation, due to this exception: {}".format(
                     e))
        return cleaned_data


