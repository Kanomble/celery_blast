from django import forms
from blast_project.py_services import list_taxonomic_files, get_taxonomic_files_tuple

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

    taxid_uploaded_file = forms.ChoiceField(
        choices=get_taxonomic_files_tuple(),
        required=False
    )


    def __init__(self, data=None, *args, **kwargs):
        super(RefseqDatabaseForm, self).__init__(data, *args, **kwargs)
        self.fields['taxid_uploaded_file'].choices = get_taxonomic_files_tuple()

    def clean(self):
        cleaned_data = super().clean()
        taxid_file = cleaned_data['taxid_file']
        taxid_uploaded_file = cleaned_data['taxid_uploaded_file']

        if taxid_file != None and taxid_uploaded_file != '':
           self.add_error('taxid_file','Just use a previously uploaded file if you dont specify a file to upload!')

        if taxid_file != None:
            for line in list_taxonomic_files()[0]:
                if taxid_file.name == line:
                    self.add_error('taxid_file',"This taxid_file already exist!")

            if taxid_file.name.endswith('.taxid') != True and taxid_file.name.endswith('.taxids') != True:
                self.add_error('taxid_file', "Are you sure that this file is a taxonomic nodes file?"
                                             "\n Your file should only contain one taxonomic nodes for each line "
                                             "AND should be named to {species}.taxid or {species}.taxids"
                               )

