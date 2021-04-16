from django import forms
from django.core.exceptions import ValidationError

class RefseqTableForm(forms.Form):
    ASSEMBLY_LEVELS = [
        ('Scaffold','Scaffold'),
        ('Chromosome','Chromosome'),
        ('Contig','Contig'),
        ('Complete Genome','Complete Genome')]

    assembly_levels = forms.MultipleChoiceField(
        required=False,
        choices=ASSEMBLY_LEVELS)
