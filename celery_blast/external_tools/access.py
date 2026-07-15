from django.db.models import Q
from django.shortcuts import get_object_or_404

from .models import EntrezSearch, QuerySequences


def get_owned_entrez_search_or_404(user, search_id):
    return get_object_or_404(EntrezSearch, id=search_id, entrez_user=user)


def get_owned_query_sequence_or_404(user, query_sequence_id):
    owned_query_sequences = QuerySequences.objects.select_related(
        'external_tool_for_query_sequence__associated_project__project_user',
        'external_tool_for_query_sequence__remote_associated_project__r_project_user',
        'multiple_sequence_alignment_task',
        'phylogenetic_tree_construction_task',
    ).filter(
        Q(external_tool_for_query_sequence__associated_project__project_user=user)
        | Q(external_tool_for_query_sequence__remote_associated_project__r_project_user=user)
    )
    return get_object_or_404(
        owned_query_sequences,
        id=query_sequence_id,
    )
