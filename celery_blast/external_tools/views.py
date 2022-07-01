from django.shortcuts import render, redirect, HttpResponse
from django.contrib.auth.decorators import login_required
from blast_project.views import failure_view, success_view
from .tasks import execute_multiple_sequence_alignment, execute_phylogenetic_tree_building,\
    execute_multiple_sequence_alignment_for_all_query_sequences, execute_fasttree_phylobuild_for_all_query_sequences,\
    entrez_search_task, download_entrez_search_associated_protein_sequences
from .models import ExternalTools, EntrezSearch, QuerySequences
from .forms import EntrezSearchForm
from .entrez_search_service import get_entrezsearch_object_with_entrezsearch_id, delete_esearch_by_id
import os
from django.http import JsonResponse
from time import sleep

@login_required(login_url='login')
def view_downloaded_sequences(request, search_id):
    try:
        entrez_search = get_entrezsearch_object_with_entrezsearch_id(search_id)
        filepath = entrez_search.fasta_file_name
        if os.path.isfile(filepath):
            with open(filepath,'r') as download_file:
                content = [str(line) for line in download_file.readlines()]
                content = ''.join(content)
            response = HttpResponse(content,content_type="text/plain")
            response['Contnt-Disposition'] = "attachment; filename={}".format(filepath.split("/")[-1])
        else:
            raise FileNotFoundError

        return response
    except Exception as e:
        return failure_view(request,e)


@login_required(login_url='login')
def search_detail_view(request, search_id):
    #get edirectpaper entry based on search_id (which is id of db object)
    #get paper content and fill context object with edirectpaper and paper content
    #return template with context
    try:
        sleep(0.5)
        entrez_search = get_entrezsearch_object_with_entrezsearch_id(search_id)
        context = {'EntrezSearch':entrez_search,'HtmlTable':entrez_search.get_paper_content()}

        if entrez_search.database == 'protein' and entrez_search.download_task_result == None:
            context['DownloadProteins'] = True
        elif entrez_search.database == 'protein' and entrez_search.download_task_result != None:
            if entrez_search.download_task_result.status == 'SUCCESS':
                context['DownloadTaskSuccess'] = 0
            elif entrez_search.download_task_result.status == 'PROGRESS':
                context['DownloadTaskSuccess'] = 1
            else:
                context['DownloadTaskSuccess'] = 2
        return render(request,'external_tools/search_details.html',context)
    except Exception as e:
        return failure_view(request,e)

@login_required(login_url='login')
def download_proteins_from_entrez_search(request, search_id):
    try:
        download_entrez_search_associated_protein_sequences.delay(search_id)
        return redirect('search_details',search_id=search_id)
    except Exception as e:
        return failure_view(request,e)

@login_required(login_url='login')
def delete_search_view(request, search_id):
    try:
        retcode = delete_esearch_by_id(search_id)
        if retcode == 0:
            #entrez_search_form = EntrezSearchForm()
            #context = {"EntrezSearches":  EntrezSearch.edirect_objects.get_all_entrez_searches_from_current_user(
            #        request.user.id),
            #           "EntrezSearchForm":entrez_search_form}
            return redirect("entrez_dashboard")
        elif retcode == 1:
            return failure_view(request, "Couldnt delete esearch object")
        else:
            return failure_view(request, "ERROR during deletion of edirectpaper with id: {}".format(search_id))

    except Exception as e:
        return failure_view(request,e)

@login_required(login_url='login')
def entrez_dashboard_view(request):
    try:
        context = {}
        if request.method == "POST":
            entrez_search_form = EntrezSearchForm(request.POST)
            if entrez_search_form.is_valid():
                database = entrez_search_form.cleaned_data['database']
                entrez_query = entrez_search_form.cleaned_data['entrez_query']
                entrez_search_task.delay(database,entrez_query, request.user.id)

                entrez_searches = EntrezSearch.edirect_objects.get_all_entrez_searches_from_current_user(
                    request.user.id)
                context['EntrezSearches'] = entrez_searches
                context['EntrezSearchForm'] = entrez_search_form

                return redirect('entrez_dashboard')

        else:
            sleep(0.5)
            entrez_searches = EntrezSearch.edirect_objects.get_all_entrez_searches_from_current_user(
                request.user.id)
            context['EntrezSearches'] = entrez_searches
            entrez_search_form = EntrezSearchForm()
            context['EntrezSearchForm'] = entrez_search_form

        return render(request,'external_tools/entrez_search_dashboard.html',context)
    except Exception as e:
        return failure_view(request,e)


#TODO documentation
@login_required(login_url='login')
def project_informations(request, project_id):
    try:
        context = {}
        qseqids = ExternalTools.objects.get_external_tools_based_on_project_id(project_id)
        context['qseqids'] = qseqids
        context['project_id'] = project_id
        return render(request,"external_tools/external_tools_dashboard.html",context)
    except Exception as e:
        return failure_view(request,e)

#TODO documentation
@login_required(login_url='login')
def perform_simple_msa(request,project_id,query_sequence_id):
    try:
        if request.method == 'POST':
            execute_multiple_sequence_alignment.delay(project_id,query_sequence_id)
            return redirect('external_project_informations',project_id=project_id)
        else:
            e = "There is no GET method for this view function"
            return failure_view(request,e)
    except Exception as e:
        return failure_view(request,e)

@login_required(login_url='login')
def perform_simple_msa_for_all_query_sequences(request,project_id):
    try:
        if request.method == 'POST':
            execute_multiple_sequence_alignment_for_all_query_sequences.delay(project_id)
            return redirect('external_project_informations',project_id=project_id)
        else:
            e = "There is no GET method for this view function"
            return failure_view(request,e)
    except Exception as e:
        return failure_view(request,e)

@login_required(login_url='login')
def perform_fasttree_phylobuild_for_all_query_sequences(request,project_id):
    try:
        if request.method == 'POST':
            execute_fasttree_phylobuild_for_all_query_sequences.delay(project_id)
            return redirect('external_project_informations',project_id=project_id)
        else:
            e = "There is no GET method for this view function"
            return failure_view(request,e)
    except Exception as e:
        return failure_view(request,e)

#TODO documentation
@login_required(login_url='login')
def perform_fasttree_phylobuild(request,project_id,query_sequence_id):
    try:
        if request.method == 'POST':
            execute_phylogenetic_tree_building.delay(project_id,query_sequence_id)
            return redirect('external_project_informations',project_id=project_id)
        else:
            e = "There is no GET method for this view function"
            return failure_view(request,e)
    except Exception as e:
        return failure_view(request,e)

@login_required(login_url='login')
def ajax_call_progress_phylo_task(request, query_sequence_id):
    try:
        if request.is_ajax and request.method == "GET":
            qseq = QuerySequences.objects.get(id=query_sequence_id)
            progress = qseq.phylogenetic_tree_construction_task.status
            return JsonResponse({"progress":progress})
        else:
            return JsonResponse({"progress":"error"})
    except Exception as e:
        return failure_view(request,e)

@login_required(login_url='login')
def ajax_call_progress_msa_task(request, query_sequence_id):
    try:
        if request.is_ajax and request.method == "GET":
            qseq = QuerySequences.objects.get(id=query_sequence_id)
            progress = qseq.multiple_sequence_alignment_task.status
            return JsonResponse({"progress":progress})
        else:
            return JsonResponse({"progress":"error"})
    except Exception as e:
        return failure_view(request,e)

@login_required(login_url='login')
def ajax_call_progress_entrezsearch_to_fasta(request,search_id):
    try:
        if request.is_ajax and request.method == "GET":
            entrez_search = get_entrezsearch_object_with_entrezsearch_id(search_id)
            progress = entrez_search.download_task_result.status
            return JsonResponse({"progress":progress})
        else:
            return JsonResponse({"progress":"error"})
    except Exception as e:
        return failure_view(request,e)

#view for phylogenetic dashboard
@login_required(login_url='login')
def phylogenetic_information(request, project_id, query_sequence_id):
    try:
        context = {}
        qseqids = ExternalTools.objects.get_external_tools_based_on_project_id(project_id)
        context['qseqids'] = qseqids
        context['query_sequence_id']= query_sequence_id
        context['project_id'] = project_id
        return render(request,"external_tools/phylogenetic_dashboard.html",context)
    except Exception as e:
        return failure_view(request,e)
