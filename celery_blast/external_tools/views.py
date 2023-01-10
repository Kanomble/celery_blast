import os
from json import loads
from time import sleep

from blast_project.py_services import get_html_results
from blast_project.views import failure_view
from django.contrib.auth.decorators import login_required
from django.core.handlers.wsgi import WSGIRequest
from django.http import JsonResponse
from django.shortcuts import render, redirect, HttpResponse
from django_celery_results.models import TaskResult

from .entrez_search_service import get_entrezsearch_object_with_entrezsearch_id, delete_esearch_by_id
from .forms import EntrezSearchForm, RpsBLASTSettingsForm
from .models import ExternalTools, EntrezSearch, QuerySequences
from .py_services import delete_cdd_search_output, check_if_cdd_search_can_get_executed
from .tasks import execute_multiple_sequence_alignment, execute_phylogenetic_tree_building, \
    execute_multiple_sequence_alignment_for_all_query_sequences, execute_fasttree_phylobuild_for_all_query_sequences, \
    download_organism_protein_sequences_task, \
    entrez_search_task, download_entrez_search_associated_protein_sequences, cdd_domain_search_with_rbhs_task


@login_required(login_url='login')
def view_downloaded_sequences(request: WSGIRequest, search_id: int):
    # creates a file download button to the search_details.html to download the result file of the download_proteins_from_entrez_search view
    # if the entrezsearch is of the protein or pubmed database
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
def view_downloaded_organism_sequences(request: WSGIRequest, search_id: int,organism_download: str):
    #creates a file download button to the search_details.html to download the result file of the download_protein_by_organisms_from_entrez_search view
    try:
        filepath = "media/esearch_output/"+ str(search_id)+"/"+str(organism_download)+".faa"
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
def search_detail_view(request: WSGIRequest, search_id: int):
    #get edirectpaper entry based on search_id (which is id of db object)
    #get paper content and fill context object with edirectpaper and paper content
    #return template with context
    try:
        sleep(0.5)
        entrez_search = get_entrezsearch_object_with_entrezsearch_id(search_id)
        context = {'EntrezSearch':entrez_search,'HtmlTable':entrez_search.get_paper_content()}
        #context['Stat_col_length'] = entrez_search.get_stat_columns_lenght()

        if entrez_search.database == 'protein':

            organism_progress_dic = {}
            email = request.user.email
            organism_download = entrez_search.get_organisms()
            for organism in organism_download:
                task_arg_str = '"('+str(search_id)+", '"+organism+"', "+ "'"+ email +"'"+')"'
                file_name = "media/esearch_output/"+str(search_id)+"/"+organism+".faa"
                if os.path.isfile(file_name):
                    task_row = TaskResult.objects.get(task_args=task_arg_str)
                    task_status = task_row.status
                    organism_progress_dic[organism] = task_status
                else:
                    organism_progress_dic[organism] = ""
            context['Organism_progress_dic'] = organism_progress_dic
            context['Organisms'] = entrez_search.get_organisms()

        if ((entrez_search.database == 'protein' or entrez_search.database == "pubmed") and entrez_search.download_task_result == None):
            context['DownloadProteins'] = True
        elif ((entrez_search.database == 'protein' or  entrez_search.database == "pubmed") and entrez_search.download_task_result != None):
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
def download_proteins_from_entrez_search(request: WSGIRequest, search_id: int):
    #starts a task that downloads the associated protein sequences in the  entrez search results page
    # if the entrezsearch is of the protein or pubmed database
    try:
        download_entrez_search_associated_protein_sequences.delay(search_id)
        return redirect('search_details',search_id=search_id)
    except Exception as e:
        return failure_view(request,e)

@login_required(login_url='login')
def download_protein_by_organisms_from_entrez_search(request: WSGIRequest, search_id: int, organism_download: str):
    # starts a task that downloads the protein sequences in the entrez search results page of an choosen organism that is selected
    try:
        email = request.user.email
        download_organism_protein_sequences_task.delay(search_id,organism_download, email)

        return redirect('search_details', search_id=search_id)
    except Exception as e:
        return failure_view(request, e)

#unfinished
@login_required(login_url='login')
def delete_search_view(request: WSGIRequest, search_id: int):
    #delete edirectpaper file entry based on search_id
    #also deletes taskresult and entrezsearch db entry by search_id
    try:
        retcode = delete_esearch_by_id(search_id)
        if retcode == 0:

            return redirect("entrez_dashboard")
        elif retcode == 1:
            return failure_view(request, "Couldnt delete esearch object")
        else:
            return failure_view(request, "ERROR during deletion of edirectpaper with id: {}".format(search_id))

    except Exception as e:
        return failure_view(request,e)

@login_required(login_url='login')
def entrez_dashboard_view(request: WSGIRequest):
    # loads an entrez dashboard overview with a search bar in the upper section and searchresults in lower sections
    #if a search is added it is saved in the entrzsearch database and a download task is started based on it
    try:
        context = {}
        if request.method == "POST":
            entrez_search_form = EntrezSearchForm(request.POST)
            if entrez_search_form.is_valid():
                database = entrez_search_form.cleaned_data['database']
                entrez_query = entrez_search_form.cleaned_data['entrez_query'] #needs an if clause to check if in db already
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
def ajax_call_progress_entrezsearch_to_fasta(request,search_id: int):
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
        context['html_results'] = ''.join(get_html_results(project_id,str(query_sequence_id)+'/results_rbhs.html'))
        return render(request,"external_tools/phylogenetic_dashboard.html",context)
    except Exception as e:
        return failure_view(request,e)

'''cdd_domain_search_dashboard
    
    Dashboard page for the conserved domain search. It holds a form for rpsblast settings, the html table with 
    domains of the query sequences and the buttons for the detail/result pages and form submit. 
    
    :param project_id
        :type int

'''
@login_required(login_url='login')
def cdd_domain_search_dashboard(request, project_id):
    try:
        if request.method == "GET":
            query_sequences_rdy_for_cdd = ExternalTools.objects.get_cdd_searchable_queries(project_id=project_id)

            rps_blast_settings_form = RpsBLASTSettingsForm(query_sequences_rdy_for_cdd)
            query_sequence_cdd_search_dict = ExternalTools.objects.check_cdd_domain_search_task_status(project_id=project_id)
            context = {"query_task_dict":query_sequence_cdd_search_dict,
                       "project_id":project_id}
            context['html_results'] = ''.join(get_html_results(project_id,'query_domains.html'))
            context['rpsblast_settingsform'] = rps_blast_settings_form
            return render(request, "external_tools/cdd_domain_search_dashboard.html", context)
        else:
            return failure_view(request, "There is no post method for this view.")
    except Exception as e:
        return failure_view(request,e)


'''execute_cdd_domain_search_for_target_query

    This function executes the celery task for searching conserved domains within the CDD database. 
    Before the rpsblast is executed, the function validates if it is "practical" to execute the search.
    If the query sequence has just one domain or if there are only two reciprocal results the function is 
    not executed. Settings for the rpsblast are saved within a django form object. 
    The form is also validated. The query sequence identifier specified via a selection widget is then used
    as input for the rpsblast. 
    
    :param project:id
        :type int

'''
@login_required(login_url='login')
def execute_cdd_domain_search_for_target_query(request, project_id:int):
    try:
        if request.method == "POST":
            query_sequences_rdy_for_cdd = ExternalTools.objects.get_cdd_searchable_queries(project_id=project_id)

            rps_blast_form = RpsBLASTSettingsForm(query_sequences_rdy_for_cdd,request.POST)
            if rps_blast_form.is_valid():
                rps_blast_task_data = rps_blast_form.cleaned_data
                query_sequence = rps_blast_task_data['query_sequence']
                if check_if_cdd_search_can_get_executed(query_sequence, project_id) == 0:
                    cdd_domain_search_with_rbhs_task.delay(project_id, rps_blast_task_data)
            return redirect('cdd_domain_search_dashboard',project_id=project_id)
        else:
            return failure_view(request, "There is no other method than post for this view.")
    except Exception as e:
        return failure_view(request, e)


'''cdd_domain_search_details_view
    
    Function for loading the output of the cdd domain search task. Detail result page for each
    query sequence and successfull domain search tasks.
    
    :param query_id
        :type str
    :param project_id
        :type int

'''
@login_required(login_url='login')
def cdd_domain_search_details_view(request, query_id:str, project_id:int):
    try:
        #query_sequence_model = QuerySequences.objects.get(query_accession_id=)
        #okeh_plot = settings.BLAST_PROJECT_DIR + str(project_id) + '/' + query_id + '/pca_bokeh_domain_plot.html'
        context = {}
        context['query_id'] = query_id
        context['project_id'] = project_id
        context['CDDSearchPCABokehPlot'] = str(project_id) + '/' + query_id + '/pca_bokeh_domain_plot.html'
        return render(request, "external_tools/cdd_domain_search_details.html", context)
    except Exception as e:
        return failure_view(request, e)

'''delete_cdd_domain_search_view

    Function for deletion of the cdd_domain_search_task celery task result database object and all 
    associated files. The files to delete are defined as strings in a list within the delete_cdd_search_output
    function.
    
    :param query_id
        :type str
    :param project_id
        :type int
'''
@login_required(login_url='login')
def delete_cdd_domain_search_view(request, query_id:str, project_id:int):
    try:
        query_sequence = ExternalTools.objects.get_associated_query_sequence(project_id,query_id)
        #the query_sequence variable is a QuerySet but should just hold one query_sequence
        if len(query_sequence) > 1:
            raise Exception("[-] There are multiple query sequences with the name {} "
                            "in the associated ExternalTools model object, this should not happen as project creation"
                            "filters for duplicate entries ...".format(query_id))
        else:
            query_sequence[0].delete_cdd_search_task_result()
            delete_cdd_search_output(query_id, project_id)
            return redirect('cdd_domain_search_dashboard',project_id=project_id)
    except Exception as e:
        return failure_view(request, e)

'''get_cdd_task_status_ajax_call

    This function returns a json object with the cdd_domain_search_task result column for the 
    specified query sequence. The json object is used for displaying the progress of the task. 
    E.g.: {"pending": false, "current": 30, "total": 100, "percent": 30.0, "description": "PROGRESS"}

    :param query_id
        :type str
    :param project_id
        :type int
    
    :returns JsonResponse

'''
@login_required(login_url='login')
def get_cdd_task_status_ajax_call(request, query_id:str, project_id:int):
    try:
        if request.is_ajax and request.method == "GET":
            query_sequence = ExternalTools.objects.get_associated_query_sequence(project_id,query_id)[0]
            data = query_sequence.cdd_domain_search_task.result
            return JsonResponse({"data":loads(data)}, status=200)
        return JsonResponse({"ERROR":"NOT OK"},status=200)
    except Exception as e:
        return JsonResponse({"error": "{}".format(e)}, status=400)