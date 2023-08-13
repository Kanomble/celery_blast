function showLoader(idLoadingDiv,idHidingDiv) {
    //declare variables
    var loadingDiv = document.getElementById(idLoadingDiv)
    var hidingDiv = document.getElementById(idHidingDiv)

    //hide selected div and show giphy loading div
    hidingDiv.style.display = "none";
    loadingDiv.style.display = "block";
};

function displayDivElement(divElement) {
  var settingsForm = document.getElementById(divElement);
  if (settingsForm.style.display === "none") {
    settingsForm.style.display = "block";
  } else {
    settingsForm.style.display = "none";
  }
}

function loadRemoteOrLocal(project_type) {
    document.getElementById('project_creation_' + project_type).style.display = "block";
    if (project_type === 'remote') {
        document.getElementById('project_creation_local').style.display = 'none';
        document.getElementById('local').checked = false;
        document.getElementById('remote').checked = true;
    } else if (project_type === 'local') {
        document.getElementById('project_creation_remote').style.display = 'none';
        document.getElementById('remote').checked = false;
        document.getElementById('local').checked = true;
    }
}

// function to switch between the two different input forms for uploading custom protein fasta files
function loadMultipleOrSingleFileUpload(file_upload_type) {
    if (document.getElementById('upload_form_master_div').style.display === 'none') {
        document.getElementById('upload_form_master_div').style.display = "block";
    }

    document.getElementById('post_form_' + file_upload_type + '_container').style.display = "block";
    if (file_upload_type === 'single') {
        document.getElementById('post_form_multiple_container').style.display = 'none';
        document.getElementById('multiple').checked = false;
        document.getElementById('single').checked = true;
    } else if (file_upload_type === 'multiple') {
        document.getElementById('post_form_single_container').style.display = 'none';
        document.getElementById('single').checked = false;
        document.getElementById('multiple').checked = true;
    }
}
// function for creating interactive elements within the progress bar of the reciprocal blast project detail page
function ajax_call_to_project_details(data, static_url, reciprocal_results_url, query_sequence_info_url) { // check if available
    var progress_bar = document.getElementById('progress_bar');
    var progress_container = document.getElementsByClassName('progress')[0];
    var progress = data.progress;

    var dropdown_content_div = "min-width:300px;box-shadow: 0px 8px 16px 0px rgba(0,0,0,0.2);"

    var new_div_base_style = "position:absolute; \
      margin-right: auto; \
      color: white; \
      border-radius: 50% 50% 50% 50%;\
      background-repeat: no-repeat;\
      background-color: white;\
      background-position: center; ";

    var new_div_big =  new_div_base_style +
        "margin-top: -20px; height: 80px; width: 80px;" +
        "background-size: 64px 64px;"

    var new_div_small = new_div_base_style +
        "margin-top: -8px; height: 50px; width: 50px;"
        + "background-size: 32px 32px;";

    var new_link_style = "\
      position: absolute;\
      text-align: center;border-radius: 50% 50% 50% 50%;\
      border: 2px solid #d1d5db;\
      background-repeat: no-repeat;\
      background-color: white;\
      background-position: center;";

    var new_link_style_big = new_link_style +
        "width:80px;height: 80px;" +
        "background-size: 64px 64px;"

    var new_link_style_small = new_link_style +
        "width:50px;height: 50px;"
        +"background-size: 32px 32px;";

    if(progress >= 70){
        var query_sequence_infos_div = new_div_small
            + "margin-left: 65%;" + "background-image: url("
            +static_url+"festival-fireworks-icon.png"+");" + "border: 2px solid #d1d5db;"

        function buildProgressButtonQueryInfoFinished(){
            var query_info_finished = document.createElement("div");

            query_info_finished.style.cssText = query_sequence_infos_div
            query_info_finished.classList.add("dropdown_menu")
            var query_info_dropdown_menu = document.createElement("div")

            var query_info_link = document.createElement("a")
            var query_info_log_link = document.createElement("a")
            var blast_tables_to_plot_link = document.createElement("a")
            query_info_link.text = "Query Information"
            query_info_log_link.text = "Query Information Log"
            blast_tables_to_plot_link.text = "BLAST Tables To Plots Log"



            query_info_link.href = query_sequence_info_url
            query_info_log_link.href = "project_details/query_sequences_to_html_table"
            blast_tables_to_plot_link.href = "project_details/blast_tables_to_plots"

            query_info_link.target = "_blank"
            query_info_log_link.target = "_blank"
            blast_tables_to_plot_link.target = "_blank"

            query_info_dropdown_menu.style.cssText = dropdown_content_div
            query_info_dropdown_menu.classList.add("dropdown_content")
            query_info_dropdown_menu.appendChild(query_info_link)
            query_info_dropdown_menu.appendChild(query_info_log_link)
            query_info_dropdown_menu.appendChild(blast_tables_to_plot_link)

            query_info_finished.appendChild(query_info_dropdown_menu)

            progress_container.appendChild(query_info_finished);
        }
        setTimeout(buildProgressButtonQueryInfoFinished, 1050)
    }

    if(progress >= 50){
        var reciprocal_results_div_style = new_div_big
            + "margin-left: 45%;" + "background-size: 64px 64px;" + "background-image: url("
            +static_url+"festival-fireworks-icon.png"+");" + "border: 2px solid #d1d5db;"

        function buildProgressButtonReciprocalBlastFinished(){
            var reciprocal_blast_finished = document.createElement("div");

            reciprocal_blast_finished.style.cssText = reciprocal_results_div_style
            reciprocal_blast_finished.classList.add("dropdown_menu")
            var dropdown_menu = document.createElement("div")
            var recblast_result_link = document.createElement("a")
            var recblast_result_table_link = document.createElement("a")
            var backward_blast_log = document.createElement("a")
            var fw_result_preparation_log = document.createElement("a")
            var build_folders_with_hit_info_for_each_qseqid_log = document.createElement("a")
            var blast_tables_to_csv_log = document.createElement("a")

            recblast_result_link.text = "Reciprocal BLAST Logfile"
            recblast_result_table_link.text = "RBH Result Table"
            backward_blast_log.text = "Backward BLAST Results Log"
            fw_result_preparation_log.text = "Forward BLAST Results Log"
            build_folders_with_hit_info_for_each_qseqid_log.text = "Subdirectories For Targets Log"
            blast_tables_to_csv_log.text = "RBH Taxonomy Inference Log"

            recblast_result_table_link.href = reciprocal_results_url
            recblast_result_link.href = "project_details/reciprocal_best_hits"
            backward_blast_log.href = "project_details/backward_blast"
            fw_result_preparation_log.href = "project_details/fw_result_processing"
            build_folders_with_hit_info_for_each_qseqid_log.href = "project_details/build_folders_with_hit_info_for_each_qseqid"
            blast_tables_to_csv_log.href = "project_details/blast_tables_to_csv"

            recblast_result_link.target = "_blank"
            recblast_result_table_link.target = "_blank"
            backward_blast_log.target = "_blank"
            fw_result_preparation_log.target = "_blank"
            build_folders_with_hit_info_for_each_qseqid_log.target = "_blank"
            blast_tables_to_csv_log.target = "_blank"

            //recblast_result_link.classList.add("btn")
            dropdown_menu.classList.add("dropdown_content")
            dropdown_menu.appendChild(blast_tables_to_csv_log)
            dropdown_menu.appendChild(build_folders_with_hit_info_for_each_qseqid_log)
            dropdown_menu.appendChild(recblast_result_table_link)
            dropdown_menu.appendChild(recblast_result_link)
            dropdown_menu.appendChild(backward_blast_log)
            dropdown_menu.appendChild(fw_result_preparation_log)

            reciprocal_blast_finished.appendChild(dropdown_menu)

            progress_container.appendChild(reciprocal_blast_finished);
        }
        setTimeout(buildProgressButtonReciprocalBlastFinished, 1050)
    }

    if(progress >= 20){
        var backward_blast_style = new_div_small
            + "margin-left: 15%;"
        var backward_blast_link_style = new_link_style_small
            + "background-image: url("
            +static_url+"festival-fireworks-icon.png"+");"


        function buildProgressButtonBackwardBlastFinished(){
            var backward_blast_finished = document.createElement("div");
            var backward_blast_link = document.createElement("a")

            backward_blast_link.href = "project_details/backward_blast"
            backward_blast_link.target = "_blank"
            backward_blast_link.classList.add("btn")
            backward_blast_link.style.cssText = backward_blast_link_style
            backward_blast_finished.appendChild(backward_blast_link)
            backward_blast_finished.style.cssText = backward_blast_style
            progress_container.appendChild(backward_blast_finished);
        }

        setTimeout(buildProgressButtonBackwardBlastFinished, 500)
    }

    if(progress >= 5){
        var new_div_style = new_div_small
            + "margin-left: 1%;"+ "background-size: 32px 32px;"
        var link_style = new_link_style_small
            + "background-image: url("
            +static_url+"festival-fireworks-icon.png"+");"
            +"background-size: 32px 32px;";

        function buildProgressButtonForwardBlastFinished(){
            var blast_finished = document.createElement("div");
            var result_link = document.createElement("a")

            result_link.href = "project_details/forward_blast"
            result_link.target = "_blank"
            result_link.classList.add("btn")
            result_link.style.cssText = link_style
            blast_finished.appendChild(result_link)
            blast_finished.style.cssText = new_div_style
            progress_container.appendChild(blast_finished);
        }

        setTimeout(buildProgressButtonForwardBlastFinished, 100)
    }

    progress_bar.style.width = String(progress) + "%";
    progress_bar.innerHTML = String(progress) + "%";

 }