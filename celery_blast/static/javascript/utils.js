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
  if (settingsForm.style.display == "none") {
    settingsForm.style.display = "block";
  } else {
    settingsForm.style.display = "none";
  }
}

function loadRemoteOrLocal(project_type) {
    document.getElementById('project_creation_' + project_type).style.display = "block";
    if (project_type == 'remote') {
        document.getElementById('project_creation_local').style.display = 'none'
        document.getElementById('local').checked = false

    } else if (project_type == 'local') {
        document.getElementById('project_creation_remote').style.display = 'none'
        document.getElementById('remote').checked = false
    }
}

function loadMultipleOrSingleFileUpload(file_upload_type) {
    if (document.getElementById('upload_form_master_div').style.display == 'none') {
        document.getElementById('upload_form_master_div').style.display = "block";
    }

    document.getElementById('post_form_' + file_upload_type + '_container').style.display = "block";
    if (file_upload_type == 'single') {
        document.getElementById('post_form_multiple_container').style.display = 'none'
        document.getElementById('multiple').checked = false
    } else if (file_upload_type == 'multiple') {
        document.getElementById('post_form_single_container').style.display = 'none'
        document.getElementById('single').checked = false
    }
}

function ajax_call_to_project_details(data, static_url, reciprocal_results_url) { // check if available
    var progress_bar = document.getElementById('progress_bar');
    var progress_container = document.getElementsByClassName('progress')[0];
    var progress = data.progress;

    var new_div_base_style = "position:absolute; \
      margin-right: auto; \
      color: white; \
      border-radius: 50% 50% 50% 50%;\
      background-repeat: no-repeat;\
      background-color: white;\
      background-position: center; ";

    var new_div_big =  new_div_base_style +
        "margin-top: 156px; height: 80px; width: 80px;" +
        "background-size: 64px 64px;"

    var new_div_small = new_div_base_style +
        "margin-top: 171px; height: 50px; width: 50px;"
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

    if(progress >= 25){
        var reciprocal_results_div_style = new_div_big
            + "margin-left: 20%;" + "background-size: 64px 64px;"
        var recblast_result_link_style = new_link_style_big
            + "background-image: url("
            +static_url+"festival-fireworks-icon.png"+");"

        function buildProgressButtonReciprocalBlastFinished(){
            var reciprocal_blast_finished = document.createElement("div");
            var recblast_result_link = document.createElement("a")
            recblast_result_link.href = "project_details/reciprocal_best_hits"
            recblast_result_link.target = "_blank"
            recblast_result_link.classList.add("btn")
            recblast_result_link.style.cssText = recblast_result_link_style
            reciprocal_blast_finished.appendChild(recblast_result_link)
            reciprocal_blast_finished.style.cssText = reciprocal_results_div_style

            progress_container.appendChild(reciprocal_blast_finished);
        }
        setTimeout(buildProgressButtonReciprocalBlastFinished, 1250)
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
            + "margin-left: 5%;"+ "background-size: 32px 32px;"
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