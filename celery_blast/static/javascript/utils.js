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
    document.getElementById('post_form_'+file_upload_type+'_container').style.display = "block";
    if(file_upload_type == 'single'){
        document.getElementById('post_form_multiple_container').style.display = 'none'
        document.getElementById('multiple').checked = false
    } else if (file_upload_type == 'multiple'){
        document.getElementById('post_form_single_container').style.display = 'none'
        document.getElementById('single').checked = false
    }
}