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

function loadTemplate(project_type) {
    document.getElementById('project_creation_' + project_type).style.display = "block";
    if (project_type == 'remote') {
        document.getElementById('project_creation_local').style.display = 'none'
        document.getElementById('local').checked = false

    } else if (project_type == 'local') {
        document.getElementById('project_creation_remote').style.display = 'none'
        document.getElementById('remote').checked = false
    }
}