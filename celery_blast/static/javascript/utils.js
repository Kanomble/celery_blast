function showLoader(idLoadingDiv,idHidingDiv) {
    //declare variables
    var loadingDiv = document.getElementById(idLoadingDiv)
    var hidingDiv = document.getElementById(idHidingDiv)

    //hide selected div and show giphy loading div
    hidingDiv.style.display = "none";
    loadingDiv.style.display = "block";
};

function displayBlastSettings(divElement) {
  var settingsForm = document.getElementById(divElement);
  if (settingsForm.style.display == "none") {
    settingsForm.style.display = "block";
  } else {
    settingsForm.style.display = "none";
  }
}