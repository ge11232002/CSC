
function validate_checkboxes(boxes, nr, max) {
  var count = 0;

  for (var j = 0; j < boxes.length; j++) {
    if (boxes[j].checked) {
      count++;
    }
  }


  if (boxes.length == undefined || count >= nr && (!max || count <= max)) {
    return true;
  }

  if (count <= nr) {
    errorMsg = "You need to select at least "+nr+" matri"+((nr == 1) ? "x" : "ces");
  } else {
    errorMsg = "We only allow a maximum of "+max+" matrices to be clustered. Use STAMP if you need more.";
  }

  alert(errorMsg);
  
  return false;
}

function Start(page, height, width) {
  OpenWin = this.open(page, "CtrlWindow", "height="+height+", width="+width+", toolbar=no,menubar=no,location=no,scrollbars=yes,resizable=yes");
}
	
function goto(page) { 
  window.location = page; 
}


