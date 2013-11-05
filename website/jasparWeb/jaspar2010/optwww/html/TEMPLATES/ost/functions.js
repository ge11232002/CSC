
function validate_checkboxes(boxes, nr) {
  var count = 0;

  for (var j = 0; j < boxes.length; j++) {
    if (boxes[j].checked) {
      count++;
    }
  }

  if (count >= nr) {
    return true;
  }

  errorMsg = "You need to select at least "+nr+" matri"+((nr == 1) ? "x" : "ces");
  alert(errorMsg);
  
  return false;
}

function Start(page, height, width) {
  OpenWin = this.open(page, "CtrlWindow", "height="+height+", width="+width+", toolbar=no,menubar=no,location=no,scrollbars=yes,resizable=yes");
}
	
function goto(page) { 
  window.location = page; 
}
