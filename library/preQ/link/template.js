window.onscroll = function() {scrollFunction()};

function scrollFunction() {
    if (document.getElementById("main").scrollTop > 20 || document.documentElement.scrollTop > 20) {
        document.getElementById("btpBtn").style.display = "block";
    } else {
        document.getElementById("btpBtn").style.display = "none";
    }
}

// When the user clicks on the button, scroll to the top of the document
function topFunction() {
    document.getElementById("main").scrollTop = 0;
    document.documentElement.scrollTop = 0;
}


  //href function
  function link(href){
    window.location.href = href;
  }