(function () {
  "use strict";
  $(document).on("shiny:inputchanged", function (event) {
    if (event.name !== "main_nav") return;
    var navbar = document.querySelector(".navbar");
    if (!navbar) return;
    var menu = navbar.querySelector(".navbar-collapse");
    var toggle = navbar.querySelector(".navbar-toggler, .navbar-toggle");
    if (toggle && menu && (menu.classList.contains("show") || menu.classList.contains("in"))) {
      toggle.click();
    }
    window.scrollTo({ top: 0, behavior: "auto" });
  });
})();
