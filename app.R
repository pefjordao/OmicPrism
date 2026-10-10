# =========================================================
# OmicPrism - Shiny App (EN)
# =========================================================

# =========================================================
# 0) Package Installation and Loading
# =========================================================

auto_install_packages <- interactive()

cran_packages <- c(
  "shiny",
  "bslib",
  "shinyjs",
  "colourpicker",
  "tidyverse",
  "DT",
  "svglite",
  "cowplot",
  "msigdbr",
  "rstatix",
  "ggpubr",
  "eulerr",
  "shinyWidgets",
  "readxl"
)

bioc_packages <- c(
  "Mfuzz",
  "Biobase",
  "fgsea",
  "GSVA",
  "clusterProfiler",
  "EnhancedVolcano"
)

install_missing_cran <- function(packages) {
  
  missing_packages <- packages[
    !vapply(
      packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]
  
  if (length(missing_packages) > 0) {
    
    message(
      "Installing missing CRAN package(s): ",
      paste(missing_packages, collapse = ", ")
    )
    
    install.packages(
      pkgs = missing_packages,
      dependencies = TRUE
    )
  }
}

install_missing_bioc <- function(packages) {
  
  missing_packages <- packages[
    !vapply(
      packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]
  
  if (length(missing_packages) > 0) {
    
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    
    message(
      "Installing missing Bioconductor package(s): ",
      paste(missing_packages, collapse = ", ")
    )
    
    BiocManager::install(
      pkgs = missing_packages,
      ask = FALSE,
      update = FALSE
    )
  }
}

if (isTRUE(auto_install_packages)) {
  
  install_missing_cran(cran_packages)
  install_missing_bioc(bioc_packages)
}

library(shiny)
library(bslib)
library(shinyjs)
library(colourpicker)
library(tidyverse)

# =========================================================
# 1) Global App Settings
# =========================================================

options(
  shiny.maxRequestSize = 500 * 1024^2
)

# =========================================================
# 2) OmicPrism Theme
# =========================================================

app_theme <- bs_theme(
  version = 5,
  bg = "#f6f7fb",
  fg = "#17182d",
  primary = "#233e66",
  base_font = "system-ui, -apple-system, BlinkMacSystemFont, Segoe UI, sans-serif",
  heading_font = "system-ui, -apple-system, BlinkMacSystemFont, Segoe UI, sans-serif"
)

custom_css <- tags$link(rel = "stylesheet", type = "text/css", href = "omicprism.css")

# =========================================================
# 4) Module Loading
# =========================================================

module_files <- c(
  "modules/shared_helpers.R",
  "modules/overlap_module.R",
  "modules/pca_module.R",
  "modules/volcano_module.R",
  "modules/heatmap_module.R",
  "modules/mfuzz_module.R",
  "modules/ora_module.R",
  "modules/ssGSEA_module.R",
  "modules/single_GSEA_module.R",
  "modules/comparative_GSEA_module.R"
)

check_module_files <- function(files) {
  
  missing_files <- files[!file.exists(files)]
  
  if (length(missing_files) > 0) {
    
    stop(
      paste0(
        "Missing module file(s):\n",
        paste(missing_files, collapse = "\n")
      ),
      call. = FALSE
    )
  }
}

source_module_files <- function(files) {
  
  check_module_files(files)
  
  # Keep each module's helpers private; exported UI/server closures retain them.
  for (file in files) {
    if (basename(file) == "shared_helpers.R") {
      sys.source(file, envir = .GlobalEnv)
    } else {
      module_env <- new.env(parent = .GlobalEnv)
      sys.source(file, envir = module_env)
      exports <- grep("(ModuleUI|ModuleServer|CompareUI|CompareServer)$",
                      ls(module_env), value = TRUE)
      for (name in exports) {
        if (is.function(module_env[[name]])) {
          assign(name, module_env[[name]], envir = .GlobalEnv)
        }
      }
    }
  }
}

source_module_files(module_files)

# =========================================================
# 4.1) Module Function Resolution
# =========================================================

resolve_module_function <- function(candidates) {
  
  for (candidate in candidates) {
    
    if (exists(candidate, mode = "function", envir = .GlobalEnv)) {
      return(get(candidate, envir = .GlobalEnv))
    }
  }
  
  stop(
    paste0(
      "None of the expected module functions were found: ",
      paste(candidates, collapse = ", ")
    ),
    call. = FALSE
  )
}

overlap_ui <- resolve_module_function(c("overlapModuleUI", "vennModuleUI"))
overlap_server <- resolve_module_function(c("overlapModuleServer", "vennModuleServer"))

pca_ui <- resolve_module_function(c("pcaModuleUI"))
pca_server <- resolve_module_function(c("pcaModuleServer"))

volcano_ui <- resolve_module_function(c("volcanoModuleUI"))
volcano_server <- resolve_module_function(c("volcanoModuleServer"))

heatmap_ui <- resolve_module_function(c("heatmapModuleUI"))
heatmap_server <- resolve_module_function(c("heatmapModuleServer"))

mfuzz_ui <- resolve_module_function(c("mfuzzModuleUI"))
mfuzz_server <- resolve_module_function(c("mfuzzModuleServer"))

ora_ui <- resolve_module_function(c("oraModuleUI"))
ora_server <- resolve_module_function(c("oraModuleServer"))

ssgsea_ui <- resolve_module_function(c("ssGSEAModuleUI", "ssgseaModuleUI"))
ssgsea_server <- resolve_module_function(c("ssGSEAModuleServer", "ssgseaModuleServer"))

single_gsea_ui <- resolve_module_function(c("singleGSEAModuleUI", "singleGseaModuleUI", "gseaModuleUI"))
single_gsea_server <- resolve_module_function(c("singleGSEAModuleServer", "singleGseaModuleServer", "gseaModuleServer"))

comparative_gsea_ui <- resolve_module_function(c("comparativeGSEAModuleUI", "comparativeGseaModuleUI", "gseaCompareUI"))
comparative_gsea_server <- resolve_module_function(c("comparativeGSEAModuleServer", "comparativeGseaModuleServer", "gseaCompareServer"))

# =========================================================
# 5) Home Page and Example Files
# =========================================================

source("ui/home.R", local = TRUE)

# =========================================================
# 6) Main UI
# =========================================================

ui <- tagList(
  
  useShinyjs(),
  custom_css,
  tags$head(
    tags$link(rel = "icon", type = "image/png", href = "omicprism-icon.png"),
    tags$meta(name = "description", content = "OmicPrism: omics visualization, time-series clustering and pathway enrichment."),
    tags$script(src = "omicprism.js")
  ),
  
  navbarPage(
    id = "main_nav",
    title = tags$span(class = "brand-lockup",
      tags$img(src = "omicprism-icon.png", class = "brand-icon", alt = "OmicPrism", width = 34, height = 40),
      tags$span(class = "brand-name", "OmicPrism")
    ),
    windowTitle = "OmicPrism",
    theme = app_theme,
    position = "static-top",
    collapsible = TRUE,
    
    home_page_ui(),
    
    tabPanel(
      title = "Overlap",
      value = "overlap",
      icon = icon("circle-half-stroke"),
      div(class = "module-workspace", overlap_ui("overlap_id"))
    ),
    
    tabPanel(
      title = "PCA",
      value = "pca",
      icon = icon("chart-simple"),
      div(class = "module-workspace", pca_ui("pca_id"))
    ),
    
    tabPanel(
      title = "Volcano Plot",
      value = "volcano",
      icon = icon("mountain"),
      div(class = "module-workspace", volcano_ui("volcano_id"))
    ),
    
    tabPanel(
      title = "Heatmap",
      value = "heatmap",
      icon = icon("table-cells"),
      div(class = "module-workspace", heatmap_ui("heatmap_id"))
    ),
    
    tabPanel(
      title = "Time-Series (Mfuzz)",
      value = "mfuzz",
      icon = icon("chart-line"),
      div(class = "module-workspace", mfuzz_ui("mfuzz_id"))
    ),
    
    navbarMenu(
      title = "Enrichment analysis",
      icon = icon("flask"),
      
      tabPanel(
        title = "ORA",
        value = "ora",
        icon = icon("list-check"),
        div(class = "module-workspace", ora_ui("ora_id"))
      ),
      
      tabPanel(
        title = "ssGSEA",
        value = "ssgsea",
        icon = icon("chart-column"),
        div(class = "module-workspace", ssgsea_ui("ssgsea_id"))
      ),
      
      tabPanel(
        title = "Single GSEA",
        value = "gsea_single",
        icon = icon("chart-line"),
        div(class = "module-workspace", single_gsea_ui("gsea_single"))
      ),
      
      tabPanel(
        title = "Comparative GSEA",
        value = "gsea_comp",
        icon = icon("layer-group"),
        div(class = "module-workspace", comparative_gsea_ui("gsea_comp"))
      )
    )
  ),
  
  tags$footer(
    class = "app-footer",
    tags$span(paste0("OmicPrism \u00a9 ", format(Sys.Date(), "%Y"), " \u00b7 v1.0")),
    tags$span("Omics visualization & enrichment analysis")
  )
)

# =========================================================
# 7) Server
# =========================================================

server <- function(input, output, session) {
  register_example_downloads(output)
  
  card_navigation <- list(
    go_overlap     = "overlap",
    go_pca         = "pca",
    go_volcano     = "volcano",
    go_heatmap     = "heatmap",
    go_mfuzz       = "mfuzz",
    go_ora         = "ora",
    go_ssgsea      = "ssgsea",
    go_gsea_single = "gsea_single",
    go_gsea_comp   = "gsea_comp"
  )
  
  purrr::iwalk(
    card_navigation,
    function(tab_value, input_id) {
      
      observeEvent(input[[input_id]], {
        
        updateNavbarPage(
          session = session,
          inputId = "main_nav",
          selected = tab_value
        )
        
      }, ignoreInit = TRUE)
    }
  )
  
  overlap_server("overlap_id")
  
  pca_server("pca_id")
  
  volcano_server("volcano_id")
  
  heatmap_server("heatmap_id")
  
  mfuzz_server("mfuzz_id")
  
  ora_server("ora_id")
  
  ssgsea_server("ssgsea_id")
  
  single_gsea_server("gsea_single")
  
  comparative_gsea_server("gsea_comp")
}

# =========================================================
# 8) Run App
# =========================================================

shinyApp(
  ui = ui,
  server = server
)
