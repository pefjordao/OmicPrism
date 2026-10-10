# =========================================================
# MODULE: ORA — Over-Representation Analysis (EN)
# =========================================================

library(shiny)
library(shinyjs)
library(tidyverse)
library(DT)
library(svglite)
library(colourpicker)
library(clusterProfiler)
library(msigdbr)
library(readxl)

# =========================================================
# Helper Functions
# =========================================================

ora_clean_column_names <- function(df) {
  
  colnames(df) <- colnames(df) %>%
    as.character() %>%
    stringr::str_trim()
  
  df
}

ora_try_read_delim <- function(path, delim, decimal_mark = ".") {
  
  tryCatch(
    {
      readr::read_delim(
        file = path,
        delim = delim,
        name_repair = "minimal",
        show_col_types = FALSE,
        col_types = readr::cols(.default = readr::col_character()),
        trim_ws = TRUE,
        locale = readr::locale(decimal_mark = decimal_mark)
      ) %>%
        omics_delimited_frame() %>%
        ora_clean_column_names()
    },
    error = function(e) NULL
  )
}

ora_try_read_whitespace <- function(path, decimal_mark = ".") {
  
  tryCatch(
    {
      readr::read_table(
        file = path,
        show_col_types = FALSE,
        col_types = readr::cols(.default = readr::col_character()),
        trim_ws = TRUE,
        locale = readr::locale(decimal_mark = decimal_mark)
      ) %>%
        omics_delimited_frame() %>%
        ora_clean_column_names()
    },
    error = function(e) NULL
  )
}

ora_try_read_excel <- function(path) {
  
  tryCatch(
    {
      readxl::read_excel(
        path = path,
        sheet = 1,
          .name_repair = "minimal"
      ) %>%
        as.data.frame() %>%
        ora_clean_column_names()
    },
    error = function(e) NULL
  )
}

ora_score_table <- function(df) {
  
  if (is.null(df) || !is.data.frame(df) || ncol(df) <= 0) {
    return(-Inf)
  }
  
  non_empty_names <- sum(!is.na(colnames(df)) & colnames(df) != "")
  
  ncol(df) * 1000 + non_empty_names
}

ora_smart_read_unchecked <- function(path, original_name = NULL, min_cols = 1) {
  
  ext_source <- original_name %||% path
  ext <- tolower(tools::file_ext(ext_source))
  
  if (ext %in% c("xls", "xlsx")) {
    
    validate(
      need(
        requireNamespace("readxl", quietly = TRUE),
        "Package 'readxl' is required to read Excel files."
      )
    )
    
    df <- ora_try_read_excel(path)
    
    validate(
      need(
        !is.null(df) && ncol(df) >= min_cols,
        "Invalid Excel file. The first sheet must contain a valid table."
      )
    )
    
    return(as.data.frame(df))
  }
  
  read_attempts <- list(
    comma_dot = ora_try_read_delim(
      path = path,
      delim = ",",
      decimal_mark = "."
    ),
    
    semicolon_dot = ora_try_read_delim(
      path = path,
      delim = ";",
      decimal_mark = "."
    ),
    
    semicolon_comma = ora_try_read_delim(
      path = path,
      delim = ";",
      decimal_mark = ","
    ),
    
    tab_dot = ora_try_read_delim(
      path = path,
      delim = "\t",
      decimal_mark = "."
    ),
    
    tab_comma = ora_try_read_delim(
      path = path,
      delim = "\t",
      decimal_mark = ","
    ),
    
    whitespace_dot = ora_try_read_whitespace(
      path = path,
      decimal_mark = "."
    ),
    
    whitespace_comma = ora_try_read_whitespace(
      path = path,
      decimal_mark = ","
    )
  )
  
  scores <- vapply(
    read_attempts,
    ora_score_table,
    numeric(1)
  )
  
  best_name <- names(which.max(scores))
  df <- read_attempts[[best_name]]
  
  validate(
    need(
      !is.null(df) && ncol(df) >= min_cols,
      paste(
        "The uploaded file could not be parsed.",
        "Use a file separated by comma, semicolon, tab, or spaces,",
        "or upload a valid Excel file with a header row."
      )
    )
  )
  
  as.data.frame(df)
}

ora_smart_read <- function(...) {
  omics_validate_table(ora_smart_read_unchecked(...))
}

as_numeric_flexible <- omics_numeric


strip_ansi_codes <- function(x) {
  
  x <- as.character(x)
  
  x <- gsub(
    pattern = "\033\\[[0-9;]*[[:alpha:]]",
    replacement = "",
    x = x,
    perl = TRUE
  )
  
  x <- gsub(
    pattern = "[\001-\010\013\014\016-\037\177]",
    replacement = "",
    x = x,
    perl = TRUE
  )
  
  x
}


guess_column <- function(cols, pattern, default = NULL) {
  
  matched <- grep(
    pattern = pattern,
    x = cols,
    value = TRUE,
    ignore.case = TRUE
  )
  
  if (length(matched) == 0) {
    default %||% cols[1]
  } else {
    matched[1]
  }
}

clean_feature_ids <- function(x) {
  
  x <- x %>%
    as.character() %>%
    stringr::str_trim()
  
  unique(x[!is.na(x) & x != ""])
}

safe_neglog10 <- function(x) {
  -log10(pmax(x, .Machine$double.xmin))
}

make_ora_overlap_diagnostics <- function(selected_ids, universe_ids, term2gene_df) {
  
  selected_ids <- clean_feature_ids(selected_ids)
  universe_ids <- clean_feature_ids(universe_ids)
  gene_set_ids <- clean_feature_ids(term2gene_df$gene)
  
  selected_exact <- length(intersect(selected_ids, gene_set_ids))
  universe_exact <- length(intersect(universe_ids, gene_set_ids))
  
  selected_case_insensitive <- length(
    intersect(
      stringr::str_to_upper(selected_ids),
      stringr::str_to_upper(gene_set_ids)
    )
  )
  
  universe_case_insensitive <- length(
    intersect(
      stringr::str_to_upper(universe_ids),
      stringr::str_to_upper(gene_set_ids)
    )
  )
  
  list(
    selected_exact = selected_exact,
    universe_exact = universe_exact,
    selected_case_insensitive = selected_case_insensitive,
    universe_case_insensitive = universe_case_insensitive
  )
}

format_ora_mapping_error <- function(diag) {
  
  paste0(
    "No gene set remained after ID mapping and set-size filtering.\n\n",
    "Check species/ID format or change Min/Max set size.\n\n",
    "Selected overlap: ",
    diag$selected_exact,
    " | Universe overlap: ",
    diag$universe_exact,
    "\nCase-insensitive selected overlap: ",
    diag$selected_case_insensitive,
    " | Case-insensitive universe overlap: ",
    diag$universe_case_insensitive
  )
}

show_ora_error <- function(message) {
  
  showModal(
    modalDialog(
      title = "ORA failed",
      div(
        style = "white-space: pre-wrap;",
        strip_ansi_codes(message)
      ),
      easyClose = TRUE,
      footer = modalButton("OK")
    )
  )
}

ratio_to_numeric <- function(x) {
  
  sapply(x, function(z) {
    
    parts <- strsplit(z, "/", fixed = TRUE)[[1]]
    
    if (length(parts) != 2) {
      return(NA_real_)
    }
    
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
}

read_universe_file <- function(path, original_name = NULL) {
  omics_read_id_list(path, original_name)
}

get_msigdb_terms <- function(species, category, subcategory = NULL) {
  
  msigdbr_args <- list(
    species = species
  )
  
  formals_names <- names(formals(msigdbr::msigdbr))
  
  if ("collection" %in% formals_names) {
    msigdbr_args$collection <- category
    if (!is.null(subcategory)) msigdbr_args$subcollection <- subcategory
  } else {
    msigdbr_args$category <- category
    if (!is.null(subcategory)) msigdbr_args$subcategory <- subcategory
  }
  
  msig <- do.call(msigdbr::msigdbr, msigdbr_args)
  
  validate(
    need(
      all(c("gs_name", "gene_symbol") %in% colnames(msig)),
      "MSigDB output does not contain expected columns: gs_name and gene_symbol."
    )
  )
  
  msig %>%
    dplyr::select(term = gs_name, gene = gene_symbol) %>%
    dplyr::distinct() %>%
    dplyr::filter(!is.na(term), !is.na(gene), term != "", gene != "")
}

read_gmt_terms <- function(path) {
  
  gmt <- clusterProfiler::read.gmt(path)
  
  validate(
    need(
      all(c("term", "gene") %in% colnames(gmt)),
      "Invalid GMT file. Expected columns: term and gene."
    )
  )
  
  gmt %>%
    dplyr::select(term, gene) %>%
    dplyr::distinct() %>%
    dplyr::filter(!is.na(term), !is.na(gene), term != "", gene != "")
}

format_gene_preview <- function(x, n = 30) {
  
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    return("")
  }
  
  genes <- unlist(strsplit(as.character(x), "/", fixed = TRUE))
  genes <- genes[!is.na(genes) & genes != ""]
  
  paste(
    head(genes, n),
    collapse = ", "
  )
}

open_ora_plot_device <- function(file, format, width, height, dpi) {
  
  if (format == "png") {
    
    png(
      filename = file,
      width = width,
      height = height,
      units = "in",
      res = dpi
    )
    
  } else if (format == "pdf") {
    
    pdf(
      file = file,
      width = width,
      height = height
    )
    
  } else if (format == "tiff") {
    
    tiff(
      filename = file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      compression = "lzw"
    )
    
  } else if (format == "svg") {
    
    svglite::svglite(
      file = file,
      width = width,
      height = height
    )
  }
}

# =========================================================
# UI
# =========================================================

oraModuleUI <- function(id) {
  
  ns <- NS(id)
  
  tagList(
    
    useShinyjs(),
    
    tags$head(
      tags$style(HTML("
        .ora-panel {
          background: #ffffff;
          border-radius: 8px;
          padding: 20px;
          border: 1px solid #eeeeee;
          margin-bottom: 20px;
        }

        .ora-help {
          color: #6b7280;
          font-size: 0.90rem;
          margin-bottom: 12px;
        }

        .well {
          background-color: #ffffff;
          border: 1px solid #e3e3df;
          border-radius: 8px;
          box-shadow: none;
        }

        table.dataTable td.dt-term-cell {
          white-space: normal !important;
          word-break: break-word;
          min-width: 320px;
          max-width: 520px;
        }

        table.dataTable td.dt-gene-cell {
          white-space: normal !important;
          word-break: break-word;
          min-width: 280px;
          max-width: 420px;
        }

        table.dataTable td.dt-center-cell {
          text-align: center;
          vertical-align: middle;
        }
      "))
    ),
    
    fluidRow(
      
      # ===================================================
      # Left Sidebar
      # ===================================================
      
      column(
        width = 4,
        
        wellPanel(
          
          # -------------------------------------------------
          # 1) Data Input
          # -------------------------------------------------
          
          h4("1. Data Input"),
          
          fileInput(
            ns("upload"),
            "Upload feature table (.txt, .csv, .tsv, .xls, .xlsx)",
            accept = c(".txt", ".csv", ".tsv", ".xls", ".xlsx")
          ),
          
          radioButtons(
            ns("input_mode"),
            "Input table type:",
            choices = c(
              "Differential results table" = "differential",
              "Feature list / cluster table" = "feature_list"
            ),
            selected = "differential"
          ),
          
          uiOutput(ns("column_selector")),
          
          uiOutput(ns("feature_group_selector")),
          
          p(
            class = "ora-help",
            "Differential results tables should contain feature IDs, log2FC, and nominal or adjusted P-values. Feature lists or cluster tables only require a feature ID column, with optional filtering by a grouping column such as cluster assignment."
          ),

          hr(),

          # -------------------------------------------------
          # 2) Feature Filtering
          # -------------------------------------------------

          h4("2. Feature Filtering"),

          uiOutput(ns("feature_filtering_ui")),

          hr(),
          
          # -------------------------------------------------
          # 3) Universe / Background
          # -------------------------------------------------
          
          h4("3. Universe / Background"),
          
          radioButtons(
            ns("universe_source"),
            "Universe source:",
            choices = c(
              "All features in uploaded table" = "all_uploaded",
              "Upload custom universe/background" = "custom"
            ),
            selected = "all_uploaded"
          ),
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'custom'", ns("universe_source")),
            fileInput(
              ns("universe_file"),
              "Upload universe/background file (.txt, .csv, .tsv, .xls, .xlsx)",
              accept = c(".txt", ".csv", ".tsv", ".xls", ".xlsx")
            )
          ),
          
          p(
            class = "ora-help",
            "The universe should be the set of features that could have been selected by the experiment and analysis. For proteomics, this is usually all quantified or tested proteins. For cluster tables, uploading a broader custom universe is recommended when available."
          ),
          
          hr(),
          
          # -------------------------------------------------
          # 4) Gene Set Database
          # -------------------------------------------------
          
          h4("4. Gene Set Database"),
          
          radioButtons(
            ns("gene_set_source"),
            "Gene set source:",
            choices = c(
              "MSigDB" = "msigdb",
              "Custom GMT file" = "gmt"
            ),
            selected = "msigdb"
          ),
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'msigdb'", ns("gene_set_source")),
            
            selectInput(
              ns("species"),
              "Species:",
              choices = c(
                "Homo sapiens",
                "Mus musculus",
                "Rattus norvegicus"
              ),
              selected = "Homo sapiens"
            ),
            
            selectInput(
              ns("category"),
              "Main category:",
              choices = msigdb_category_choices(),
              selected = "H"
            ),
            
            selectInput(
              ns("subcategory"),
              "Subcategory:",
              choices = NULL
            )
          ),
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'gmt'", ns("gene_set_source")),
            
            fileInput(
              ns("gmt_file"),
              "Upload GMT file",
              accept = c(".gmt")
            )
          ),
          
          hr(),
          
          # -------------------------------------------------
          # 5) ORA Parameters
          # -------------------------------------------------
          
          h4("5. ORA Parameters"),
          
          fluidRow(
            column(
              6,
              numericInput(
                ns("min_size"),
                "Min set size:",
                value = 10,
                min = 1
              )
            ),
            column(
              6,
              numericInput(
                ns("max_size"),
                "Max set size:",
                value = 500,
                min = 2
              )
            )
          ),
          
          selectInput(
            ns("padj_method"),
            "P-value adjustment method:",
            choices = c(
              "BH",
              "BY",
              "bonferroni",
              "holm",
              "hochberg",
              "hommel"
            ),
            selected = "BH"
          ),
          
          radioButtons(
            ns("result_p_source"),
            "Use for ORA result filtering:",
            choices = c(
              "Adjusted P-value / FDR" = "adjusted",
              "Nominal P-value" = "nominal"
            ),
            selected = "adjusted"
          ),
          
          numericInput(
            ns("result_p_cutoff"),
            "ORA result P-value/FDR cutoff:",
            value = 0.05,
            min = 0,
            max = 1,
            step = 0.01
          ),

          actionButton(
            ns("run_ora"),
            "Run ORA",
            class = "btn-primary",
            style = "width:100%; font-weight:bold;"
          ),

          hr(),

          # -------------------------------------------------
          # 6) Plot Options
          # -------------------------------------------------
          
          h4("6. Plot Options"),

          textInput(
            ns("plot_title"),
            "Plot Title",
            value = "ORA"
          ),

          selectInput(
            ns("plot_type"),
            "Plot type:",
            choices = c(
              "Dot plot" = "dot",
              "Bar plot" = "bar"
            ),
            selected = "dot"
          ),
          
          sliderInput(
            ns("top_n"),
            "Top N terms:",
            min = 5,
            max = 50,
            value = 20,
            step = 1
          ),
          
          sliderInput(
            ns("font_size"),
            "Font size:",
            min = 8,
            max = 18,
            value = 12,
            step = 1
          ),
          
          sliderInput(
            ns("wrap_width"),
            "Term label wrap width:",
            min = 20,
            max = 80,
            value = 45,
            step = 5
          ),
          
          hr(),
          
          # -------------------------------------------------
          # 7) Colors
          # -------------------------------------------------
          
          h4("7. Colors"),
          
          colourInput(
            ns("col_low"),
            "Low significance",
            "#2A78D6"
          ),

          colourInput(
            ns("col_high"),
            "High significance",
            "#E34948"
          ),
          
          hr(),
          
          # -------------------------------------------------
          # 8) Export
          # -------------------------------------------------
          
          h4("8. Export"),
          
          fluidRow(
            column(
              6,
              numericInput(
                ns("exp_width"),
                "Width (cm)",
                value = 24,
                min = 5,
                max = 60
              )
            ),
            column(
              6,
              numericInput(
                ns("exp_height"),
                "Height (cm)",
                value = 14,
                min = 5,
                max = 60
              )
            )
          ),

          fluidRow(
            column(
              6,
              selectInput(
                ns("exp_fmt"),
                "Format",
                choices = c("pdf", "png", "tiff", "svg"),
                selected = "svg"
              )
            ),
            column(
              6,
              conditionalPanel(
                condition = sprintf("input['%s'] == 'png' || input['%s'] == 'tiff'", ns("exp_fmt"), ns("exp_fmt")),
                numericInput(
                  ns("exp_dpi"),
                  "DPI",
                  value = 300,
                  min = 72
                )
              )
            )
          ),
          
          downloadButton(
            ns("download_plot"),
            "Download plot",
            class = "btn-success",
            style = "width:100%"
          ),
          
          br(),
          br(),
          
          downloadButton(
            ns("download_table"),
            "Download ORA results table",
            class = "btn-info",
            style = "width:100%"
          ),
          
          br(),
          br(),
          
          downloadButton(
            ns("download_features"),
            "Download selected features",
            class = "btn-secondary",
            style = "width:100%"
          )
        )
      ),
      
      # ===================================================
      # Right Output Area
      # ===================================================
      
      column(
        width = 8,
        
        div(
          class = "ora-panel",
          h4("Selection Summary"),
          verbatimTextOutput(ns("selection_summary"))
        ),
        
        div(
          class = "ora-panel",
          h4("ORA Plot"),
          plotOutput(ns("ora_plot"), height = "700px")
        ),
        
        div(
          class = "ora-panel",
          tabsetPanel(
            tabPanel(
              "ORA Results",
              br(),
              DTOutput(ns("ora_table"))
            ),
            tabPanel(
              "Selected Features",
              br(),
              DTOutput(ns("selected_features_table"))
            )
          )
        )
      )
    )
  )
}

# =========================================================
# SERVER
# =========================================================

oraModuleServer <- function(id) {
  
  moduleServer(id, function(input, output, session) {
    
    # =====================================================
    # 1) Read Uploaded Table
    # =====================================================
    
    raw_data <- reactive({
      
      req(input$upload)
      
      ora_smart_read(
        path = input$upload$datapath,
        original_name = input$upload$name,
        min_cols = 1
      )
    })
    
    # =====================================================
    # 2) Column Selectors
    # =====================================================
    
    output$column_selector <- renderUI({
      
      req(input$upload, input$input_mode)
      
      df <- raw_data()
      cols <- colnames(df)
      
      selected_id <- guess_column(
        cols,
        "gene|symbol|protein|uniprot|id|feature",
        default = cols[1]
      )
      
      if (input$input_mode == "feature_list") {
        
        selected_filter_col <- guess_column(
          cols,
          "cluster|group|class|category|term|set|direction",
          default = "__none__"
        )
        
        tagList(
          selectInput(
            session$ns("col_id"),
            "Gene/protein/feature ID column:",
            choices = cols,
            selected = selected_id
          ),
          
          selectInput(
            session$ns("feature_filter_col"),
            "Optional grouping/filter column:",
            choices = c("Use all rows" = "__none__", cols),
            selected = selected_filter_col
          )
        )
        
      } else {
        
        selected_logfc <- guess_column(
          cols,
          "log2FoldChange|logFC|log2FC|avg_log2FC|FC",
          default = cols[min(2, length(cols))]
        )
        
        selected_p <- guess_column(
          cols,
          "^p$|P.Value|pvalue|p.value|p_val|pval",
          default = cols[min(3, length(cols))]
        )
        
        selected_padj <- guess_column(
          cols,
          "padj|adj.P.Val|FDR|qvalue|q.value|adj",
          default = "__none__"
        )
        
        tagList(
          selectInput(
            session$ns("col_id"),
            "Gene/protein/feature ID column:",
            choices = cols,
            selected = selected_id
          ),
          
          selectInput(
            session$ns("col_logfc"),
            "log2FC column:",
            choices = cols,
            selected = selected_logfc
          ),
          
          selectInput(
            session$ns("col_p_nominal"),
            "Nominal P-value column:",
            choices = cols,
            selected = selected_p
          ),
          
          selectInput(
            session$ns("col_p_adj"),
            "Adjusted P-value / FDR column:",
            choices = c("Not available" = "__none__", cols),
            selected = selected_padj
          )
        )
      }
    })
    
    output$feature_group_selector <- renderUI({
      
      req(input$upload, input$input_mode)
      
      if (input$input_mode != "feature_list") {
        return(NULL)
      }
      
      req(input$feature_filter_col)
      
      if (input$feature_filter_col == "__none__") {
        return(NULL)
      }
      
      df <- raw_data()
      
      validate(
        need(
          input$feature_filter_col %in% colnames(df),
          "Selected grouping/filter column not found."
        )
      )
      
      vals <- df[[input$feature_filter_col]] %>%
        as.character() %>%
        stringr::str_trim()
      
      vals <- sort(unique(vals[!is.na(vals) & vals != ""]))
      
      validate(
        need(
          length(vals) > 0,
          "The selected grouping/filter column does not contain valid values."
        )
      )
      
      selectizeInput(
        session$ns("feature_filter_values"),
        "Values to include:",
        choices = vals,
        selected = vals[1],
        multiple = TRUE,
        options = list(
          plugins = list("remove_button")
        )
      )
    })
    
    output$feature_filtering_ui <- renderUI({
      
      req(input$input_mode)
      
      if (input$input_mode == "feature_list") {
        
        tagList(
          p(
            class = "ora-help",
            "Feature-list mode uses the selected feature IDs directly. If a grouping/filter column is selected, only rows matching the selected value(s) are used for ORA. No P-value or log2FC filtering is applied."
          )
        )
        
      } else {
        
        tagList(
          radioButtons(
            session$ns("p_source"),
            "Use for feature filtering:",
            choices = c(
              "Adjusted P-value / FDR" = "adjusted",
              "Nominal P-value" = "nominal"
            ),
            selected = "adjusted"
          ),
          
          numericInput(
            session$ns("feature_p_cutoff"),
            "Feature P-value/FDR cutoff:",
            value = 0.05,
            min = 0,
            max = 1,
            step = 0.01
          ),
          
          numericInput(
            session$ns("logfc_cutoff"),
            "Absolute log2FC cutoff:",
            value = 1,
            min = 0,
            step = 0.1
          ),
          
          radioButtons(
            session$ns("direction"),
            "Direction:",
            choices = c(
              "Both directions" = "both",
              "Upregulated only" = "up",
              "Downregulated only" = "down"
            ),
            selected = "both"
          )
        )
      }
    })
    
    # =====================================================
    # 3) MSigDB Subcategory Mapping
    # =====================================================
    
    observeEvent(input$category, {
      
      subcats <- msigdb_subcategory_choices(input$category)

      updateSelectInput(
        session,
        "subcategory",
        choices = subcats,
        selected = ""
      )

    }, ignoreInit = FALSE)
    
    get_safe_subcategory <- reactive({
      
      allowed <- input$category %in% c("C2", "C3", "C4", "C5", "C7")
      
      if (allowed && !is.null(input$subcategory) && input$subcategory != "") {
        input$subcategory
      } else {
        NULL
      }
    })
    
    # =====================================================
    # 4) Prepared Input Table
    # =====================================================
    
    prepared_data <- reactive({
      
      req(
        input$input_mode,
        input$col_id
      )
      
      df <- raw_data()
      
      validate(
        need(input$col_id %in% colnames(df), "ID column not found.")
      )
      
      if (input$input_mode == "feature_list") {
        
        out <- df %>%
          mutate(
            .feature_id = as.character(.data[[input$col_id]]) %>%
              stringr::str_trim(),
            .input_group = if (
              !is.null(input$feature_filter_col) &&
              input$feature_filter_col != "__none__" &&
              input$feature_filter_col %in% colnames(df)
            ) {
              as.character(.data[[input$feature_filter_col]]) %>%
                stringr::str_trim()
            } else {
              "All"
            },
            .log2fc = NA_real_,
            .p_nominal = NA_real_,
            .p_adjusted = NA_real_
          ) %>%
          filter(
            !is.na(.feature_id),
            .feature_id != ""
          )
        
        validate(
          need(
            nrow(out) > 0,
            "No valid feature IDs were found in the selected column."
          )
        )
        
        return(out)
      }
      
      req(
        input$col_logfc,
        input$col_p_nominal,
        input$col_p_adj
      )
      
      validate(
        need(input$col_logfc %in% colnames(df), "log2FC column not found."),
        need(input$col_p_nominal %in% colnames(df), "Nominal P-value column not found.")
      )
      
      if (input$p_source == "adjusted") {
        validate(
          need(
            input$col_p_adj != "__none__",
            "You selected Adjusted P-value/FDR filtering, but no Adjusted P-value column was selected."
          ),
          need(
            input$col_p_adj %in% colnames(df),
            "Adjusted P-value/FDR column not found."
          )
        )
      }
      
      omics_probability(as_numeric_flexible(df[[input$col_p_nominal]]), "Nominal P-values")
      if (input$col_p_adj != "__none__") {
        omics_probability(as_numeric_flexible(df[[input$col_p_adj]]), "Adjusted P-values")
      }
      fc <- as_numeric_flexible(df[[input$col_logfc]])
      validate(need(all(is.na(fc) | is.finite(fc)), "Log2FC values must be finite."))
      omics_number_parameter(input$feature_p_cutoff, "Feature P-value cutoff", 0, 1)
      omics_number_parameter(input$logfc_cutoff, "Log2FC cutoff", 0, Inf)

      omics_excluded_rows(sum(is.na(df[[input$col_id]]) | trimws(as.character(df[[input$col_id]])) == "" |
        is.na(fc) | is.na(as_numeric_flexible(df[[input$col_p_nominal]]))), "ORA input")

      out <- df %>%
        mutate(
          .feature_id = as.character(.data[[input$col_id]]) %>% stringr::str_trim(),
          .log2fc = as_numeric_flexible(.data[[input$col_logfc]]),
          .p_nominal = as_numeric_flexible(.data[[input$col_p_nominal]]),
          .p_adjusted = if (input$col_p_adj != "__none__") {
            as_numeric_flexible(.data[[input$col_p_adj]])
          } else {
            NA_real_
          },
          .input_group = case_when(
            .log2fc > 0 ~ "up",
            .log2fc < 0 ~ "down",
            TRUE ~ "unchanged"
          )
        ) %>%
        filter(
          !is.na(.feature_id),
          .feature_id != "",
          !is.na(.log2fc),
          !is.na(.p_nominal)
        )
      
      validate(
        need(
          nrow(out) > 0,
          "No valid rows after removing missing ID, log2FC, or P-value values."
        )
      )
      
      out
    })
    
    # =====================================================
    # 5) Selected Feature List
    # =====================================================
    
    selected_features <- reactive({
      
      df <- prepared_data()
      
      if (input$input_mode == "feature_list") {
        
        if (
          !is.null(input$feature_filter_col) &&
          input$feature_filter_col != "__none__"
        ) {
          
          validate(need(length(input$feature_filter_values) > 0,
            "Select at least one value for the grouping/filter column."))
          
          df <- df %>%
            filter(.input_group %in% input$feature_filter_values)
        }
        
        df <- df %>%
          distinct(.feature_id, .keep_all = TRUE)
        
        validate(
          need(
            nrow(df) > 0,
            "No selected features using the current list or grouping/filter selection."
          )
        )
        
        return(df)
      }
      
      p_vector <- if (input$p_source == "adjusted") {
        df$.p_adjusted
      } else {
        df$.p_nominal
      }
      
      df <- df %>%
        mutate(
          .p_filter = p_vector,
          .regulation = case_when(
            .log2fc > 0 ~ "up",
            .log2fc < 0 ~ "down",
            TRUE ~ "unchanged"
          )
        ) %>%
        filter(!is.na(.p_filter)) %>%
        filter(.p_filter <= input$feature_p_cutoff)
      
      if (input$direction == "up") {
        
        df <- df %>%
          filter(.log2fc >= input$logfc_cutoff)
        
      } else if (input$direction == "down") {
        
        df <- df %>%
          filter(.log2fc <= -input$logfc_cutoff)
        
      } else {
        
        df <- df %>%
          filter(abs(.log2fc) >= input$logfc_cutoff)
      }
      
      df <- df %>%
        arrange(.p_filter, desc(abs(.log2fc))) %>%
        distinct(.feature_id, .keep_all = TRUE)
      
      validate(
        need(
          nrow(df) > 0,
          "No selected features using the current P-value/FDR, log2FC, and direction thresholds."
        )
      )
      
      df
    })
    
    selected_features_display <- reactive({
      
      req(selected_features())
      
      if (input$input_mode == "feature_list") {
        
        selected_features() %>%
          transmute(
            FeatureID = .feature_id,
            Selected_group = .input_group
          )
        
      } else {
        
        selected_features() %>%
          transmute(
            FeatureID = .feature_id,
            log2FC = round(.log2fc, 4),
            `P-value used` = signif(.p_filter, 4),
            `Nominal P-value` = signif(.p_nominal, 4),
            `Adjusted P-value` = ifelse(
              is.na(.p_adjusted),
              NA,
              signif(.p_adjusted, 4)
            ),
            Direction = .regulation
          )
      }
    })
    
    # =====================================================
    # 6) Universe / Background
    # =====================================================
    
    universe_features <- reactive({
      
      if (input$universe_source == "custom") {
        
        req(input$universe_file)
        
        universe <- read_universe_file(
          path = input$universe_file$datapath,
          original_name = input$universe_file$name
        )
        
      } else {
        
        df <- prepared_data()
        universe <- clean_feature_ids(df$.feature_id)
      }
      
      validate(
        need(
          length(universe) > 0,
          "Universe/background is empty."
        )
      )
      
      universe
    })
    
    # =====================================================
    # 7) Gene Set Collection
    # =====================================================
    
    term2gene <- reactive({
      
      if (input$gene_set_source == "msigdb") {
        
        get_msigdb_terms(
          species = input$species,
          category = input$category,
          subcategory = get_safe_subcategory()
        )
        
      } else {
        
        req(input$gmt_file)
        read_gmt_terms(input$gmt_file$datapath)
      }
    })
    
    # =====================================================
    # 8) Selection Summary
    # =====================================================
    
    output$selection_summary <- renderPrint({
      
      req(input$upload)
      
      # Evaluate reactive validation before S4 generics can wrap its condition.
      prepared <- prepared_data()
      selected <- selected_features()
      universe <- universe_features()
      total_features <- length(clean_feature_ids(prepared$.feature_id))
      selected_n <- nrow(selected)
      universe_n <- length(universe)
      
      cat("Input mode:", input$input_mode, "\n")
      cat("Total features in uploaded table:", total_features, "\n")
      cat("Selected features for ORA:", selected_n, "\n")
      cat("Universe/background size:", universe_n, "\n")
      
      if (input$input_mode == "feature_list") {
        
        cat("Feature ID column:", input$col_id, "\n")
        
        if (
          !is.null(input$feature_filter_col) &&
          input$feature_filter_col != "__none__"
        ) {
          
          cat("Grouping/filter column:", input$feature_filter_col, "\n")
          cat("Selected value(s):", paste(input$feature_filter_values, collapse = ", "), "\n")
          
        } else {
          
          cat("Grouping/filter column: none; all rows used\n")
        }
        
      } else {
        
        selected_up <- sum(selected_features()$.log2fc > 0, na.rm = TRUE)
        selected_down <- sum(selected_features()$.log2fc < 0, na.rm = TRUE)
        
        cat("Selected upregulated features:", selected_up, "\n")
        cat("Selected downregulated features:", selected_down, "\n")
        cat("Filtering P-value source:", input$p_source, "\n")
        cat("Feature P-value/FDR cutoff:", input$feature_p_cutoff, "\n")
        cat("log2FC cutoff:", input$logfc_cutoff, "\n")
        cat("Direction:", input$direction, "\n")
      }
      
      cat("Gene set source:", input$gene_set_source, "\n")
      
      if (input$gene_set_source == "msigdb") {
        cat("MSigDB species:", input$species, "\n")
        cat("MSigDB category:", input$category, "\n")
        cat("MSigDB subcategory:", get_safe_subcategory() %||% "None", "\n")
      }
    })
    
    # =====================================================
    # 9) Run ORA
    # =====================================================
    
    is_running <- reactiveVal(FALSE)
    
    observe({
      shinyjs::toggleState("run_ora", !is_running())
    })
    
    ora_res_state <- reactive({
      list(selected = selected_features()$.feature_id, universe = universe_features(),
        source = input$gene_set_source, gmt = input$gmt_file, species = input$species,
        category = input$category, subcategory = get_safe_subcategory(),
        min_size = input$min_size, max_size = input$max_size, method = input$padj_method)
    })

    ora_res <- omics_event_reactive(input$run_ora, {
      omics_set_sizes(input$min_size, input$max_size)

      
      is_running(TRUE)
      on.exit(is_running(FALSE), add = TRUE)
      
      withProgress(message = "Running ORA...", value = 0, {
        
        incProgress(0.15, detail = "Preparing selected features")
        
        selected_ids <- clean_feature_ids(selected_features()$.feature_id)
        universe_ids <- universe_features()
        
        gene_ids <- intersect(selected_ids, universe_ids)
        
        validate(
          need(
            length(gene_ids) > 0,
            "No selected features overlap with the universe/background."
          )
        )
        
        incProgress(0.35, detail = "Loading gene sets")
        
        t2g <- term2gene()
        
        validate(
          need(
            nrow(t2g) > 0,
            "No gene sets available for the selected database."
          )
        )
        
        overlap_universe <- intersect(universe_ids, t2g$gene)
        overlap_selected <- intersect(gene_ids, t2g$gene)
        
        overlap_diag <- make_ora_overlap_diagnostics(
          selected_ids = gene_ids,
          universe_ids = universe_ids,
          term2gene_df = t2g
        )
        
        if (length(overlap_universe) == 0 || length(overlap_selected) == 0) {
          
          error_msg <- format_ora_mapping_error(overlap_diag)
          show_ora_error(error_msg)
          
          validate(
            need(FALSE, error_msg)
          )
        }
        
        incProgress(0.65, detail = "Performing enrichment test")
        
        enr <- tryCatch(
          {
            clusterProfiler::enricher(
              gene = gene_ids,
              universe = universe_ids,
              TERM2GENE = t2g,
              pAdjustMethod = input$padj_method,
              pvalueCutoff = 1,
              qvalueCutoff = 1,
              minGSSize = input$min_size,
              maxGSSize = input$max_size
            )
          },
          error = function(e) {
            e
          }
        )
        
        if (inherits(enr, "error")) {
          
          raw_error <- strip_ansi_codes(conditionMessage(enr))
          
          if (grepl("No gene set left after mapping and filtering", raw_error, ignore.case = TRUE)) {
            error_msg <- format_ora_mapping_error(overlap_diag)
          } else {
            error_msg <- paste0(
              "ORA failed.\n\n",
              raw_error
            )
          }
          
          show_ora_error(error_msg)
          
          validate(
            need(FALSE, error_msg)
          )
        }
        
        if (is.null(enr)) {
          
          error_msg <- format_ora_mapping_error(overlap_diag)
          show_ora_error(error_msg)
          
          validate(
            need(FALSE, error_msg)
          )
        }
        
        res <- as.data.frame(enr)
        
        if (nrow(res) == 0) {
          
          error_msg <- format_ora_mapping_error(overlap_diag)
          show_ora_error(error_msg)
          
          validate(
            need(FALSE, error_msg)
          )
        }
        
        incProgress(0.9, detail = "Formatting results")
        
        direction_label <- if (input$input_mode == "feature_list") {
          
          if (
            !is.null(input$feature_filter_col) &&
            input$feature_filter_col != "__none__"
          ) {
            
            paste(
              input$feature_filter_col,
              paste(input$feature_filter_values, collapse = "|"),
              sep = ":"
            )
            
          } else {
            
            "feature_list"
          }
          
        } else {
          
          input$direction
        }
        
        res %>%
          as_tibble() %>%
          mutate(
            GeneRatio_numeric = ratio_to_numeric(GeneRatio),
            BgRatio_numeric = ratio_to_numeric(BgRatio),
            neg_log10_pvalue = safe_neglog10(pvalue),
            neg_log10_padjust = safe_neglog10(p.adjust),
            Direction = direction_label,
            InputFeatureN = length(gene_ids),
            UniverseN = length(universe_ids)
          ) %>%
          arrange(p.adjust, pvalue)
      })
    }, state = ora_res_state, label = "ORA")
    
    # =====================================================
    # 10) Filtered ORA Results
    # =====================================================
    
    ora_table_data <- reactive({
      omics_number_parameter(input$result_p_cutoff, "Result P-value cutoff", 0, 1)
      
      req(ora_res())
      
      metric <- if (input$result_p_source == "adjusted") {
        "p.adjust"
      } else {
        "pvalue"
      }
      
      ora_res() %>%
        filter(.data[[metric]] <= input$result_p_cutoff) %>%
        arrange(.data[[metric]])
    })
    
    ora_display_table <- reactive({
      
      req(ora_table_data())
      
      ora_table_data() %>%
        mutate(
          Gene_preview = purrr::map_chr(
            geneID,
            ~ format_gene_preview(.x, n = 30)
          )
        ) %>%
        transmute(
          Term = Description,
          GeneRatio = GeneRatio,
          BgRatio = BgRatio,
          Count = Count,
          P = signif(pvalue, 3),
          `Adj. P` = signif(p.adjust, 3),
          `-Log10(Adj. P)` = round(neg_log10_padjust, 3),
          `Gene preview` = Gene_preview
        )
    })
    
    # =====================================================
    # 11) ORA Plot
    # =====================================================
    
    create_ora_plot <- reactive({
      omics_number_parameter(input$top_n, "Number of displayed pathways", 1, .Machine$integer.max, TRUE)
      
      df <- ora_table_data()
      
      validate(
        need(
          nrow(df) > 0,
          "No ORA terms pass the current result P-value/FDR cutoff."
        )
      )
      
      df <- df %>%
        slice_head(n = input$top_n) %>%
        mutate(
          Term = stringr::str_wrap(Description, width = input$wrap_width),
          Term = forcats::fct_reorder(Term, GeneRatio_numeric)
        )
      
      color_metric <- if (input$result_p_source == "adjusted") {
        "neg_log10_padjust"
      } else {
        "neg_log10_pvalue"
      }
      
      color_label <- if (input$result_p_source == "adjusted") {
        expression(-log[10]("Adjusted P-value"))
      } else {
        expression(-log[10]("Nominal P-value"))
      }
      
      if (input$plot_type == "dot") {
        
        ggplot(
          df,
          aes(
            x = GeneRatio_numeric,
            y = Term
          )
        ) +
          geom_point(
            aes(
              size = Count,
              color = .data[[color_metric]]
            ),
            alpha = 0.9
          ) +
          scale_color_gradient(
            low = input$col_low,
            high = input$col_high
          ) +
          scale_size_continuous(range = c(3, 9)) +
          theme_bw(base_size = input$font_size) +
          labs(
            title = input$plot_title,
            x = "Gene ratio",
            y = NULL,
            size = "Count",
            color = color_label
          ) +
          theme(
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = input$font_size, color = "black"),
            axis.title = element_text(size = input$font_size, color = "black"),
            plot.title = element_text(size = input$font_size + 2, face = "bold", color = "black", hjust = 0.5),
            legend.title = element_text(size = input$font_size, color = "black"),
            legend.text = element_text(size = input$font_size, color = "black")
          )
        
      } else {
        
        ggplot(
          df,
          aes(
            x = GeneRatio_numeric,
            y = Term,
            fill = .data[[color_metric]]
          )
        ) +
          geom_col() +
          scale_fill_gradient(
            low = input$col_low,
            high = input$col_high
          ) +
          theme_bw(base_size = input$font_size) +
          labs(
            title = input$plot_title,
            x = "Gene ratio",
            y = NULL,
            fill = color_label
          ) +
          theme(
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = input$font_size, color = "black"),
            axis.title = element_text(size = input$font_size, color = "black"),
            plot.title = element_text(size = input$font_size + 2, face = "bold", color = "black", hjust = 0.5),
            legend.title = element_text(size = input$font_size, color = "black"),
            legend.text = element_text(size = input$font_size, color = "black")
          )
      }
    })
    
    output$ora_plot <- renderPlot({
      create_ora_plot()
    })
    
    # =====================================================
    # 12) Tables
    # =====================================================
    
    output$ora_table <- renderDT({
      
      datatable(
        ora_display_table(),
        rownames = FALSE,
        escape = TRUE,
        class = "cell-border stripe hover compact",
        options = list(
          pageLength = 5,
          scrollX = TRUE,
          autoWidth = FALSE,
          columnDefs = list(
            list(
              targets = 0,
              width = "380px",
              className = "dt-term-cell"
            ),
            list(
              targets = 7,
              width = "420px",
              className = "dt-gene-cell"
            ),
            list(
              targets = c(1, 2, 3, 4, 5, 6),
              className = "dt-center-cell"
            )
          )
        )
      )
    })
    
    output$selected_features_table <- renderDT({
      
      datatable(
        selected_features_display(),
        rownames = FALSE,
        escape = TRUE,
        class = "cell-border stripe hover compact",
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          autoWidth = FALSE,
          columnDefs = list(
            list(
              targets = 0,
              width = "260px",
              className = "dt-term-cell"
            ),
            list(
              targets = "_all",
              className = "dt-center-cell"
            )
          )
        )
      )
    })
    
    # =====================================================
    # 13) Downloads
    # =====================================================
    
    output$download_table <- downloadHandler(
      
      filename = function() {
        paste0("ora_results_", Sys.Date(), ".txt")
      },
      
      content = function(file) {
        write.table(
          ora_display_table(),
          file = file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )
      }
    )
    
    output$download_features <- downloadHandler(
      
      filename = function() {
        paste0("ora_selected_features_", Sys.Date(), ".txt")
      },
      
      content = function(file) {
        write.table(
          selected_features_display(),
          file = file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )
      }
    )
    
    output$download_plot <- downloadHandler(
      
      filename = function() {
        paste0("ora_plot_", Sys.Date(), "_", input$exp_width, "x", input$exp_height, ".", input$exp_fmt)
      },
      
      content = function(file) {
        
        width_in <- input$exp_width * 0.393701
        height_in <- input$exp_height * 0.393701
        
        open_ora_plot_device(
          file = file,
          format = input$exp_fmt,
          width = width_in,
          height = height_in,
          dpi = input$exp_dpi
        )
        
        print(create_ora_plot())
        dev.off()
      }
    )
  })
}
