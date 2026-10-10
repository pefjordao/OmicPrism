# =========================================================
# MODULE: Comparative GSEA (EN)
# =========================================================

library(shiny)
library(shinyjs)
library(tidyverse)
library(DT)
library(svglite)
library(colourpicker)
library(msigdbr)
library(fgsea)
library(readxl)

# =========================================================
# Helper Functions
# =========================================================

gsea_compare_clean_column_names <- function(df) {
  
  colnames(df) <- colnames(df) %>%
    as.character() %>%
    stringr::str_trim()
  
  df
}

gsea_compare_try_read_delim <- function(path, delim, decimal_mark = ".") {
  
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
        gsea_compare_clean_column_names()
    },
    error = function(e) NULL
  )
}

gsea_compare_try_read_whitespace <- function(path, decimal_mark = ".") {
  
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
        gsea_compare_clean_column_names()
    },
    error = function(e) NULL
  )
}

gsea_compare_try_read_excel <- function(path) {
  
  tryCatch(
    {
      readxl::read_excel(
        path = path,
        sheet = 1,
          .name_repair = "minimal"
      ) %>%
        as.data.frame() %>%
        gsea_compare_clean_column_names()
    },
    error = function(e) NULL
  )
}

gsea_compare_score_table <- function(df) {
  
  if (is.null(df) || !is.data.frame(df) || ncol(df) <= 0) {
    return(-Inf)
  }
  
  non_empty_names <- sum(!is.na(colnames(df)) & colnames(df) != "")
  
  ncol(df) * 1000 + non_empty_names
}

gsea_compare_smart_read_unchecked <- function(path, original_name = NULL, min_cols = 2) {
  
  ext_source <- original_name %||% path
  ext <- tolower(tools::file_ext(ext_source))
  
  if (ext %in% c("xls", "xlsx")) {
    
    validate(
      need(
        requireNamespace("readxl", quietly = TRUE),
        "Package 'readxl' is required to read Excel files."
      )
    )
    
    df <- gsea_compare_try_read_excel(path)
    
    validate(
      need(
        !is.null(df) && ncol(df) >= min_cols,
        "Invalid Excel file. The first sheet must contain a valid table with at least two columns."
      )
    )
    
    return(as.data.frame(df))
  }
  
  read_attempts <- list(
    comma_dot = gsea_compare_try_read_delim(path, ",", "."),
    semicolon_dot = gsea_compare_try_read_delim(path, ";", "."),
    semicolon_comma = gsea_compare_try_read_delim(path, ";", ","),
    tab_dot = gsea_compare_try_read_delim(path, "\t", "."),
    tab_comma = gsea_compare_try_read_delim(path, "\t", ","),
    whitespace_dot = gsea_compare_try_read_whitespace(path, "."),
    whitespace_comma = gsea_compare_try_read_whitespace(path, ",")
  )
  
  scores <- vapply(
    read_attempts,
    gsea_compare_score_table,
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

gsea_compare_smart_read <- function(...) {
  omics_validate_table(gsea_compare_smart_read_unchecked(...))
}

as_numeric_flexible <- omics_numeric

open_plot_device <- function(file, fmt, width, height, dpi) {
  
  if (fmt == "png") {
    
    png(
      filename = file,
      width = width,
      height = height,
      units = "in",
      res = dpi
    )
    
  } else if (fmt == "pdf") {
    
    pdf(
      file = file,
      width = width,
      height = height
    )
    
  } else if (fmt == "tiff") {
    
    tiff(
      filename = file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      compression = "lzw"
    )
    
  } else if (fmt == "svg") {
    
    svglite::svglite(
      file = file,
      width = width,
      height = height
    )
  }
}

normalize_column_name <- function(x) {
  
  x %>%
    as.character() %>%
    stringr::str_trim() %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "")
}

match_column_by_normalized_name <- function(target_col, candidate_cols, fallback = NULL) {
  
  target_clean <- normalize_column_name(target_col)
  candidate_clean <- normalize_column_name(candidate_cols)
  
  idx <- which(candidate_clean == target_clean)
  
  if (length(idx) > 0) {
    candidate_cols[idx[1]]
  } else {
    fallback %||% candidate_cols[1]
  }
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

normalize_feature_ids <- function(x, method = "keep") {
  
  x <- x %>%
    as.character() %>%
    stringr::str_trim()
  
  if (method == "uppercase") {
    x <- stringr::str_to_upper(x)
  }
  
  if (method == "lowercase") {
    x <- stringr::str_to_lower(x)
  }
  
  x
}

safe_neglog10 <- function(x) {
  -log10(pmax(x, .Machine$double.xmin))
}

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

gsea_compare_clean_feature_ids <- function(x) {
  
  x <- x %>%
    as.character() %>%
    stringr::str_trim()
  
  unique(x[!is.na(x) & x != ""])
}

make_gsea_compare_overlap_diagnostics <- function(ranked_ids, pathways) {
  
  ranked_ids <- gsea_compare_clean_feature_ids(ranked_ids)
  gene_set_ids <- gsea_compare_clean_feature_ids(unlist(pathways, use.names = FALSE))
  
  exact_overlap <- length(intersect(ranked_ids, gene_set_ids))
  
  case_insensitive_overlap <- length(
    intersect(
      stringr::str_to_upper(ranked_ids),
      stringr::str_to_upper(gene_set_ids)
    )
  )
  
  list(
    exact_overlap = exact_overlap,
    case_insensitive_overlap = case_insensitive_overlap
  )
}

format_gsea_compare_mapping_error <- function(diag, group_label = "Ranked list") {
  
  paste0(
    "No gene set remained after ID mapping and set-size filtering.\n\n",
    "Check species/ID format or change Min/Max set size.\n\n",
    group_label,
    " exact overlap: ",
    diag$exact_overlap,
    " | Case-insensitive overlap: ",
    diag$case_insensitive_overlap
  )
}

show_gsea_compare_error <- function(message) {
  
  showModal(
    modalDialog(
      title = "Comparative GSEA failed",
      div(
        style = "white-space: pre-wrap;",
        strip_ansi_codes(message)
      ),
      easyClose = TRUE,
      footer = modalButton("OK")
    )
  )
}

make_pathway_html <- function(x) {
  
  x <- htmltools::htmlEscape(as.character(x))
  x <- stringr::str_replace_all(x, "_", "_<wbr>")
  x
}

format_numeric_table <- function(x, digits = 3) {
  
  ifelse(
    is.na(x),
    NA_character_,
    formatC(x, digits = digits, format = "fg", flag = "#")
  )
}

get_msigdb_collection <- function(species, category, subcategory = NULL) {
  
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
      nrow(msig) > 0,
      "No MSigDB gene sets returned for this selection."
    ),
    need(
      "gs_name" %in% colnames(msig),
      "MSigDB output does not contain the expected column: gs_name."
    )
  )
  
  msig
}

prepare_ranked_list <- function(df, col_id, col_stat, id_case = "keep") {
  
  validate(
    need(
      col_id %in% colnames(df),
      paste0("ID column not found: ", col_id)
    ),
    need(
      col_stat %in% colnames(df),
      paste0("Rank/statistic column not found: ", col_stat)
    )
  )
  
  ids <- df[[col_id]] %>%
    normalize_feature_ids(method = id_case)
  
  stats <- as_numeric_flexible(df[[col_stat]])
  validate(need(col_id != col_stat, "ID and rank columns must be different."))
  validate(need(all(is.na(stats) | is.finite(stats)), "Rank statistics must be finite."))
  
  omics_excluded_rows(sum(is.na(ids) | ids == "" | is.na(stats)), "Comparative GSEA ranks")

  ranks_df <- tibble(
    gene = ids,
    stat = stats
  ) %>%
    filter(
      !is.na(gene),
      gene != "",
      !is.na(stat),
      is.finite(stat)
    ) %>%
    group_by(gene) %>%
    summarise(
      stat = stat[which.max(abs(stat))],
      .groups = "drop"
    ) %>%
    arrange(desc(stat))
  
  validate(
    need(
      nrow(ranks_df) > 0,
      "No valid ranked genes after removing missing IDs or statistics."
    )
  )
  
  ranks <- ranks_df$stat
  names(ranks) <- ranks_df$gene
  
  sort(ranks, decreasing = TRUE)
}

prepare_pathways <- function(msig, id_column, id_case = "keep") {
  
  validate(
    need(
      id_column %in% colnames(msig),
      paste0("Selected MSigDB identifier column not found: ", id_column)
    )
  )
  
  msig_clean <- msig %>%
    transmute(
      pathway = gs_name,
      gene = normalize_feature_ids(.data[[id_column]], method = id_case)
    ) %>%
    filter(
      !is.na(pathway),
      pathway != "",
      !is.na(gene),
      gene != ""
    ) %>%
    distinct()
  
  pathways <- split(
    msig_clean$gene,
    msig_clean$pathway
  )
  
  pathways <- lapply(pathways, unique)
  pathways <- pathways[lengths(pathways) > 0]
  
  validate(
    need(
      length(pathways) > 0,
      "No valid pathways remained after filtering MSigDB identifiers."
    )
  )
  
  pathways
}

# =========================================================
# UI
# =========================================================

gseaCompareUI <- function(id) {
  
  ns <- NS(id)
  
  tagList(
    
    useShinyjs(),
    
    tags$head(
      tags$style(HTML("
        .gsea-compare-panel {
          background: #ffffff;
          border-radius: 8px;
          padding: 20px;
          border: 1px solid #eeeeee;
          margin-bottom: 20px;
        }

        .gsea-compare-help {
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

        table.dataTable td.gsea-pathway-cell {
          white-space: normal !important;
          word-break: break-word;
          min-width: 280px;
          max-width: 560px;
          line-height: 1.25;
        }

        table.dataTable td.gsea-leading-cell {
          white-space: normal !important;
          word-break: break-word;
          min-width: 220px;
          max-width: 420px;
          line-height: 1.25;
        }

        table.dataTable td.dt-center,
        table.dataTable th.dt-center {
          text-align: center;
        }
      "))
    ),
    
    fluidRow(
      
      # ===================================================
      # Left Sidebar
      # ===================================================
      
      column(
        width = 4,
        
        # -------------------------------------------------
        # 1) Data Input
        # -------------------------------------------------
        
        wellPanel(
          h4("1. Data Input"),
          
          fileInput(
            ns("upload_1"),
            "Upload ranked list 1 (.txt, .csv, .tsv, .xls, .xlsx)",
            accept = c(".txt", ".csv", ".tsv", ".xls", ".xlsx")
          ),
          
          fileInput(
            ns("upload_2"),
            "Upload ranked list 2 (.txt, .csv, .tsv, .xls, .xlsx)",
            accept = c(".txt", ".csv", ".tsv", ".xls", ".xlsx")
          ),
          
          h5("Group names"),
          
          textInput(
            ns("group1_name"),
            "Group 1 name:",
            value = "Group 1"
          ),
          
          textInput(
            ns("group2_name"),
            "Group 2 name:",
            value = "Group 2"
          )
        ),
        
        # -------------------------------------------------
        # 2) Column Mapping
        # -------------------------------------------------
        
        wellPanel(
          h4("2. Column Mapping"),
          
          uiOutput(ns("column_selector")),
          
          p(
            class = "gsea-compare-help",
            "Columns are selected independently for each file. Equivalent names such as gene_symbol, Gene_symbol, gene.symbol, and Gene Symbol are auto-matched when possible."
          )
        ),
        
        # -------------------------------------------------
        # 3) MSigDB Selection
        # -------------------------------------------------
        
        wellPanel(
          h4("3. MSigDB Selection"),
          
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
        
        # -------------------------------------------------
        # 4) GSEA Parameters
        # -------------------------------------------------
        
        wellPanel(
          h4("4. GSEA Parameters"),
          fluidRow(
            column(
              6,
              numericInput(
                ns("min_size"),
                "Min set size:",
                value = 15,
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
          
          numericInput(
            ns("padj_cutoff"),
            "Adjusted P-value cutoff:",
            value = 0.05,
            min = 0,
            max = 1,
            step = 0.01
          ),
          
          actionButton(
            ns("run_compare"),
            "Run comparative GSEA",
            class = "btn-primary",
            style = "width:100%; font-weight:bold;"
          )
        ),
        
        # -------------------------------------------------
        # 5) Plot Options
        # -------------------------------------------------
        
        wellPanel(
          h4("5. Plot Options"),

          textInput(
            ns("plot_title"),
            "Plot Title",
            value = "Comparative GSEA"
          ),

          selectInput(
            ns("plot_type"),
            "Plot type:",
            choices = c(
              "Bubble plot" = "bubble",
              "Bar plot" = "bar"
            ),
            selected = "bubble"
          ),
          
          selectInput(
            ns("selection_strategy"),
            "Pathway selection strategy:",
            choices = c(
              "Balanced top pathways per group" = "balanced",
              "Top pathways significant in either group" = "either",
              "Shared significant pathways only" = "shared"
            ),
            selected = "balanced"
          ),
          
          p(
            class = "gsea-compare-help",
            "Balanced selection takes significant pathways from each ranked list before displaying both groups side by side."
          ),
          
          sliderInput(
            ns("top_n"),
            "Top N pathways:",
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
            "Pathway label wrap width:",
            min = 20,
            max = 80,
            value = 45,
            step = 5
          )
        ),
        
        # -------------------------------------------------
        # 6) Colors
        # -------------------------------------------------
        
        wellPanel(
          h4("6. Colors"),
          
          colourInput(
            ns("col_neg"),
            "Negative NES",
            "#2A78D6"
          ),

          colourInput(
            ns("col_mid"),
            "Midpoint (NES = 0)",
            "#F0EFEC"
          ),

          colourInput(
            ns("col_pos"),
            "Positive NES",
            "#E34948"
          )
        ),
        
        # -------------------------------------------------
        # 7) Export
        # -------------------------------------------------
        
        wellPanel(
          h4("7. Export"),
          
          fluidRow(
            column(
              6,
              numericInput(
                ns("exp_width"),
                "Width (cm)",
                value = 26,
                min = 5,
                max = 60
              )
            ),

            column(
              6,
              numericInput(
                ns("exp_height"),
                "Height (cm)",
                value = 15,
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
                  min = 72,
                  max = 1200,
                  step = 50
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
            "Download selected pathways table",
            class = "btn-info",
            style = "width:100%"
          ),
          
          br(),
          br(),
          
          downloadButton(
            ns("download_all_results"),
            "Download all GSEA results",
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
          class = "gsea-compare-panel",
          h4("Input / Mapping Summary"),
          verbatimTextOutput(ns("input_summary"))
        ),
        
        div(
          class = "gsea-compare-panel",
          h4("Comparative GSEA Plot"),
          plotOutput(ns("compare_plot"), height = "700px")
        ),
        
        div(
          class = "gsea-compare-panel",
          h4("Selected pathways"),
          DTOutput(ns("compare_table"))
        )
      )
    )
  )
}

# =========================================================
# SERVER
# =========================================================

gseaCompareServer <- function(id) {
  
  moduleServer(id, function(input, output, session) {
    
    # =====================================================
    # 1) Read Uploaded Data
    # =====================================================
    
    raw_data_1 <- reactive({
      
      req(input$upload_1)
      
      gsea_compare_smart_read(
        path = input$upload_1$datapath,
        original_name = input$upload_1$name,
        min_cols = 2
      )
    })
    
    raw_data_2 <- reactive({
      
      req(input$upload_2)
      
      gsea_compare_smart_read(
        path = input$upload_2$datapath,
        original_name = input$upload_2$name,
        min_cols = 2
      )
    })
    
    # =====================================================
    # 2) Column Selectors
    # =====================================================
    
    output$column_selector <- renderUI({
      
      req(input$upload_1, input$upload_2)
      
      df1 <- raw_data_1()
      df2 <- raw_data_2()
      
      cols1 <- colnames(df1)
      cols2 <- colnames(df2)
      
      selected_id_1 <- guess_column(
        cols1,
        "gene|symbol|protein|uniprot|id|feature",
        default = cols1[1]
      )
      
      selected_stat_1 <- guess_column(
        cols1,
        "stat|t$|t_value|t.value|wald|score|rank|logFC|log2FC|log2FoldChange",
        default = cols1[min(2, length(cols1))]
      )
      
      selected_id_2 <- match_column_by_normalized_name(
        target_col = selected_id_1,
        candidate_cols = cols2,
        fallback = guess_column(
          cols2,
          "gene|symbol|protein|uniprot|id|feature",
          default = cols2[1]
        )
      )
      
      selected_stat_2 <- match_column_by_normalized_name(
        target_col = selected_stat_1,
        candidate_cols = cols2,
        fallback = guess_column(
          cols2,
          "stat|t$|t_value|t.value|wald|score|rank|logFC|log2FC|log2FoldChange",
          default = cols2[min(2, length(cols2))]
        )
      )
      
      tagList(
        
        h5("Ranked list 1 columns"),
        
        selectInput(
          session$ns("col_id_1"),
          "Gene/protein/feature ID column — list 1:",
          choices = cols1,
          selected = selected_id_1
        ),
        
        selectInput(
          session$ns("col_stat_1"),
          "Rank/statistic column — list 1:",
          choices = cols1,
          selected = selected_stat_1
        ),
        
        hr(),
        
        h5("Ranked list 2 columns"),
        
        selectInput(
          session$ns("col_id_2"),
          "Gene/protein/feature ID column — list 2:",
          choices = cols2,
          selected = selected_id_2
        ),
        
        selectInput(
          session$ns("col_stat_2"),
          "Rank/statistic column — list 2:",
          choices = cols2,
          selected = selected_stat_2
        )
      )
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
    # 4) Running State
    # =====================================================
    
    is_running <- reactiveVal(FALSE)
    
    observe({
      shinyjs::toggleState("run_compare", condition = !is_running())
    })
    
    # =====================================================
    # 5) Input Summary
    # =====================================================
    
    output$input_summary <- renderPrint({
      
      if (is.null(input$upload_1) || is.null(input$upload_2)) {
        cat("Upload both ranked lists to inspect column mapping.\n")
        return()
      }
      
      req(
        input$col_id_1,
        input$col_stat_1,
        input$col_id_2,
        input$col_stat_2
      )
      
      df1 <- raw_data_1()
      df2 <- raw_data_2()
      
      cat("Ranked list 1:", input$upload_1$name, "\n")
      cat("Rows:", nrow(df1), " | Columns:", ncol(df1), "\n")
      cat("ID column:", input$col_id_1, "\n")
      cat("Rank column:", input$col_stat_1, "\n\n")
      
      cat("Ranked list 2:", input$upload_2$name, "\n")
      cat("Rows:", nrow(df2), " | Columns:", ncol(df2), "\n")
      cat("ID column:", input$col_id_2, "\n")
      cat("Rank column:", input$col_stat_2, "\n\n")
      
      cat("MSigDB species:", input$species, "\n")
      cat("MSigDB category:", input$category, "\n")
      cat("MSigDB subcategory:", get_safe_subcategory() %||% "None", "\n")
      cat("GSEA algorithm: fgseaMultilevel\n")
      cat("Pathway selection strategy:", input$selection_strategy %||% "balanced", "\n")
    })
    
    # =====================================================
    # 6) Run Comparative GSEA
    # =====================================================
    
    comparison_raw_state <- reactive({
      list(data1 = raw_data_1(), data2 = raw_data_2(),
        id1 = input$col_id_1, stat1 = input$col_stat_1,
        id2 = input$col_id_2, stat2 = input$col_stat_2,
        species = input$species, category = input$category,
        subcategory = get_safe_subcategory(), min_size = input$min_size, max_size = input$max_size)
    })

    comparison_raw <- omics_event_reactive(input$run_compare, {
      omics_set_sizes(input$min_size, input$max_size)

      
      req(
        input$upload_1,
        input$upload_2,
        input$col_id_1,
        input$col_stat_1,
        input$col_id_2,
        input$col_stat_2
      )
      
      is_running(TRUE)
      on.exit(is_running(FALSE), add = TRUE)
      
      withProgress(message = "Running comparative GSEA...", value = 0, {
        
        incProgress(0.10, detail = "Checking input files")
        
        df1 <- raw_data_1()
        df2 <- raw_data_2()
        
        validate(
          need(input$col_id_1 %in% colnames(df1), "ID column not found in ranked list 1."),
          need(input$col_stat_1 %in% colnames(df1), "Rank column not found in ranked list 1."),
          need(input$col_id_2 %in% colnames(df2), "ID column not found in ranked list 2."),
          need(input$col_stat_2 %in% colnames(df2), "Rank column not found in ranked list 2.")
        )
        
        incProgress(0.25, detail = "Preparing ranked lists")
        
        ranks1 <- prepare_ranked_list(
          df = df1,
          col_id = input$col_id_1,
          col_stat = input$col_stat_1,
          id_case = "keep"
        )
        
        ranks2 <- prepare_ranked_list(
          df = df2,
          col_id = input$col_id_2,
          col_stat = input$col_stat_2,
          id_case = "keep"
        )
        
        incProgress(0.45, detail = "Loading MSigDB gene sets")
        
        msig <- get_msigdb_collection(
          species = input$species,
          category = input$category,
          subcategory = get_safe_subcategory()
        )
        
        pathways <- prepare_pathways(
          msig = msig,
          id_column = "gene_symbol",
          id_case = "keep"
        )
        
        diag_1 <- make_gsea_compare_overlap_diagnostics(
          ranked_ids = names(ranks1),
          pathways = pathways
        )
        
        diag_2 <- make_gsea_compare_overlap_diagnostics(
          ranked_ids = names(ranks2),
          pathways = pathways
        )
        
        overlap_1 <- diag_1$exact_overlap
        overlap_2 <- diag_2$exact_overlap
        
        if (overlap_1 == 0) {
          
          error_msg <- format_gsea_compare_mapping_error(
            diag = diag_1,
            group_label = input$group1_name %||% "Ranked list 1"
          )
          
          show_gsea_compare_error(error_msg)
          
          validate(
            need(FALSE, error_msg)
          )
        }
        
        if (overlap_2 == 0) {
          
          error_msg <- format_gsea_compare_mapping_error(
            diag = diag_2,
            group_label = input$group2_name %||% "Ranked list 2"
          )
          
          show_gsea_compare_error(error_msg)
          
          validate(
            need(FALSE, error_msg)
          )
        }
        
        run_gsea <- function(ranks, diag, group_label) {
          
          res <- tryCatch(
            {
              fgsea::fgseaMultilevel(
                pathways = pathways,
                stats = ranks,
                minSize = input$min_size,
                maxSize = input$max_size,
                BPPARAM = BiocParallel::SerialParam()
              ) %>%
                as_tibble()
            },
            error = function(e) {
              e
            }
          )
          
          if (inherits(res, "error")) {
            
            raw_error <- strip_ansi_codes(conditionMessage(res))
            
            if (grepl("No gene set left after mapping and filtering", raw_error, ignore.case = TRUE)) {
              error_msg <- format_gsea_compare_mapping_error(
                diag = diag,
                group_label = group_label
              )
            } else {
              error_msg <- paste0(
                "Comparative GSEA failed.\n\n",
                raw_error
              )
            }
            
            show_gsea_compare_error(error_msg)
            
            validate(
              need(FALSE, error_msg)
            )
          }
          
          if (nrow(res) == 0) {
            
            error_msg <- format_gsea_compare_mapping_error(
              diag = diag,
              group_label = group_label
            )
            
            show_gsea_compare_error(error_msg)
            
            validate(
              need(FALSE, error_msg)
            )
          }
          
          res
        }
        
        incProgress(0.65, detail = paste0("Computing ", input$group1_name))
        
        res1 <- run_gsea(
          ranks = ranks1,
          diag = diag_1,
          group_label = input$group1_name %||% "Ranked list 1"
        ) %>%
          mutate(
            InternalGroup = "G1",
            InputGenes = length(ranks1),
            OverlapWithGeneSets = overlap_1
          )
        
        incProgress(0.85, detail = paste0("Computing ", input$group2_name))
        
        res2 <- run_gsea(
          ranks = ranks2,
          diag = diag_2,
          group_label = input$group2_name %||% "Ranked list 2"
        ) %>%
          mutate(
            InternalGroup = "G2",
            InputGenes = length(ranks2),
            OverlapWithGeneSets = overlap_2
          )
        
        bind_rows(res1, res2) %>%
          mutate(
            pathway_clean = pathway,
            neg_log10_padj = safe_neglog10(padj),
            neg_log10_pval = safe_neglog10(pval),
            leadingEdge = purrr::map_chr(
              leadingEdge,
              ~ paste(.x, collapse = ", ")
            )
          ) %>%
          arrange(padj, desc(abs(NES)))
      })
    }, state = comparison_raw_state, label = "Comparative GSEA")
    
    # =====================================================
    # 7) Apply Display Names Dynamically
    # =====================================================
    
    comparison_data <- reactive({
      labels <- trimws(c(input$group1_name, input$group2_name))
      validate(need(length(labels) == 2 && all(!is.na(labels) & nzchar(labels)) &&
                    !anyDuplicated(labels), "Group names must be nonempty and different."))
      
      req(comparison_raw())
      
      comparison_raw() %>%
        mutate(
          Group = if_else(
            InternalGroup == "G1",
            trimws(input$group1_name),
            trimws(input$group2_name)
          )
        )
    })
    
    # =====================================================
    # 8) Selected Pathways
    # =====================================================
    
    selected_pathways <- reactive({
      
      req(comparison_data())
      
      df <- comparison_data()
      strategy <- input$selection_strategy %||% "balanced"
      cutoff <- input$padj_cutoff %||% 0.05
      omics_number_parameter(cutoff, "Adjusted P-value cutoff", 0, 1)
      omics_number_parameter(input$top_n, "Number of displayed pathways", 1, .Machine$integer.max, TRUE)
      
      pathway_stats <- df %>%
        group_by(pathway, pathway_clean) %>%
        summarise(
          min_padj = min(padj, na.rm = TRUE),
          max_abs_nes = max(abs(NES), na.rm = TRUE),
          n_significant_groups = sum(!is.na(padj) & is.finite(padj) & padj <= cutoff),
          .groups = "drop"
        )
      
      if (strategy == "balanced") {
        
        n_groups <- max(1, dplyr::n_distinct(df$InternalGroup))
        n_per_group <- max(1, ceiling(input$top_n / n_groups))
        
        pathways_keep <- df %>%
          filter(
            !is.na(padj),
            is.finite(padj),
            padj <= cutoff
          ) %>%
          group_by(InternalGroup) %>%
          arrange(
            padj,
            desc(abs(NES)),
            .by_group = TRUE
          ) %>%
          slice_head(n = n_per_group) %>%
          ungroup() %>%
          distinct(pathway, pathway_clean) %>%
          left_join(
            pathway_stats,
            by = c("pathway", "pathway_clean")
          ) %>%
          arrange(
            min_padj,
            desc(n_significant_groups),
            desc(max_abs_nes)
          )
        
      } else if (strategy == "shared") {
        
        expected_groups <- max(1, dplyr::n_distinct(df$InternalGroup))
        
        pathways_keep <- pathway_stats %>%
          filter(
            is.finite(min_padj),
            n_significant_groups >= expected_groups
          ) %>%
          arrange(
            min_padj,
            desc(max_abs_nes)
          ) %>%
          slice_head(n = input$top_n)
        
      } else {
        
        pathways_keep <- pathway_stats %>%
          filter(
            is.finite(min_padj),
            min_padj <= cutoff
          ) %>%
          arrange(
            min_padj,
            desc(max_abs_nes)
          ) %>%
          slice_head(n = input$top_n)
      }
      
      validate(
        need(
          nrow(pathways_keep) > 0,
          "No pathways pass the current Adjusted P-value cutoff for the selected strategy."
        )
      )
      
      pathways_keep
    })
    
    selected_df <- reactive({
      
      keep <- selected_pathways()
      
      pathway_levels <- keep %>%
        mutate(
          pathway_plot = stringr::str_wrap(
            pathway_clean,
            width = input$wrap_width
          )
        ) %>%
        pull(pathway_plot)
      
      comparison_data() %>%
        filter(pathway %in% keep$pathway) %>%
        mutate(
          pathway_plot = stringr::str_wrap(
            pathway_clean,
            width = input$wrap_width
          ),
          pathway_plot = factor(
            pathway_plot,
            levels = rev(pathway_levels)
          ),
          Group = factor(
            Group,
            levels = trimws(c(input$group1_name, input$group2_name))
          ),
          Significant = !is.na(padj) & padj <= input$padj_cutoff,
          Significance = if_else(
            Significant,
            paste0("adj. P-value <= ", input$padj_cutoff),
            paste0("adj. P-value > ", input$padj_cutoff)
          ),
          Significance_label = if_else(
            Significant,
            "*",
            ""
          ),
          NES_fill = if_else(
            Significant,
            NES,
            NA_real_
          )
        )
    })
    
    selected_table_display <- reactive({
      
      req(selected_df())
      
      selected_df() %>%
        transmute(
          Pathway = make_pathway_html(pathway),
          Group = as.character(Group),
          NES = format_numeric_table(NES, digits = 3),
          P = format_numeric_table(pval, digits = 3),
          `Adj. P` = format_numeric_table(padj, digits = 3),
          Significance = as.character(Significance),
          `Selection strategy` = input$selection_strategy %||% "balanced",
          `-Log10(Adj. P)` = format_numeric_table(neg_log10_padj, digits = 3),
          Size = size,
          `Input genes` = InputGenes,
          `Gene-set overlap` = OverlapWithGeneSets,
          `Leading edge preview` = stringr::str_trunc(leadingEdge, width = 120)
        )
    })
    
    selected_table_export <- reactive({
      
      req(selected_df())
      
      selected_df() %>%
        transmute(
          Pathway = pathway,
          Group = as.character(Group),
          NES = NES,
          `P-value` = pval,
          `Adjusted P-value` = padj,
          Significance = as.character(Significance),
          SelectionStrategy = input$selection_strategy %||% "balanced",
          `-Log10(Adjusted P-value)` = neg_log10_padj,
          Size = size,
          `Input genes` = InputGenes,
          `Gene-set overlap` = OverlapWithGeneSets,
          LeadingEdge = leadingEdge
        )
    })
    
    # =====================================================
    # 9) Comparative Plot
    # =====================================================
    
    create_plot <- reactive({
      
      df <- selected_df()
      
      validate(
        need(
          nrow(df) > 0,
          "No selected pathways available for plotting."
        )
      )
      
      if (input$plot_type == "bubble") {
        
        ggplot(
          df,
          aes(
            x = Group,
            y = pathway_plot
          )
        ) +
          geom_point(
            aes(
              size = neg_log10_padj,
              fill = NES
            ),
            shape = 21,
            color = "grey35",
            stroke = 0.45,
            alpha = 0.90
          ) +
          geom_text(
            aes(
              label = Significance_label
            ),
            color = "black",
            fontface = "bold",
            size = input$font_size * 0.45,
            vjust = 0.35,
            show.legend = FALSE
          ) +
          scale_fill_gradient2(
            low = input$col_neg,
            mid = input$col_mid,
            high = input$col_pos,
            midpoint = 0
          ) +
          scale_size_continuous(
            range = c(3, 10)
          ) +
          theme_bw(base_size = input$font_size) +
          labs(
            title = input$plot_title,
            x = NULL,
            y = NULL,
            fill = "NES",
            size = expression(-log[10]("Adjusted P-value")),
            caption = paste0("* adj. P-value <= ", input$padj_cutoff)
          ) +
          theme(
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = input$font_size, color = "black"),
            axis.title = element_text(size = input$font_size, color = "black"),
            plot.title = element_text(size = input$font_size + 2, face = "bold", color = "black", hjust = 0.5),
            plot.caption = element_text(size = input$font_size, hjust = 0, color = "black"),
            legend.title = element_text(size = input$font_size, color = "black"),
            legend.text = element_text(size = input$font_size, color = "black")
          )
        
      } else {
        
        ggplot(
          df,
          aes(
            x = NES,
            y = pathway_plot,
            fill = NES_fill
          )
        ) +
          geom_col(width = 0.75) +
          facet_wrap(~ Group) +
          geom_vline(
            xintercept = 0,
            linetype = "dashed",
            linewidth = 0.4,
            color = "grey40"
          ) +
          scale_fill_gradient2(
            low = input$col_neg,
            mid = input$col_mid,
            high = input$col_pos,
            midpoint = 0,
            na.value = "grey80"
          ) +
          theme_bw(base_size = input$font_size) +
          labs(
            title = input$plot_title,
            x = "Normalized Enrichment Score (NES)",
            y = NULL,
            fill = "NES",
            caption = paste0("Grey bars: adj. P-value > ", input$padj_cutoff)
          ) +
          theme(
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = input$font_size, color = "black"),
            axis.title = element_text(size = input$font_size, color = "black"),
            plot.title = element_text(size = input$font_size + 2, face = "bold", color = "black", hjust = 0.5),
            plot.caption = element_text(size = input$font_size, hjust = 0, color = "black"),
            strip.background = element_rect(fill = "grey95", color = "grey40"),
            strip.text = element_text(size = input$font_size, face = "bold", color = "black"),
            legend.title = element_text(size = input$font_size, color = "black"),
            legend.text = element_text(size = input$font_size, color = "black")
          )
      }
    })
    
    # =====================================================
    # 10) Render Outputs
    # =====================================================
    
    output$compare_plot <- renderPlot({
      create_plot()
    })
    
    output$compare_table <- renderDT({
      
      datatable(
        selected_table_display(),
        escape = FALSE,
        rownames = FALSE,
        filter = "top",
        class = "stripe hover compact",
        options = list(
          pageLength = 5,
          scrollX = TRUE,
          autoWidth = FALSE,
          columnDefs = list(
            list(
              targets = 0,
              className = "gsea-pathway-cell"
            ),
            list(
              targets = 11,
              className = "gsea-leading-cell"
            ),
            list(
              targets = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
              className = "dt-center"
            )
          )
        )
      )
    })
    
    # =====================================================
    # 11) Downloads
    # =====================================================
    
    output$download_plot <- downloadHandler(
      
      filename = function() {
        paste0(
          "comparative_gsea_plot_",
          Sys.Date(),
          "_",
          input$exp_width, "x", input$exp_height,
          ".",
          input$exp_fmt
        )
      },
      
      content = function(file) {
        
        width_in <- input$exp_width / 2.54
        height_in <- input$exp_height / 2.54
        
        open_plot_device(
          file = file,
          fmt = input$exp_fmt,
          width = width_in,
          height = height_in,
          dpi = input$exp_dpi
        )
        
        print(create_plot())
        dev.off()
      }
    )
    
    output$download_table <- downloadHandler(
      
      filename = function() {
        paste0(
          "comparative_gsea_selected_pathways_",
          Sys.Date(),
          ".txt"
        )
      },
      
      content = function(file) {
        
        write.table(
          selected_table_export(),
          file = file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )
      }
    )
    
    output$download_all_results <- downloadHandler(
      
      filename = function() {
        paste0(
          "comparative_gsea_all_results_",
          Sys.Date(),
          ".txt"
        )
      },
      
      content = function(file) {
        
        write.table(
          comparison_data(),
          file = file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )
      }
    )
    
    # =====================================================
    # 12) Return Reactive Objects
    # =====================================================
    
    return(
      list(
        comparison_data = comparison_data,
        selected_df = selected_df
      )
    )
  })
}
