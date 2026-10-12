# =========================================================
# MODULE: Single GSEA (EN)
# =========================================================
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

gsea_clean_column_names <- function(df) {
  
  colnames(df) <- colnames(df) %>%
    as.character() %>%
    stringr::str_trim()
  
  df
}

gsea_try_read_delim <- function(path, delim, decimal_mark = ".") {
  
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
        gsea_clean_column_names()
    },
    error = function(e) NULL
  )
}

gsea_try_read_whitespace <- function(path, decimal_mark = ".") {
  
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
        gsea_clean_column_names()
    },
    error = function(e) NULL
  )
}

gsea_try_read_excel <- function(path) {
  
  tryCatch(
    {
      readxl::read_excel(
        path = path,
        sheet = 1,
          .name_repair = "minimal"
      ) %>%
        as.data.frame() %>%
        gsea_clean_column_names()
    },
    error = function(e) NULL
  )
}

gsea_score_table <- function(df) {
  
  if (is.null(df) || !is.data.frame(df) || ncol(df) <= 0) {
    return(-Inf)
  }
  
  non_empty_names <- sum(!is.na(colnames(df)) & colnames(df) != "")
  
  ncol(df) * 1000 + non_empty_names
}

gsea_smart_read_unchecked <- function(path, original_name = NULL, min_cols = 2) {
  
  ext_source <- original_name %||% path
  ext <- tolower(tools::file_ext(ext_source))
  
  if (ext %in% c("xls", "xlsx")) {
    
    validate(
      need(
        requireNamespace("readxl", quietly = TRUE),
        "Package 'readxl' is required to read Excel files."
      )
    )
    
    df <- gsea_try_read_excel(path)
    
    validate(
      need(
        !is.null(df) && ncol(df) >= min_cols,
        "Invalid Excel file. The first sheet must contain a valid table with at least two columns."
      )
    )
    
    return(as.data.frame(df))
  }
  
  read_attempts <- list(
    comma_dot = gsea_try_read_delim(
      path = path,
      delim = ",",
      decimal_mark = "."
    ),
    
    semicolon_dot = gsea_try_read_delim(
      path = path,
      delim = ";",
      decimal_mark = "."
    ),
    
    semicolon_comma = gsea_try_read_delim(
      path = path,
      delim = ";",
      decimal_mark = ","
    ),
    
    tab_dot = gsea_try_read_delim(
      path = path,
      delim = "	",
      decimal_mark = "."
    ),
    
    tab_comma = gsea_try_read_delim(
      path = path,
      delim = "	",
      decimal_mark = ","
    ),
    
    whitespace_dot = gsea_try_read_whitespace(
      path = path,
      decimal_mark = "."
    ),
    
    whitespace_comma = gsea_try_read_whitespace(
      path = path,
      decimal_mark = ","
    )
  )
  
  scores <- vapply(
    read_attempts,
    gsea_score_table,
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

gsea_smart_read <- function(...) {
  omics_validate_table(gsea_smart_read_unchecked(...))
}

as_numeric_flexible <- omics_numeric

gsea_open_plot_device <- function(file, fmt, width, height, dpi) {
  
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

gsea_clean_feature_ids <- function(x) {
  
  x <- x %>%
    as.character() %>%
    stringr::str_trim()
  
  unique(x[!is.na(x) & x != ""])
}

make_gsea_overlap_diagnostics <- function(ranked_ids, pathways) {
  
  ranked_ids <- gsea_clean_feature_ids(ranked_ids)
  gene_set_ids <- gsea_clean_feature_ids(unlist(pathways, use.names = FALSE))
  
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

format_gsea_mapping_error <- function(diag) {
  
  paste0(
    "No gene set remained after ID mapping and set-size filtering.\n\n",
    "Check species/ID format or change Min/Max set size.\n\n",
    "Exact overlap: ",
    diag$exact_overlap,
    " | Case-insensitive overlap: ",
    diag$case_insensitive_overlap
  )
}

show_gsea_error <- function(message) {
  
  showModal(
    modalDialog(
      title = "GSEA failed",
      div(
        style = "white-space: pre-wrap;",
        strip_ansi_codes(message)
      ),
      easyClose = TRUE,
      footer = modalButton("OK")
    )
  )
}

format_leading_edge_preview <- function(x, n = 20) {
  
  if (is.null(x) || length(x) == 0) {
    return("")
  }
  
  paste(
    head(as.character(x), n),
    collapse = ", "
  )
}

get_msigdb_collection <- function(species, category, subcategory = NULL) {
  omics_get_msigdb(species, category, subcategory)
}

# =========================================================
# UI
# =========================================================

gseaModuleUI <- function(id) {
  
  ns <- NS(id)
  
  tagList(
    
    useShinyjs(),
    
    tags$head(
      tags$style(HTML("
        .gsea-panel {
          background: #ffffff;
          border-radius: 8px;
          padding: 20px;
          border: 1px solid #eeeeee;
          margin-bottom: 20px;
        }

        .gsea-help {
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

        table.dataTable td.dt-pathway-cell {
          white-space: normal !important;
          word-break: break-word;
          min-width: 320px;
          max-width: 520px;
        }

        table.dataTable td.dt-leading-cell {
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
            "Upload ranked gene/protein list (.txt, .csv, .tsv, .xls, .xlsx)",
            accept = c(".txt", ".csv", ".tsv", ".xls", ".xlsx")
          ),
          
          uiOutput(ns("column_selector")),
          
          p(
            class = "gsea-help",
            "The ranking column should be a continuous statistic, such as t-statistic, Wald statistic, signed -log10(p), or log2FC."
          ),
          
          hr(),
          
          # -------------------------------------------------
          # 2) MSigDB Selection
          # -------------------------------------------------
          
          h4("2. Gene Set Source"),
          radioButtons(
            ns("gene_set_source"), "Gene set source:",
            choices = c("MSigDB" = "msigdb", "Custom GMT file" = "gmt"),
            selected = "msigdb"
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'msigdb'", ns("gene_set_source")),
            omics_msigdb_ui(ns)
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'gmt'", ns("gene_set_source")),
            fileInput(ns("gmt_file"), "Upload custom gene sets (.gmt)", accept = c(".gmt")),
            p(class = "help-block", "Input IDs must match the GMT exactly. No gene/protein IDs are converted.")
          ),
          
          hr(),
          
          # -------------------------------------------------
          # 3) GSEA Parameters
          # -------------------------------------------------
          
          h4("3. GSEA Parameters"),
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
            ns("seed"),
            "Random seed:",
            value = 1234,
            step = 1
          ),

          numericInput(
            ns("padj_cutoff"),
            "Adjusted P-value cutoff:",
            value = 0.05,
            min = 0,
            max = 1,
            step = 0.01
          ),
          
          br(),
          
          actionButton(
            ns("run_gsea"),
            "Run GSEA",
            class = "btn-primary",
            style = "width:100%; font-weight:bold;"
          ),
          
          hr(),
          
          # -------------------------------------------------
          # 4) Plot Options
          # -------------------------------------------------
          
          h4("4. Plot Options"),

          textInput(
            ns("plot_title"),
            "Plot Title",
            value = "Single GSEA"
          ),

          selectInput(
            ns("plot_type"),
            "Plot type:",
            choices = c(
              "Bar plot" = "bar",
              "Bubble plot" = "bubble"
            ),
            selected = "bubble"
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
          ),
          
          hr(),
          
          # -------------------------------------------------
          # 5) Colors
          # -------------------------------------------------
          
          h4("5. Colors"),
          
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
          ),
          
          hr(),
          
          # -------------------------------------------------
          # 6) Export
          # -------------------------------------------------
          
          h4("6. Export"),
          
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
            "Download results table",
            class = "btn-info",
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
          class = "gsea-panel",
          h4("GSEA Plot"),
          plotOutput(ns("gsea_plot"), height = "700px")
        ),
        
        div(
          class = "gsea-panel",
          h4("Results table"),
          DTOutput(ns("gsea_table"))
        )
      )
    )
  )
}

# =========================================================
# SERVER
# =========================================================

gseaModuleServer <- function(id) {
  
  moduleServer(id, function(input, output, session) {
    
    # =====================================================
    # 1) Read Uploaded Data
    # =====================================================
    
    raw_data <- reactive({
      
      req(input$upload)
      
      gsea_smart_read(
        path = input$upload$datapath,
        original_name = input$upload$name,
        min_cols = 2
      )
    })
    
    # =====================================================
    # 2) Column Selectors
    # =====================================================
    
    output$column_selector <- renderUI({
      
      req(input$upload)
      
      df <- raw_data()
      cols <- colnames(df)
      
      selected_id <- guess_column(
        cols,
        "gene|symbol|protein|uniprot|id|feature",
        default = cols[1]
      )
      
      selected_stat <- guess_column(
        cols,
        "stat|t$|t_value|t.value|wald|score|rank|logFC|log2FC|log2FoldChange",
        default = cols[2]
      )
      
      tagList(
        selectInput(
          session$ns("col_id"),
          "Gene symbol / ID column:",
          choices = cols,
          selected = selected_id
        ),
        
        selectInput(
          session$ns("col_stat"),
          "Rank column:",
          choices = cols,
          selected = selected_stat
        )
      )
    })
    
    # =====================================================
    # 3) MSigDB Subcategory Mapping
    # =====================================================
    
    get_safe_subcategory <- omics_msigdb_server(input, session, "gene_set_source")
    
    # =====================================================
    # 4) Prepare Ranked List
    # =====================================================
    
    ranked_list <- reactive({
      
      req(input$col_id, input$col_stat)
      
      df <- raw_data()
      
      validate(
        need(input$col_id %in% colnames(df), "ID column not found."),
        need(input$col_stat %in% colnames(df), "Rank column not found.")
      )
      
      ids <- df[[input$col_id]] %>%
        as.character() %>%
        stringr::str_trim()
      
      stats <- as_numeric_flexible(df[[input$col_stat]])
      validate(need(input$col_id != input$col_stat, "ID and rank columns must be different."))
      validate(need(all(is.na(stats) | is.finite(stats)), "Rank statistics must be finite."))
      
      omics_excluded_rows(sum(is.na(ids) | ids == "" | is.na(stats)), "GSEA ranks")

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
    })
    
    # =====================================================
    # 5) Run GSEA
    # =====================================================
    
    is_running <- reactiveVal(FALSE)
    
    observe({
      shinyjs::toggleState("run_gsea", !is_running())
    })
    
    gene_sets <- reactive({
      if ((input$gene_set_source %||% "msigdb") == "gmt") {
        req(input$gmt_file)
        return(omics_read_gmt(input$gmt_file$datapath))
      }
      msig <- get_msigdb_collection(
        species = input$species, category = input$category,
        subcategory = get_safe_subcategory()
      )
      lapply(split(msig$gene_symbol, msig$gs_name), unique)
    })

    gsea_res_state <- reactive({
      list(ranks = ranked_list(), source = input$gene_set_source %||% "msigdb", gmt = input$gmt_file, species = input$species, category = input$category,
        subcategory = get_safe_subcategory(), min_size = input$min_size, max_size = input$max_size, seed = input$seed)
    })

    gsea_res <- omics_event_reactive(input$run_gsea, {
      omics_set_sizes(input$min_size, input$max_size)
      omics_number_parameter(input$seed, "Random seed", 0, .Machine$integer.max, TRUE)


      is_running(TRUE)
      on.exit(is_running(FALSE), add = TRUE)
      
      withProgress(message = "Running GSEA...", value = 0, {
        
        incProgress(0.20, detail = "Preparing ranked list")
        
        ranks <- ranked_list()
        
        incProgress(0.45, detail = "Loading gene sets")
        pathways <- gene_sets()
        
        validate(
          need(
            length(pathways) > 0,
            "No pathways available for the selected gene set source."
          )
        )
        
        omics_validate_gene_set_ids(names(ranks), pathways, "Ranked list")

        overlap_diag <- make_gsea_overlap_diagnostics(
          ranked_ids = names(ranks),
          pathways = pathways
        )
        
        if (overlap_diag$exact_overlap == 0) {
          
          error_msg <- format_gsea_mapping_error(overlap_diag)
          show_gsea_error(error_msg)
          
          validate(
            need(FALSE, error_msg)
          )
        }
        
        incProgress(0.70, detail = "Running fgseaMultilevel")

        set.seed(input$seed)

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
            error_msg <- format_gsea_mapping_error(overlap_diag)
          } else {
            error_msg <- paste0(
              "GSEA failed.\n\n",
              raw_error
            )
          }
          
          show_gsea_error(error_msg)
          
          validate(
            need(FALSE, error_msg)
          )
        }
        
        if (nrow(res) == 0) {
          
          error_msg <- format_gsea_mapping_error(overlap_diag)
          show_gsea_error(error_msg)
          
          validate(
            need(FALSE, error_msg)
          )
        }
        
        incProgress(0.90, detail = "Formatting results")
        
        res %>%
          mutate(
            neg_log10_padj = safe_neglog10(padj),
            neg_log10_pval = safe_neglog10(pval)
          ) %>%
          arrange(padj, desc(abs(NES)))
      })
    }, state = gsea_res_state, label = "GSEA")
    
    # =====================================================
    # 6) Significant / Filtered Results
    # =====================================================
    
    gsea_table_data <- reactive({
      
      req(gsea_res())
      omics_number_parameter(input$padj_cutoff, "Adjusted P-value cutoff", 0, 1)
      
      gsea_res() %>%
        filter(padj <= input$padj_cutoff) %>%
        arrange(padj, desc(abs(NES)))
    })
    
    gsea_display_table <- reactive({
      
      req(gsea_table_data())
      
      gsea_table_data() %>%
        mutate(
          leadingEdge_preview = purrr::map_chr(
            leadingEdge,
            ~ format_leading_edge_preview(.x, n = 20)
          )
        ) %>%
        transmute(
          Pathway = pathway,
          NES = round(NES, 3),
          ES = round(ES, 3),
          P = signif(pval, 3),
          `Adj. P` = signif(padj, 3),
          `-Log10(Adj. P)` = round(neg_log10_padj, 3),
          Size = size,
          `Leading edge preview` = leadingEdge_preview
        )
    })
    
    # Full-precision export is separate from the compact on-screen preview.
    gsea_export_table <- reactive({
      gsea_table_data() %>%
        transmute(Pathway = pathway, NES = NES, ES = ES, P = pval,
          `Adj. P` = padj, `-Log10(Adj. P)` = neg_log10_padj, Size = size,
          `Leading edge` = purrr::map_chr(leadingEdge, ~ paste(.x, collapse = ", ")))
    })

    # =====================================================
    # 7) GSEA Plot
    # =====================================================
    
    create_gsea_plot <- reactive({
      omics_number_parameter(input$top_n, "Number of displayed pathways", 1, .Machine$integer.max, TRUE)
      
      df <- gsea_table_data()
      
      validate(
        need(
          nrow(df) > 0,
          "No pathways pass the current Adjusted P-value cutoff."
        )
      )
      
      df <- df %>%
        slice_head(n = input$top_n) %>%
        mutate(
          pathway_plot = pathway %>%
            stringr::str_replace_all("_", "_ ") %>%
            stringr::str_wrap(width = input$wrap_width),
          pathway_plot = forcats::fct_reorder(
            pathway_plot,
            NES
          )
        )
      
      if (input$plot_type == "bar") {
        
        ggplot(
          df,
          aes(
            x = NES,
            y = pathway_plot,
            fill = NES
          )
        ) +
          geom_col(width = 0.75) +
          scale_fill_gradient2(
            low = input$col_neg,
            mid = input$col_mid,
            high = input$col_pos,
            midpoint = 0,
            limits = c(-1, 1) * max(abs(df$NES), .Machine$double.eps)
          ) +
          theme_bw(base_size = input$font_size) +
          labs(
            title = input$plot_title,
            x = "Normalized Enrichment Score (NES)",
            y = NULL,
            fill = "NES"
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
            x = NES,
            y = pathway_plot
          )
        ) +
          geom_point(
            aes(
              size = neg_log10_padj,
              color = NES
            ),
            alpha = 0.90
          ) +
          geom_vline(
            xintercept = 0,
            linetype = "dashed",
            linewidth = 0.4,
            color = "grey40"
          ) +
          scale_color_gradient2(
            low = input$col_neg,
            mid = input$col_mid,
            high = input$col_pos,
            midpoint = 0,
            limits = c(-1, 1) * max(abs(df$NES), .Machine$double.eps)
          ) +
          scale_size_continuous(
            range = c(3, 10)
          ) +
          theme_bw(base_size = input$font_size) +
          labs(
            title = input$plot_title,
            x = "Normalized Enrichment Score (NES)",
            y = NULL,
            color = "NES",
            size = expression(-log[10]("Adjusted P-value"))
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
    
    # =====================================================
    # 8) Render Plot
    # =====================================================
    
    output$gsea_plot <- renderPlot({
      create_gsea_plot()
    })
    
    # =====================================================
    # 9) Results Table
    # =====================================================
    
    output$gsea_table <- renderDT({
      
      datatable(
        gsea_display_table(),
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
              className = "dt-pathway-cell"
            ),
            list(
              targets = 7,
              width = "420px",
              className = "dt-leading-cell"
            ),
            list(
              targets = c(1, 2, 3, 4, 5, 6),
              className = "dt-center-cell"
            )
          )
        )
      )
    })
    
    # =====================================================
    # 10) Download Results Table
    # =====================================================
    
    output$download_table <- downloadHandler(
      
      filename = function() {
        paste0(
          "single_gsea_results_",
          Sys.Date(),
          ".txt"
        )
      },
      
      content = function(file) {
        
        write.table(
          gsea_export_table(),
          file = file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )
      }
    )
    
    # =====================================================
    # 11) Download Plot
    # =====================================================
    
    output$download_plot <- downloadHandler(
      
      filename = function() {
        paste0(
          "single_gsea_plot_",
          Sys.Date(),
          "_",
          input$exp_width, "x", input$exp_height,
          ".",
          input$exp_fmt
        )
      },
      
      content = function(file) {
        
        width_in <- input$exp_width * 0.393701
        height_in <- input$exp_height * 0.393701
        
        if (input$exp_fmt == "png") {
          
          png(
            filename = file,
            width = width_in,
            height = height_in,
            units = "in",
            res = input$exp_dpi
          )
          
        } else if (input$exp_fmt == "pdf") {
          
          pdf(
            file = file,
            width = width_in,
            height = height_in
          )
          
        } else if (input$exp_fmt == "tiff") {
          
          tiff(
            filename = file,
            width = width_in,
            height = height_in,
            units = "in",
            res = input$exp_dpi,
            compression = "lzw"
          )
          
        } else if (input$exp_fmt == "svg") {
          
          svglite::svglite(
            file = file,
            width = width_in,
            height = height_in
          )
        }
        
        print(create_gsea_plot())
        dev.off()
      }
    )
  })
}
