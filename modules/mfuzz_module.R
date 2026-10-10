# =========================================================
# EzOmicsViz — Mfuzz Module
# =========================================================

library(shiny)
library(shinyjs)
library(tidyverse)
library(DT)
library(Mfuzz)
library(Biobase)
library(svglite)
library(colourpicker)
library(cowplot)
library(readxl)

# =========================================================
# Helper Functions
# =========================================================

mfuzz_membership_defaults <- c(
  low = "#FFEA46",
  mid = "#7C7B78",
  high = "#00204D"
)

mfuzz_membership_colours <- function(low, mid, high) {
  selected <- unname(c(low, mid, high))
  if (identical(toupper(selected), unname(mfuzz_membership_defaults))) {
    return(viridisLite::cividis(256, direction = -1))
  }
  selected
}

mfuzz_clean_column_names <- function(df) {
  
  colnames(df) <- colnames(df) %>%
    as.character() %>%
    stringr::str_trim()
  
  df
}

mfuzz_smart_read_unchecked <- function(path, original_name = NULL) {
  
  ext_source <- original_name %||% path
  ext <- tolower(tools::file_ext(ext_source))
  
  if (ext %in% c("xlsx", "xls")) {
    
    validate(
      need(
        requireNamespace("readxl", quietly = TRUE),
        "Package 'readxl' is required to read Excel files."
      )
    )
    
    df <- tryCatch(
      {
        readxl::read_excel(
          path = path,
          sheet = 1,
          .name_repair = "minimal"
        ) %>%
          as.data.frame() %>%
          mfuzz_clean_column_names()
      },
      error = function(e) NULL
    )
    
    validate(
      need(
        !is.null(df) && ncol(df) > 1,
        "Invalid Excel file. The first sheet must contain a valid table with at least two columns."
      )
    )
    
    return(df)
  }
  
  df <- tryCatch(
    read.csv(
      file = path,
      header = TRUE,
      check.names = FALSE,
        colClasses = "character"
    ) %>%
      mfuzz_clean_column_names(),
    error = function(e) NULL
  )
  
  if (is.null(df) || ncol(df) <= 1) {
    
    df <- tryCatch(
      read.delim(
        file = path,
        sep = ";",
        header = TRUE,
        check.names = FALSE,
        colClasses = "character"
      ) %>%
        mfuzz_clean_column_names(),
      error = function(e) NULL
    )
  }
  
  if (is.null(df) || ncol(df) <= 1) {
    
    df <- tryCatch(
      read.table(
        file = path,
        sep = "\t",
        header = TRUE,
        check.names = FALSE,
        colClasses = "character"
      ) %>%
        mfuzz_clean_column_names(),
      error = function(e) NULL
    )
  }
  
  if (is.null(df) || ncol(df) <= 1) {
    
    df <- tryCatch(
      read.table(
        file = path,
        header = TRUE,
        check.names = FALSE,
        colClasses = "character"
      ) %>%
        mfuzz_clean_column_names(),
      error = function(e) NULL
    )
  }
  
  validate(
    need(
      !is.null(df) && ncol(df) > 1,
      "Invalid input file. Upload a valid .txt, .csv, .tsv, .xls, or .xlsx file."
    )
  )
  
  as.data.frame(df)
}

mfuzz_smart_read <- function(...) {
  omics_validate_table(mfuzz_smart_read_unchecked(...))
}

detect_replicate_groups <- function(sample_names) {
  sub("_(\\d+)$", "", sample_names)
}

row_var_base <- function(mat) {
  apply(mat, 1, var, na.rm = TRUE)
}

standardize_rows <- function(mat) {
  
  scaled <- t(scale(t(mat)))
  scaled[!is.finite(scaled)] <- 0
  
  scaled
}

make_eset <- function(mat) {
  Biobase::ExpressionSet(assayData = as.matrix(mat))
}

collapse_duplicate_features_highest_mean <- function(expr_df, feature_ids) {
  
  feature_ids <- as.character(feature_ids) %>%
    stringr::str_trim()
  
  sample_cols <- colnames(expr_df)
  
  tmp <- expr_df %>%
    mutate(.FeatureID = feature_ids) %>%
    filter(
      !is.na(.FeatureID),
      .FeatureID != ""
    )
  
  validate(
    need(
      nrow(tmp) > 0,
      "No valid feature IDs found after removing empty or missing IDs."
    )
  )
  
  tmp_mat <- as.matrix(tmp[, sample_cols, drop = FALSE])
  storage.mode(tmp_mat) <- "numeric"
  
  row_mean <- rowMeans(tmp_mat, na.rm = TRUE)
  row_mean[!is.finite(row_mean)] <- -Inf
  
  tmp$.row_mean <- row_mean
  
  tmp <- tmp %>%
    group_by(.FeatureID) %>%
    slice_max(
      order_by = .row_mean,
      n = 1,
      with_ties = FALSE
    ) %>%
    ungroup() %>%
    select(-.row_mean)
  
  mat <- as.matrix(tmp[, sample_cols, drop = FALSE])
  rownames(mat) <- tmp$.FeatureID
  storage.mode(mat) <- "numeric"
  
  mat
}

# =========================================================
# UI
# =========================================================

mfuzzModuleUI <- function(id) {
  
  ns <- NS(id)
  
  tagList(
    
    useShinyjs(),
    
    tags$head(
      tags$style(HTML("
        .mfuzz-panel {
          background: #ffffff;
          border-radius: 8px;
          padding: 20px;
          border: 1px solid #eeeeee;
          margin-bottom: 20px;
        }

        .mfuzz-help {
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
      "))
    ),
    
    fluidRow(
      
      column(
        width = 4,
        
        wellPanel(
          
          h4("1. Data Input"),
          
          fileInput(
            ns("expr_file"),
            "Upload log2 expression matrix (.txt, .csv, .tsv, .xls, .xlsx)",
            accept = c(".txt", ".csv", ".tsv", ".xls", ".xlsx")
          ),
          
          checkboxInput(
            ns("transpose"),
            "Matrix is transposed (samples in rows)",
            value = FALSE
          ),
          
          radioButtons(
            ns("replicate_mode"),
            "Replicate handling:",
            choices = c(
              "Auto-detect replicates from sample names" = "auto",
              "Samples already averaged by condition" = "mean"
            ),
            selected = "auto"
          ),
          
          uiOutput(ns("gene_column_ui")),
          uiOutput(ns("condition_order_ui")),
          
          p(
            class = "mfuzz-help",
            "Input values must already be preprocessed and log2-transformed when appropriate. The module does not apply log transformation or imputation. Duplicated feature identifiers are handled automatically by keeping the row with the highest mean abundance across samples."
          ),
          
          p(
            class = "mfuzz-help",
            "For automatic replicate detection, sample names should end with replicate suffixes such as Control_1, Control_2, Treatment_1."
          ),
          
          hr(),
          
          h4("2. Pre-processing"),
          
          sliderInput(
            ns("var_quantile"),
            "Variance filter — keep features above quantile:",
            min = 0,
            max = 0.95,
            value = 0.25,
            step = 0.05
          ),
          
          checkboxInput(
            ns("standardize"),
            "Standardize each feature before Mfuzz",
            value = TRUE
          ),
          
          p(
            class = "mfuzz-help",
            "Mfuzz usually works better with standardized rows because it clusters temporal/profile patterns rather than absolute abundance."
          ),
          
          hr(),
          
          h4("3. Mfuzz Parameters"),
          
          actionButton(
            ns("estimate_m_btn"),
            "Estimate m",
            class = "btn-info",
            style = "width:100%;"
          ),
          
          br(),
          br(),
          
          numericInput(
            ns("m_value"),
            "Fuzzifier m:",
            value = 1.25,
            min = 1.01,
            step = 0.05
          ),

          p(
            class = "mfuzz-help",
            "The fuzzifier m controls how soft (fuzzy) the Mfuzz clustering is. Values close to 1 approach hard/crisp clustering, where each feature is assigned to a single cluster; higher values produce softer, more overlapping cluster memberships. Use 'Estimate m' for a data-driven default."
          ),

          numericInput(
            ns("seed"),
            "Random seed:",
            value = 1234,
            step = 1
          ),

          sliderInput(
            ns("dmin_max_clusters"),
            "Max clusters to inspect in Dmin:",
            min = 2,
            max = 30,
            value = 10,
            step = 1
          ),
          
          actionButton(
            ns("run_dmin_btn"),
            "Run Dmin",
            class = "btn-warning",
            style = "width:100%;"
          ),
          
          br(),
          br(),
          
          numericInput(
            ns("n_clusters"),
            "Number of clusters:",
            value = 6,
            min = 2,
            max = 20,
            step = 1
          ),
          
          sliderInput(
            ns("membership_cutoff"),
            "Membership cutoff for cluster table:",
            min = 0,
            max = 1,
            value = 0.5,
            step = 0.05
          ),
          
          actionButton(
            ns("run_mfuzz_btn"),
            "Run Mfuzz clustering",
            class = "btn-primary",
            style = "width:100%; font-weight:bold;"
          ),
          
          hr(),
          
          h4("4. Plot Colors"),
          
          h5("Cluster profile"),
          
          colourInput(
            ns("profile_col_low"),
            "Low membership",
            mfuzz_membership_defaults[["low"]]
          ),

          colourInput(
            ns("profile_col_mid"),
            "Mid membership",
            mfuzz_membership_defaults[["mid"]]
          ),

          colourInput(
            ns("profile_col_high"),
            "High membership",
            mfuzz_membership_defaults[["high"]]
          ),
          
          colourInput(
            ns("profile_mean_col"),
            "Mean profile line",
            "#000000"
          ),
          
          br(),
          
          h5("Heatmap"),
          
          colourInput(
            ns("heatmap_col_low"),
            "Low expression",
            "#2A78D6"
          ),

          colourInput(
            ns("heatmap_col_mid"),
            "Mid expression",
            "#F0EFEC"
          ),

          colourInput(
            ns("heatmap_col_high"),
            "High expression",
            "#E34948"
          ),
          
          hr(),
          
          h4("5. Plot Export"),

          textInput(
            ns("plot_title"),
            "Plot Title",
            value = ""
          ),

          sliderInput(
            ns("font_size"),
            "Font size:",
            min = 8,
            max = 18,
            value = 12,
            step = 1
          ),

          fluidRow(
            column(
              6,
              numericInput(
                ns("plot_width"),
                "Width (cm)",
                value = 24,
                min = 5,
                max = 60
              )
            ),

            column(
              6,
              numericInput(
                ns("plot_height"),
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
                ns("plot_format"),
                "Format",
                choices = c("pdf", "svg", "png", "tiff"),
                selected = "svg"
              )
            ),

            column(
              6,
              conditionalPanel(
                condition = sprintf("input['%s'] == 'png' || input['%s'] == 'tiff'", ns("plot_format"), ns("plot_format")),
                numericInput(
                  ns("plot_dpi"),
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
            ns("download_combined_plot"),
            "Download combined plot",
            class = "btn-success",
            style = "width:100%"
          ),
          
          hr(),
          
          h4("6. Export Tables"),
          
          downloadButton(
            ns("download_clusters"),
            "Download cluster assignments",
            class = "btn-success",
            style = "width:100%"
          ),
          
          br(),
          br(),
          
          downloadButton(
            ns("download_membership"),
            "Download membership matrix",
            class = "btn-success",
            style = "width:100%"
          ),
          
          br(),
          br(),
          
          downloadButton(
            ns("download_processed"),
            "Download processed matrix",
            class = "btn-success",
            style = "width:100%"
          )
        )
      ),
      
      column(
        width = 8,
        
        div(
          class = "mfuzz-panel",
          h4("Matrix Summary"),
          verbatimTextOutput(ns("matrix_summary"))
        ),
        
        fluidRow(
          column(
            width = 6,
            div(
              class = "mfuzz-panel",
              h4("Estimated m"),
              verbatimTextOutput(ns("m_estimate_text"))
            )
          ),
          
          column(
            width = 6,
            div(
              class = "mfuzz-panel",
              h4("Clustering Summary"),
              verbatimTextOutput(ns("cluster_summary"))
            )
          )
        ),
        
        div(
          class = "mfuzz-panel",
          h4("Dmin Plot"),
          plotOutput(ns("dmin_plot"), height = "350px")
        ),
        
        div(
          class = "mfuzz-panel",
          h4("Cluster Profiles + Heatmap"),
          uiOutput(ns("combined_cluster_panel_ui"))
        ),
        
        div(
          class = "mfuzz-panel",
          h4("Cluster Assignment Table"),
          DTOutput(ns("cluster_table"))
        )
      )
    )
  )
}

# =========================================================
# SERVER
# =========================================================

mfuzzModuleServer <- function(id) {
  
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    # =====================================================
    # 1) Read Raw Data
    # =====================================================
    
    raw_df <- reactive({
      
      req(input$expr_file)
      
      mfuzz_smart_read(
        path = input$expr_file$datapath,
        original_name = input$expr_file$name
      )
    })
    
    # =====================================================
    # 2) Feature Column Selector
    # =====================================================
    
    output$gene_column_ui <- renderUI({
      
      req(raw_df())
      
      cols <- colnames(raw_df())
      
      selectInput(
        session$ns("gene_col"),
        "Feature ID column:",
        choices = cols,
        selected = cols[1]
      )
    })
    
    # =====================================================
    # 3) Initial Matrix
    # =====================================================
    
    initial_matrix <- reactive({
      req(raw_df(), input$gene_col)
      df <- raw_df()
      validate(need(input$gene_col %in% colnames(df), "Selected feature column not found."))
      if (isTRUE(input$transpose)) {
        sample_ids <- omics_validate_ids(df[[input$gene_col]])
        values <- omics_numeric_frame(df[, setdiff(names(df), input$gene_col), drop = FALSE])
        mat <- t(as.matrix(values))
        colnames(mat) <- sample_ids
        # Keep the selected ID column valid; no selector mutation inside a reactive.
        out <- data.frame(rownames(mat), as.data.frame(mat, check.names = FALSE),
                          check.names = FALSE, stringsAsFactors = FALSE)
        colnames(out)[1] <- input$gene_col
        return(omics_validate_table(out))
      }
      df
    })
    
    # =====================================================
    # 4) Condition Order UI
    # =====================================================
    
    output$condition_order_ui <- renderUI({
      
      req(initial_matrix(), input$gene_col)
      
      df <- initial_matrix()
      sample_cols <- setdiff(colnames(df), input$gene_col)
      
      if (input$replicate_mode == "auto") {
        choices <- unique(detect_replicate_groups(sample_cols))
      } else {
        choices <- sample_cols
      }
      
      selectizeInput(
        session$ns("condition_order"),
        "Order of conditions:",
        choices = choices,
        selected = choices,
        multiple = TRUE,
        options = list(
          plugins = list("remove_button")
        )
      )
    })
    
    # =====================================================
    # 5) Duplicate Summary
    # =====================================================
    
    duplicate_summary <- reactive({
      
      req(initial_matrix(), input$gene_col)
      
      df <- initial_matrix()
      
      ids <- as.character(df[[input$gene_col]]) %>%
        stringr::str_trim()
      
      ids <- ids[!is.na(ids) & ids != ""]
      
      tibble(
        total_ids = length(ids),
        unique_ids = length(unique(ids)),
        duplicated_rows = sum(duplicated(ids)),
        duplicated_feature_names = sum(table(ids) > 1)
      )
    })
    
    # =====================================================
    # 6) Processed Matrix
    # =====================================================
    
    processed_matrix <- reactive({
      
      req(
        initial_matrix(),
        input$gene_col,
        input$condition_order
      )
      
      df <- initial_matrix()
      
      validate(
        need(
          input$gene_col %in% colnames(df),
          "Feature column not found in matrix."
        )
      )
      
      feature_ids <- df[[input$gene_col]]
      
      expr_df <- df[, setdiff(colnames(df), input$gene_col), drop = FALSE]
      
      expr_df <- omics_numeric_frame(expr_df)
      
      validate(
        need(
          ncol(expr_df) >= 1,
          "No numeric sample columns were found."
        )
      )
      
      # ---------------------------------------------------
      # 6.1) Duplicated feature handling
      # ---------------------------------------------------
      
      mat <- collapse_duplicate_features_highest_mean(
        expr_df = expr_df,
        feature_ids = feature_ids
      )
      
      storage.mode(mat) <- "numeric"
      
      validate(
        need(
          all(is.finite(mat)),
          "The input matrix contains missing or non-finite values. Please upload a complete preprocessed/log2 matrix with finite numeric values."
        )
      )
      
      # ---------------------------------------------------
      # 6.2) Replicate handling
      # ---------------------------------------------------
      
      if (input$replicate_mode == "auto") {
        
        groups <- detect_replicate_groups(colnames(mat))
        uniq_groups <- unique(groups)
        
        mat_avg <- sapply(uniq_groups, function(g) {
          rowMeans(mat[, groups == g, drop = FALSE], na.rm = TRUE)
        })
        
        mat <- as.matrix(mat_avg)
        rownames(mat) <- rownames(mat_avg)
      }
      
      mat <- as.matrix(mat)
      storage.mode(mat) <- "numeric"
      
      # ---------------------------------------------------
      # 6.3) Condition ordering
      # ---------------------------------------------------
      
      validate(
        need(
          length(input$condition_order) >= 3,
          "Mfuzz requires at least 3 ordered conditions."
        )
      )
      
      validate(
        need(
          all(input$condition_order %in% colnames(mat)),
          "Some selected conditions are not present in the matrix."
        )
      )
      
      validate(need(!anyDuplicated(input$condition_order), "Condition order must contain unique conditions."))
      mat <- mat[, input$condition_order, drop = FALSE]
      
      validate(
        need(
          all(is.finite(mat)),
          "The ordered matrix contains missing or non-finite values."
        )
      )
      
      # ---------------------------------------------------
      # 6.4) Variance filtering
      # ---------------------------------------------------
      
      rv <- row_var_base(mat)
      rv[!is.finite(rv)] <- NA_real_
      
      validate(
        need(
          any(!is.na(rv) & rv > 0),
          "No variable features were found. Mfuzz requires features with variation across conditions."
        )
      )
      
      omics_number_parameter(input$var_quantile, "Variance quantile", 0, 1)
      cutoff <- suppressWarnings(
        stats::quantile(
          rv,
          probs = input$var_quantile,
          na.rm = TRUE
        )
      )
      
      mat <- mat[!is.na(rv) & rv > cutoff, , drop = FALSE]
      
      validate(
        need(
          nrow(mat) > 10,
          "Too few features remain after variance filtering."
        )
      )
      
      # ---------------------------------------------------
      # 6.5) Row standardization
      # ---------------------------------------------------
      
      if (isTRUE(input$standardize)) {
        mat <- standardize_rows(mat)
      }
      
      validate(
        need(
          all(is.finite(mat)),
          "The processed matrix contains missing or non-finite values after standardization."
        )
      )
      
      validate(
        need(
          nrow(mat) > 10,
          "Too few features remain after filtering."
        ),
        need(
          ncol(mat) >= 3,
          "Mfuzz requires at least 3 ordered conditions."
        )
      )
      
      mat
    })
    
    # =====================================================
    # 7) Matrix Summary
    # =====================================================
    
    output$matrix_summary <- renderPrint({
      
      req(processed_matrix(), duplicate_summary())
      
      mat <- processed_matrix()
      dup <- duplicate_summary()
      
      cat("Input matrix requirement: preprocessed/log2 matrix provided by user\n")
      cat("Duplicated feature handling: keep row with highest mean across samples\n\n")
      
      cat("Raw feature IDs:", dup$total_ids, "\n")
      cat("Unique raw feature IDs:", dup$unique_ids, "\n")
      cat("Duplicated rows:", dup$duplicated_rows, "\n")
      cat("Feature names with duplicates:", dup$duplicated_feature_names, "\n\n")
      
      cat("Processed features:", nrow(mat), "\n")
      cat("Conditions:", ncol(mat), "\n")
      cat("Replicate handling:", input$replicate_mode, "\n")
      cat("Variance quantile filter:", input$var_quantile, "\n")
      cat("Row standardization:", input$standardize, "\n\n")
      
      cat("Condition order:\n")
      cat(paste(colnames(mat), collapse = " -> "), "\n\n")
      
      cat("Processed matrix range:\n")
      print(range(mat, na.rm = TRUE))
    })
    
    # =====================================================
    # 8) Estimate m
    # =====================================================
    
    m_estimated_state <- reactive({
      processed_matrix()
    })

    m_estimated <- omics_event_reactive(input$estimate_m_btn, {

      
      req(processed_matrix())
      
      eset <- make_eset(processed_matrix())
      
      withProgress(message = "Estimating m...", value = 0.5, {
        Mfuzz::mestimate(eset)
      })
    }, state = m_estimated_state, label = "Estimate m")
    
    observeEvent(m_estimated(), {
      
      updateNumericInput(
        session,
        "m_value",
        value = round(m_estimated(), 3)
      )
    })
    
    output$m_estimate_text <- renderPrint({
      
      if (is.null(m_estimated())) {
        
        cat("Click 'Estimate m' to compute the fuzzifier.\n")
        
      } else {
        
        cat("Estimated m:", round(m_estimated(), 4), "\n")
      }
    })
    
    # =====================================================
    # 9) Dmin
    # =====================================================
    
    dmin_res_state <- reactive({
      list(matrix = processed_matrix(), m = input$m_value, seed = input$seed, max_clusters = input$dmin_max_clusters)
    })

    dmin_res <- omics_event_reactive(input$run_dmin_btn, {
      omics_number_parameter(input$m_value, "Fuzzifier m", 1, Inf)
      validate(need(input$m_value > 1, "Fuzzifier m must be greater than 1."))
      omics_number_parameter(input$seed, "Random seed", 0, .Machine$integer.max, TRUE)
      omics_number_parameter(input$dmin_max_clusters, "Maximum clusters", 2, nrow(processed_matrix()), TRUE)

      
      req(processed_matrix(), input$m_value)
      
      eset <- make_eset(processed_matrix())
      crange <- seq(2, input$dmin_max_clusters)
      
      withProgress(message = "Running Dmin...", value = 0.5, {

        set.seed(input$seed)
        
        dmin_vals <- Mfuzz::Dmin(
          eset,
          m = input$m_value,
          crange = crange,
          repeats = 3,
          visu = FALSE
        )
      })
      
      tibble(
        clusters = crange,
        Dmin = dmin_vals
      )
    }, state = dmin_res_state, label = "Dmin")
    
    output$dmin_plot <- renderPlot({
      
      req(dmin_res())
      
      ggplot(
        dmin_res(),
        aes(
          x = clusters,
          y = Dmin
        )
      ) +
        geom_line() +
        geom_point() +
        theme_bw(base_size = input$font_size) +
        theme(
          axis.text = element_text(size = input$font_size, color = "black"),
          axis.title = element_text(size = input$font_size, color = "black"),
          plot.title = element_text(size = input$font_size + 2, face = "bold", color = "black", hjust = 0.5)
        ) +
        labs(
          title = "Dmin across number of clusters",
          x = "Number of clusters",
          y = "Dmin"
        )
    })
    
    # =====================================================
    # 10) Run Mfuzz
    # =====================================================
    
    mfuzz_res_state <- reactive({
      list(matrix = processed_matrix(), m = input$m_value, seed = input$seed, clusters = input$n_clusters)
    })

    mfuzz_res <- omics_event_reactive(input$run_mfuzz_btn, {
      omics_number_parameter(input$m_value, "Fuzzifier m", 1, Inf)
      validate(need(input$m_value > 1, "Fuzzifier m must be greater than 1."))
      omics_number_parameter(input$seed, "Random seed", 0, .Machine$integer.max, TRUE)
      omics_number_parameter(input$n_clusters, "Number of clusters", 2, nrow(processed_matrix()), TRUE)

      
      req(
        processed_matrix(),
        input$m_value,
        input$n_clusters
      )
      
      validate(
        need(
          input$n_clusters >= 2,
          "Number of clusters must be >= 2."
        )
      )
      
      eset <- make_eset(processed_matrix())
      
      withProgress(message = "Running Mfuzz clustering...", value = 0.5, {

        set.seed(input$seed)
        
        cl <- Mfuzz::mfuzz(
          eset,
          c = input$n_clusters,
          m = input$m_value
        )
      })
      
      cl
    }, state = mfuzz_res_state, label = "Mfuzz")
    
    output$cluster_summary <- renderPrint({
      
      req(mfuzz_res())
      
      cl <- mfuzz_res()
      
      cat("Number of clusters:", length(unique(cl$cluster)), "\n")
      cat("Features clustered:", length(cl$cluster), "\n\n")
      print(table(cl$cluster))
    })
    
    # =====================================================
    # 11) Membership and Cluster Tables
    # =====================================================
    
    membership_df <- reactive({
      
      req(mfuzz_res())
      
      cl <- mfuzz_res()
      
      mem <- as.data.frame(cl$membership)
      colnames(mem) <- paste0("Cluster_", seq_len(ncol(mem)))
      
      mem %>%
        rownames_to_column("Feature")
    })
    
    cluster_table_df <- reactive({
      
      req(mfuzz_res())
      
      cl <- mfuzz_res()
      
      tibble(
        Feature = names(cl$cluster),
        Assigned_cluster = cl$cluster,
        Max_membership = apply(cl$membership, 1, max, na.rm = TRUE)
      ) %>%
        filter(Max_membership >= input$membership_cutoff) %>%
        arrange(
          Assigned_cluster,
          desc(Max_membership),
          Feature
        )
    })
    
    output$cluster_table <- renderDT({
      
      req(cluster_table_df())
      
      datatable(
        cluster_table_df(),
        options = list(
          pageLength = 15,
          scrollX = TRUE
        ),
        rownames = FALSE
      )
    })
    
    # =====================================================
    # 12) Combined Data for Cluster Profiles + Heatmap
    # =====================================================
    
    cluster_plot_data <- reactive({
      
      req(processed_matrix(), mfuzz_res())
      
      mat <- processed_matrix()
      cl <- mfuzz_res()
      
      cluster_levels_num <- sort(unique(cl$cluster))
      cluster_sizes <- table(cl$cluster)
      
      cluster_labels <- paste0(
        "Cluster ",
        cluster_levels_num,
        "\n",
        "n = ",
        as.integer(cluster_sizes[as.character(cluster_levels_num)])
      )
      
      gene_info <- tibble(
        Feature = rownames(mat),
        Cluster_num = unname(cl$cluster[rownames(mat)]),
        membership_max = apply(
          cl$membership[rownames(mat), , drop = FALSE],
          1,
          max,
          na.rm = TRUE
        )
      ) %>%
        arrange(
          Cluster_num,
          desc(membership_max),
          Feature
        ) %>%
        group_by(Cluster_num) %>%
        mutate(Row_in_cluster = row_number()) %>%
        ungroup() %>%
        mutate(
          Cluster = factor(
            Cluster_num,
            levels = cluster_levels_num,
            labels = cluster_labels
          )
        )

      validate(
        need(
          !anyNA(gene_info$Cluster_num),
          "Clustering results are out of sync with the current preprocessing settings — click 'Run Mfuzz clustering' again."
        )
      )

      profile_df <- as.data.frame(mat) %>%
        rownames_to_column("Feature") %>%
        pivot_longer(
          cols = -Feature,
          names_to = "Condition",
          values_to = "Expression"
        ) %>%
        left_join(
          gene_info %>%
            select(
              Feature,
              Cluster,
              membership_max
            ),
          by = "Feature"
        ) %>%
        mutate(
          Condition = factor(
            Condition,
            levels = colnames(mat)
          )
        )
      
      cluster_means <- profile_df %>%
        group_by(
          Cluster,
          Condition
        ) %>%
        summarise(
          mean_expr = mean(Expression, na.rm = TRUE),
          .groups = "drop"
        )
      
      heatmap_df <- as.data.frame(mat) %>%
        rownames_to_column("Feature") %>%
        pivot_longer(
          cols = -Feature,
          names_to = "Condition",
          values_to = "Expression"
        ) %>%
        left_join(
          gene_info %>%
            select(
              Feature,
              Cluster,
              Row_in_cluster
            ),
          by = "Feature"
        ) %>%
        mutate(
          Condition = factor(
            Condition,
            levels = colnames(mat)
          )
        ) %>%
        group_by(Cluster) %>%
        mutate(
          Row_plot = dplyr::n() - Row_in_cluster + 1L
        ) %>%
        ungroup()
      
      list(
        mat = mat,
        cl = cl,
        gene_info = gene_info,
        profile_df = profile_df,
        cluster_means = cluster_means,
        heatmap_df = heatmap_df
      )
    })
    
    # =====================================================
    # 13) Cluster Profile Plot
    # =====================================================
    
    cluster_profile_plot_obj <- reactive({
      
      req(cluster_plot_data())
      
      plot_data <- cluster_plot_data()
      
      ggplot() +
        geom_line(
          data = plot_data$profile_df,
          aes(
            x = Condition,
            y = Expression,
            group = Feature,
            color = membership_max
          ),
          alpha = 0.35,
          linewidth = 0.4
        ) +
        geom_line(
          data = plot_data$cluster_means,
          aes(
            x = Condition,
            y = mean_expr,
            group = Cluster
          ),
          color = input$profile_mean_col,
          linewidth = 1.2
        ) +
        scale_color_gradientn(
          colours = mfuzz_membership_colours(
            input$profile_col_low,
            input$profile_col_mid,
            input$profile_col_high
          ),
          limits = c(0, 1),
          breaks = c(0, 0.5, 1),
          labels = c("0", "0.5", "1"),
          name = "Membership"
        ) +
        facet_grid(
          Cluster ~ .,
          switch = "y"
        ) +
        theme_bw(base_size = input$font_size) +
        theme(
          strip.placement = "outside",
          strip.background = element_rect(
            fill = "grey95",
            colour = "grey60"
          ),
          strip.text.y.left = element_text(
            angle = 0,
            size = input$font_size,
            color = "black"
          ),
          panel.grid = element_blank(),
          panel.spacing.y = grid::unit(0.5, "lines"),
          axis.text.x = element_text(angle = 0, size = input$font_size, color = "black"),
          axis.text.y = element_text(size = input$font_size, color = "black"),
          axis.title = element_text(size = input$font_size, color = "black"),
          plot.title = element_text(size = input$font_size + 2, face = "bold", color = "black", hjust = 0.5),
          legend.position = "top",
          legend.title = element_text(size = input$font_size, color = "black"),
          legend.text = element_text(size = input$font_size, color = "black"),
          plot.margin = margin(5.5, 5.5, 5.5, 5.5)
        ) +
        labs(
          title = input$plot_title,
          x = NULL,
          y = if (isTRUE(input$standardize)) "Standardized expression" else "Expression"
        )
    })
    
    # =====================================================
    # 14) Heatmap Plot
    # =====================================================
    
    heatmap_plot_obj <- reactive({
      
      req(cluster_plot_data())
      
      plot_data <- cluster_plot_data()
      heatmap_df <- plot_data$heatmap_df
      
      fill_scale <- scale_fill_gradient2(
        low = input$heatmap_col_low,
        mid = input$heatmap_col_mid,
        high = input$heatmap_col_high,
        midpoint = 0,
        name = "Expression"
      )
      
      ggplot(
        heatmap_df,
        aes(
          x = Condition,
          y = Row_plot,
          fill = Expression
        )
      ) +
        geom_tile(
          width = 0.98,
          height = 0.98
        ) +
        fill_scale +
        facet_grid(
          Cluster ~ .,
          scales = "free_y",
          space = "fixed"
        ) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_bw(base_size = input$font_size) +
        theme(
          strip.text.y = element_blank(),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          panel.spacing.y = grid::unit(0.5, "lines"),
          axis.title = element_blank(),
          axis.text.x = element_text(size = input$font_size, color = "black"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "top",
          legend.title = element_text(size = input$font_size, color = "black"),
          legend.text = element_text(size = input$font_size, color = "black"),
          plot.margin = margin(5.5, 5.5, 5.5, 0)
        )
    })
    
    # =====================================================
    # 15) Combined Panel UI
    # =====================================================
    
    output$combined_cluster_panel_ui <- renderUI({
      
      req(mfuzz_res())
      
      n_clusters <- length(sort(unique(mfuzz_res()$cluster)))
      panel_height <- max(700, 120 * n_clusters + 120)
      
      plotOutput(
        ns("combined_cluster_panel"),
        height = paste0(panel_height, "px")
      )
    })
    
    # =====================================================
    # 16) Combined Aligned Panel
    # =====================================================
    
    combined_cluster_panel_obj <- reactive({
      
      req(
        cluster_profile_plot_obj(),
        heatmap_plot_obj()
      )
      
      profile_plot_main <- cluster_profile_plot_obj() +
        theme(legend.position = "none")
      
      heatmap_plot_main <- heatmap_plot_obj() +
        theme(legend.position = "none")
      
      profile_legend <- cowplot::get_legend(
        cluster_profile_plot_obj() +
          theme(
            legend.position = "top",
            legend.justification = "center"
          )
      )
      
      heatmap_legend <- cowplot::get_legend(
        heatmap_plot_obj() +
          theme(
            legend.position = "top",
            legend.justification = "center"
          )
      )
      
      legends_row <- cowplot::plot_grid(
        profile_legend,
        heatmap_legend,
        nrow = 1,
        rel_widths = c(1.15, 1)
      )
      
      panels_row <- cowplot::plot_grid(
        profile_plot_main,
        heatmap_plot_main,
        nrow = 1,
        rel_widths = c(1.15, 1),
        align = "h",
        axis = "tb",
        greedy = FALSE
      )
      
      cowplot::plot_grid(
        legends_row,
        panels_row,
        ncol = 1,
        rel_heights = c(0.12, 1)
      )
    })
    
    output$combined_cluster_panel <- renderPlot({
      
      req(combined_cluster_panel_obj())
      
      combined_cluster_panel_obj()
      
    }, res = 120)
    
    # =====================================================
    # 17) Plot Export Helper
    # =====================================================
    
    open_plot_device <- function(file, format, width, height, dpi) {
      
      if (format == "pdf") {
        
        pdf(
          file = file,
          width = width,
          height = height
        )
        
      } else if (format == "svg") {
        
        svglite::svglite(
          file = file,
          width = width,
          height = height
        )
        
      } else if (format == "png") {
        
        png(
          filename = file,
          width = width,
          height = height,
          units = "in",
          res = dpi
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
      }
    }
    
    # =====================================================
    # 18) Downloads
    # =====================================================
    
    output$download_combined_plot <- downloadHandler(
      
      filename = function() {
        paste0(
          "mfuzz_cluster_profiles_heatmap_",
          Sys.Date(),
          "_",
          input$plot_width, "x", input$plot_height,
          ".",
          input$plot_format
        )
      },
      
      content = function(file) {
        
        width_in <- input$plot_width * 0.393701
        height_in <- input$plot_height * 0.393701
        
        open_plot_device(
          file = file,
          format = input$plot_format,
          width = width_in,
          height = height_in,
          dpi = input$plot_dpi
        )
        
        print(combined_cluster_panel_obj())
        dev.off()
      }
    )
    
    output$download_clusters <- downloadHandler(
      
      filename = function() {
        paste0(
          "mfuzz_cluster_assignments_",
          Sys.Date(),
          ".txt"
        )
      },
      
      content = function(file) {
        
        write.table(
          cluster_table_df(),
          file = file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )
      }
    )
    
    output$download_membership <- downloadHandler(
      
      filename = function() {
        paste0(
          "mfuzz_membership_matrix_",
          Sys.Date(),
          ".txt"
        )
      },
      
      content = function(file) {
        
        write.table(
          membership_df(),
          file = file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )
      }
    )
    
    output$download_processed <- downloadHandler(
      
      filename = function() {
        paste0(
          "mfuzz_processed_matrix_",
          Sys.Date(),
          ".txt"
        )
      },
      
      content = function(file) {
        
        out <- as.data.frame(processed_matrix()) %>%
          rownames_to_column("Feature")
        
        write.table(
          out,
          file = file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )
      }
    )
    
    # =====================================================
    # 19) Return Reactive Objects
    # =====================================================
    
    return(
      list(
        processed_matrix = processed_matrix,
        mfuzz_res = mfuzz_res,
        membership_df = membership_df,
        cluster_table_df = cluster_table_df
      )
    )
  })
}
