# =========================================================
# EzOmicsViz — PCA / Sample Explorer Module
# =========================================================

library(shiny)
library(shinyjs)
library(tidyverse)
library(DT)
library(svglite)
library(colourpicker)
library(cowplot)
library(readxl)

options(shiny.maxRequestSize = 500 * 1024^2)

# =========================================================
# Helper Functions
# =========================================================

pca_smart_read_unchecked <- function(path, original_name = NULL) {
  
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
          as.data.frame()
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
    ),
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
      ),
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
      ),
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
      ),
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

pca_smart_read <- function(...) {
  omics_validate_table(pca_smart_read_unchecked(...))
}

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

collapse_duplicate_features <- function(expr_df, feature_ids, method = "highest_mean") {
  
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
  
  mat
}

row_var_base <- function(mat) {
  apply(mat, 1, var, na.rm = TRUE)
}

default_pca_group_colors <- function(n) {
  
  base_cols <- c(
    "#2A78D6",
    "#EB6834",
    "#1BAF7A",
    "#EDA100",
    "#E87BA4",
    "#008300",
    "#4A3AA7",
    "#E34948"
  )

  rep(base_cols, length.out = n)
}

# =========================================================
# UI
# =========================================================

pcaModuleUI <- function(id) {
  
  ns <- NS(id)
  
  tagList(
    
    useShinyjs(),
    
    tags$head(
      tags$style(HTML("
        .pca-panel {
          background: #ffffff;
          border-radius: 8px;
          padding: 20px;
          border: 1px solid #eeeeee;
          margin-bottom: 20px;
        }

        .pca-help {
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

        .pca-busy-indicator {
          display: none;
          position: fixed;
          right: 24px;
          bottom: 24px;
          z-index: 99999;
          background: #346b84;
          color: #ffffff;
          padding: 12px 18px;
          border-radius: 10px;
          box-shadow: 0 8px 24px rgba(16, 42, 67, 0.25);
          font-weight: 600;
          font-size: 0.95rem;
        }

        html.shiny-busy .pca-busy-indicator {
          display: block;
        }
      "))
    ),
    
    div(
      class = "pca-busy-indicator",
      "Processing PCA data..."
    ),
    
    fluidRow(
      
      column(
        width = 4,
        
        wellPanel(
          
          h4("1. Data Input"),
          
          fileInput(
            ns("expr_file"),
            "Upload normalized/log2 expression matrix (.txt, .csv, .tsv, .xls, .xlsx)",
            accept = c(".txt", ".csv", ".tsv", ".xls", ".xlsx")
          ),

          uiOutput(ns("id_column_ui")),

          p(
            class = "pca-help",
            "Input values must already be normalized and log2-transformed when appropriate. The module does not normalize or log-transform values. Duplicated feature identifiers are handled automatically by keeping the row with the highest mean abundance across samples."
          ),

          hr(),
          
          h4("2. Pre-processing"),
          
          selectInput(
            ns("feature_filter"),
            "Feature filtering before PCA:",
            choices = c(
              "Use all features" = "all",
              "Use top variable features" = "top_variable"
            ),
            selected = "top_variable"
          ),
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'top_variable'", ns("feature_filter")),
            
            numericInput(
              ns("top_variable_features"),
              "Number of top variable features:",
              value = 5000,
              min = 50,
              step = 50
            )
          ),
          
          checkboxInput(
            ns("center_pca"),
            "Center features before PCA",
            value = TRUE
          ),
          
          checkboxInput(
            ns("scale_pca"),
            "Scale features before PCA",
            value = TRUE
          ),
          
          p(
            class = "pca-help",
            "Features with missing, non-finite, or zero-variance values are removed before PCA."
          ),
          
          hr(),
          
          h4("3. Sample Annotation"),
          
          fileInput(
            ns("anno_file"),
            "Upload sample annotation (.txt, .csv, .tsv, .xls, .xlsx)",
            accept = c(".txt", ".csv", ".tsv", ".xls", ".xlsx")
          ),
          
          uiOutput(ns("annotation_ui")),
          
          p(
            class = "pca-help",
            "Optional. Annotation can be used to color, shape, label, and order samples in the PCA and correlation heatmap."
          ),
          
          hr(),
          
          h4("4. PCA Plot Options"),

          textInput(
            ns("plot_title"),
            "Plot Title",
            value = "PCA sample overview"
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
                ns("pc_x"),
                "X-axis PC:",
                value = 1,
                min = 1,
                step = 1
              )
            ),
            
            column(
              6,
              numericInput(
                ns("pc_y"),
                "Y-axis PC:",
                value = 2,
                min = 1,
                step = 1
              )
            )
          ),
          
          sliderInput(
            ns("point_size"),
            "Point size:",
            min = 1,
            max = 8,
            value = 3,
            step = 0.5
          ),
          
          checkboxInput(
            ns("show_sample_labels"),
            "Show sample labels in PCA plot",
            value = FALSE
          ),
          
          sliderInput(
            ns("label_size"),
            "Sample label size:",
            min = 2,
            max = 8,
            value = 3,
            step = 0.5
          ),
          
          checkboxInput(
            ns("show_ellipse"),
            "Show group ellipses when possible",
            value = TRUE
          ),
          
          sliderInput(
            ns("ellipse_level"),
            "Ellipse confidence level:",
            min = 0.50,
            max = 0.99,
            value = 0.95,
            step = 0.01
          ),
          
          checkboxInput(
            ns("show_corr_sample_names"),
            "Show sample names in correlation heatmap",
            value = FALSE
          ),
          
          selectInput(
            ns("corr_method"),
            "Sample correlation method:",
            choices = c(
              "Pearson" = "pearson",
              "Spearman" = "spearman"
            ),
            selected = "pearson"
          ),
          
          hr(),
          
          h4("5. Colors"),
          
          colourInput(
            ns("single_point_col"),
            "Point color when no annotation is selected",
            "#2A78D6"
          ),
          
          uiOutput(ns("pca_group_colors_ui")),
          
          colourInput(
            ns("corr_col_low"),
            "Correlation low",
            "#2A78D6"
          ),

          colourInput(
            ns("corr_col_mid"),
            "Correlation mid",
            "#F0EFEC"
          ),

          colourInput(
            ns("corr_col_high"),
            "Correlation high",
            "#E34948"
          ),
          
          hr(),
          
          h4("6. Export"),

          h5("PCA plot"),

          fluidRow(
            column(
              6,
              numericInput(
                ns("exp_width"),
                "Width (cm)",
                value = 14,
                min = 5,
                max = 60
              )
            ),

            column(
              6,
              numericInput(
                ns("exp_height"),
                "Height (cm)",
                value = 12,
                min = 5,
                max = 60
              )
            )
          ),

          h5("Scree plot"),

          fluidRow(
            column(
              6,
              numericInput(
                ns("exp_width_scree"),
                "Width (cm)",
                value = 14,
                min = 5,
                max = 60
              )
            ),

            column(
              6,
              numericInput(
                ns("exp_height_scree"),
                "Height (cm)",
                value = 10,
                min = 5,
                max = 60
              )
            )
          ),

          h5("Correlation heatmap"),

          fluidRow(
            column(
              6,
              numericInput(
                ns("exp_width_corr"),
                "Width (cm)",
                value = 15,
                min = 5,
                max = 60
              )
            ),

            column(
              6,
              numericInput(
                ns("exp_height_corr"),
                "Height (cm)",
                value = 13,
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
            ns("download_pca_plot"),
            "Download PCA plot",
            class = "btn-success",
            style = "width:100%"
          ),
          
          br(),
          br(),
          
          downloadButton(
            ns("download_scree_plot"),
            "Download scree plot",
            class = "btn-success",
            style = "width:100%"
          ),
          
          br(),
          br(),
          
          downloadButton(
            ns("download_corr_plot"),
            "Download correlation heatmap",
            class = "btn-success",
            style = "width:100%"
          ),
          
          br(),
          br(),
          
          downloadButton(
            ns("download_scores"),
            "Download PCA scores",
            class = "btn-info",
            style = "width:100%"
          ),
          
          br(),
          br(),
          
          downloadButton(
            ns("download_loadings"),
            "Download PCA loadings",
            class = "btn-info",
            style = "width:100%"
          ),
          
          br(),
          br(),
          
          downloadButton(
            ns("download_processed_matrix"),
            "Download processed matrix",
            class = "btn-secondary",
            style = "width:100%"
          )
        )
      ),
      
      column(
        width = 8,
        
        div(
          class = "pca-panel",
          h4("Input Summary"),
          verbatimTextOutput(ns("input_summary"))
        ),
        
        div(
          class = "pca-panel",
          h4("PCA Plot"),
          plotOutput(ns("pca_plot"), height = "620px")
        ),
        
        fluidRow(
          column(
            width = 6,
            div(
              class = "pca-panel",
              h4("Scree Plot"),
              plotOutput(ns("scree_plot"), height = "360px")
            )
          ),
          
          column(
            width = 6,
            div(
              class = "pca-panel",
              h4("Sample Correlation Heatmap"),
              plotOutput(ns("sample_corr_heatmap"), height = "360px")
            )
          )
        ),
        
        div(
          class = "pca-panel",
          h4("PCA Tables"),
          tabsetPanel(
            tabPanel(
              "Scores",
              br(),
              DTOutput(ns("scores_table"))
            ),
            tabPanel(
              "Loadings",
              br(),
              DTOutput(ns("loadings_table"))
            ),
            tabPanel(
              "Processed matrix",
              br(),
              DTOutput(ns("processed_matrix_table"))
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

pcaModuleServer <- function(id) {
  
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    # =====================================================
    # 1) Raw Data
    # =====================================================
    
    raw_data <- reactive({
      
      req(input$expr_file)
      
      pca_smart_read(
        path = input$expr_file$datapath,
        original_name = input$expr_file$name
      )
    })
    
    # =====================================================
    # 2) ID Column UI
    # =====================================================
    
    output$id_column_ui <- renderUI({
      
      req(raw_data())
      
      cols <- colnames(raw_data())
      
      selectInput(
        ns("id_col"),
        "Gene/protein/feature ID column:",
        choices = cols,
        selected = cols[1]
      )
    })
    
    # =====================================================
    # 3) Annotation UI
    # =====================================================
    
    output$annotation_ui <- renderUI({
      
      req(input$anno_file)
      
      anno <- pca_smart_read(
        path = input$anno_file$datapath,
        original_name = input$anno_file$name
      )
      
      cols <- colnames(anno)
      
      tagList(
        selectInput(
          ns("anno_sample_col"),
          "Sample column:",
          choices = cols,
          selected = cols[1]
        ),
        
        selectInput(
          ns("color_col"),
          "Color samples by:",
          choices = c("None" = "__none__", cols),
          selected = if (length(cols) >= 2) cols[2] else "__none__"
        ),
        
        selectInput(
          ns("shape_col"),
          "Shape samples by:",
          choices = c("None" = "__none__", cols),
          selected = "__none__"
        )
      )
    })
    
    # =====================================================
    # 4) Raw Expression Matrix
    # =====================================================
    
    raw_expr_matrix <- reactive({
      
      req(
        raw_data(),
        input$id_col
      )
      
      df <- raw_data()
      
      validate(
        need(
          input$id_col %in% colnames(df),
          "Selected ID column not found."
        )
      )
      
      feature_ids <- df[[input$id_col]]
      
      expr_df <- df[, setdiff(colnames(df), input$id_col), drop = FALSE]
      
      expr_df <- omics_numeric_frame(expr_df)
      
      validate(
        need(
          ncol(expr_df) >= 2,
          "No numeric sample columns were found after removing the feature ID column."
        )
      )
      
      mat <- collapse_duplicate_features(
        expr_df = expr_df,
        feature_ids = feature_ids,
        method = "highest_mean"
      )
      
      storage.mode(mat) <- "numeric"
      
      validate(
        need(
          nrow(mat) >= 2,
          "Too few features were found."
        ),
        need(
          ncol(mat) >= 3,
          "PCA requires at least 3 samples."
        )
      )
      
      mat
    })
    
    # =====================================================
    # 5) Duplicate Summary
    # =====================================================
    
    duplicate_summary <- reactive({
      
      req(raw_data(), input$id_col)
      
      df <- raw_data()
      
      ids <- df[[input$id_col]]
      
      ids <- as.character(ids) %>%
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
        raw_expr_matrix(),
        input$feature_filter
      )
      
      mat <- raw_expr_matrix()
      
      keep_some_finite <- rowSums(is.finite(mat)) > 0
      mat <- mat[keep_some_finite, , drop = FALSE]
      
      keep_complete <- apply(mat, 1, function(x) {
        all(is.finite(x))
      })
      
      mat <- mat[keep_complete, , drop = FALSE]
      
      validate(
        need(
          nrow(mat) >= 2,
          "Too few features remain after removing missing or non-finite values."
        )
      )
      
      rv <- row_var_base(mat)
      rv[!is.finite(rv)] <- NA_real_
      
      keep_var <- !is.na(rv) & rv > 0
      mat <- mat[keep_var, , drop = FALSE]
      rv <- rv[keep_var]
      
      validate(
        need(
          nrow(mat) >= 2,
          "Too few non-zero variance features remain for PCA."
        )
      )
      
      if (input$feature_filter == "top_variable") {
        omics_number_parameter(input$top_variable_features, "Number of variable features", 2, .Machine$integer.max, TRUE)
        
        top_n <- min(
          input$top_variable_features,
          nrow(mat)
        )
        
        keep_idx <- order(
          rv,
          decreasing = TRUE
        )[seq_len(top_n)]
        
        mat <- mat[keep_idx, , drop = FALSE]
      }
      
      validate(
        need(
          nrow(mat) >= 2,
          "Too few features remain after variance filtering."
        ),
        need(
          ncol(mat) >= 3,
          "PCA requires at least 3 samples."
        )
      )
      
      mat
    })
    
    # =====================================================
    # 7) Sample Annotation
    # =====================================================
    
    sample_annotation <- reactive({
      
      if (is.null(input$anno_file)) {
        return(NULL)
      }
      
      req(
        processed_matrix(),
        input$anno_sample_col
      )
      
      anno <- pca_smart_read(
        path = input$anno_file$datapath,
        original_name = input$anno_file$name
      )
      
      validate(
        need(
          input$anno_sample_col %in% colnames(anno),
          "Selected sample column not found in annotation."
        )
      )
      
      anno <- omics_annotation(
        anno, input$anno_sample_col, colnames(processed_matrix()),
        required = c(input$color_col, input$shape_col)
      )

      anno <- anno %>%
        mutate(
          Sample = as.character(.data[[input$anno_sample_col]])
        ) %>%
        filter(
          !is.na(Sample),
          Sample != ""
        ) %>%
        distinct(
          Sample,
          .keep_all = TRUE
        )
      
      samples <- colnames(processed_matrix())
      
      out <- tibble(Sample = samples) %>%
        left_join(
          anno,
          by = "Sample"
        )
      
      out$Sample <- factor(
        out$Sample,
        levels = samples
      )
      
      out
    })
    
    sample_order <- reactive({
      
      anno <- sample_annotation()
      
      if (is.null(anno)) {
        return(colnames(processed_matrix()))
      }
      
      color_col <- input$color_col %||% "__none__"
      
      if (color_col != "__none__" && color_col %in% colnames(anno)) {
        
        anno %>%
          arrange(.data[[color_col]], Sample) %>%
          pull(Sample) %>%
          as.character()
        
      } else {
        
        as.character(anno$Sample)
      }
    })
    
    # =====================================================
    # 8) PCA Group Colors
    # =====================================================
    
    pca_color_column <- reactive({
      
      anno <- sample_annotation()
      
      if (is.null(anno)) {
        return(NULL)
      }
      
      color_col <- input$color_col %||% "__none__"
      shape_col <- input$shape_col %||% "__none__"
      
      if (color_col != "__none__" && color_col %in% colnames(anno)) {
        return(color_col)
      }
      
      if (shape_col != "__none__" && shape_col %in% colnames(anno)) {
        return(shape_col)
      }
      
      NULL
    })
    
    pca_group_levels <- reactive({
      
      anno <- sample_annotation()
      color_col <- pca_color_column()
      
      if (is.null(anno) || is.null(color_col)) {
        return(NULL)
      }
      
      vals <- anno[[color_col]]
      
      color_is_discrete <- is.character(vals) ||
        is.factor(vals) ||
        dplyr::n_distinct(vals, na.rm = TRUE) <= 12
      
      if (!isTRUE(color_is_discrete)) {
        return(NULL)
      }
      
      vals <- as.character(vals)
      vals <- vals[!is.na(vals) & vals != ""]
      
      if (length(vals) == 0) {
        return(NULL)
      }
      
      sort(unique(vals))
    })
    
    output$pca_group_colors_ui <- renderUI({
      
      levs <- pca_group_levels()
      color_col <- pca_color_column()
      
      if (is.null(levs) || length(levs) == 0 || is.null(color_col)) {
        return(NULL)
      }
      
      default_cols <- default_pca_group_colors(length(levs))
      
      tagList(
        h5("PCA group colors"),
        
        p(
          class = "pca-help",
          paste0("Colors used for: ", color_col)
        ),
        
        lapply(seq_along(levs), function(i) {
          
          colourInput(
            inputId = ns(paste0("pca_group_col_", i)),
            label = levs[i],
            value = default_cols[i]
          )
        })
      )
    })
    
    pca_group_colors <- reactive({
      
      levs <- pca_group_levels()
      
      if (is.null(levs) || length(levs) == 0) {
        return(NULL)
      }
      
      default_cols <- default_pca_group_colors(length(levs))
      
      cols <- vapply(
        seq_along(levs),
        function(i) {
          input[[paste0("pca_group_col_", i)]] %||% default_cols[i]
        },
        character(1)
      )
      
      names(cols) <- levs
      
      cols
    })
    
    # =====================================================
    # 9) PCA Object
    # =====================================================
    
    pca_obj <- reactive({
      
      req(processed_matrix())
      
      mat <- processed_matrix()
      pca_input <- t(mat)
      
      validate(
        need(
          nrow(pca_input) >= 3,
          "PCA requires at least 3 samples."
        ),
        need(
          ncol(pca_input) >= 2,
          "PCA requires at least 2 features."
        )
      )
      
      center_pca <- isTRUE(input$center_pca)
      scale_pca <- isTRUE(input$scale_pca)
      
      stats::prcomp(
        pca_input,
        center = center_pca,
        scale. = scale_pca
      )
    })
    
    pca_variance <- reactive({
      
      req(pca_obj())
      
      pca <- pca_obj()
      
      var_explained <- pca$sdev^2 / sum(pca$sdev^2)
      
      tibble(
        PC = paste0("PC", seq_along(var_explained)),
        PC_number = seq_along(var_explained),
        Variance = var_explained,
        Variance_percent = 100 * var_explained,
        Cumulative_percent = 100 * cumsum(var_explained)
      )
    })
    
    pca_scores <- reactive({
      
      req(pca_obj())
      
      scores <- as.data.frame(pca_obj()$x) %>%
        rownames_to_column("Sample")
      
      anno <- sample_annotation()
      
      if (!is.null(anno)) {
        
        scores <- scores %>%
          left_join(
            anno %>%
              mutate(Sample = as.character(Sample)),
            by = "Sample"
          )
      }
      
      scores
    })
    
    pca_loadings <- reactive({
      
      req(pca_obj())
      
      as.data.frame(pca_obj()$rotation) %>%
        rownames_to_column("Feature")
    })
    
    # =====================================================
    # 10) Input Summary
    # =====================================================
    
    output$input_summary <- renderPrint({
      
      req(processed_matrix(), duplicate_summary())
      
      raw_mat <- raw_expr_matrix()
      mat <- processed_matrix()
      dup <- duplicate_summary()
      var_df <- pca_variance()
      
      cat("Raw features after duplicate handling:", nrow(raw_mat), "\n")
      cat("Processed features used for PCA:", nrow(mat), "\n")
      cat("Samples:", ncol(mat), "\n\n")
      
      cat("Raw feature IDs:", dup$total_ids, "\n")
      cat("Unique raw feature IDs:", dup$unique_ids, "\n")
      cat("Duplicated rows:", dup$duplicated_rows, "\n")
      cat("Feature names with duplicates:", dup$duplicated_feature_names, "\n\n")
      
      cat("Feature filter:", input$feature_filter, "\n")
      
      if (input$feature_filter == "top_variable") {
        omics_number_parameter(input$top_variable_features, "Number of variable features", 2, .Machine$integer.max, TRUE)
        cat("Top variable features requested:", input$top_variable_features, "\n")
      }
      
      cat("PCA center:", isTRUE(input$center_pca), "\n")
      cat("PCA scale:", isTRUE(input$scale_pca), "\n\n")
      
      cat("Processed matrix range:\n")
      print(range(mat, na.rm = TRUE))
      
      cat("\nVariance explained:\n")
      print(
        var_df %>%
          slice_head(n = min(6, nrow(var_df))) %>%
          select(PC, Variance_percent, Cumulative_percent)
      )
    })
    
    # =====================================================
    # 11) PCA Plot
    # =====================================================
    
    create_pca_plot <- reactive({
      
      req(
        pca_scores(),
        pca_variance(),
        input$pc_x,
        input$pc_y
      )
      
      df <- pca_scores()
      var_df <- pca_variance()
      
      max_pc <- nrow(var_df)
      omics_number_parameter(input$pc_x, "X principal component", 1, max_pc, TRUE)
      omics_number_parameter(input$pc_y, "Y principal component", 1, max_pc, TRUE)
      
      validate(
        need(
          input$pc_x <= max_pc,
          paste0("X-axis PC exceeds available PCs. Maximum available PC: ", max_pc)
        ),
        need(
          input$pc_y <= max_pc,
          paste0("Y-axis PC exceeds available PCs. Maximum available PC: ", max_pc)
        ),
        need(
          input$pc_x != input$pc_y,
          "X-axis PC and Y-axis PC must be different."
        )
      )
      
      pcx <- paste0("PC", input$pc_x)
      pcy <- paste0("PC", input$pc_y)
      
      x_lab <- paste0(
        pcx,
        " (",
        round(var_df$Variance_percent[input$pc_x], 1),
        "%)"
      )
      
      y_lab <- paste0(
        pcy,
        " (",
        round(var_df$Variance_percent[input$pc_y], 1),
        "%)"
      )
      
      color_col <- input$color_col %||% "__none__"
      shape_col <- input$shape_col %||% "__none__"
      
      use_color <- color_col != "__none__" && color_col %in% colnames(df)
      use_shape <- shape_col != "__none__" && shape_col %in% colnames(df)
      
      if (!isTRUE(use_color) && isTRUE(use_shape)) {
        color_col <- shape_col
        use_color <- TRUE
      }
      
      color_is_discrete <- FALSE
      
      if (isTRUE(use_color)) {
        
        color_values <- df[[color_col]]
        
        color_is_discrete <- is.character(color_values) ||
          is.factor(color_values) ||
          dplyr::n_distinct(color_values, na.rm = TRUE) <= 12
        
        if (isTRUE(color_is_discrete)) {
          
          levs <- pca_group_levels()
          
          if (!is.null(levs)) {
            df[[color_col]] <- factor(df[[color_col]], levels = levs)
          } else {
            df[[color_col]] <- factor(df[[color_col]])
          }
        }
      }
      
      if (isTRUE(use_shape)) {
        df[[shape_col]] <- factor(df[[shape_col]])
      }
      
      if (isTRUE(use_color) && isTRUE(use_shape)) {
        
        p <- ggplot(
          df,
          aes(
            x = .data[[pcx]],
            y = .data[[pcy]],
            color = .data[[color_col]],
            shape = .data[[shape_col]]
          )
        )
        
      } else if (isTRUE(use_color)) {
        
        p <- ggplot(
          df,
          aes(
            x = .data[[pcx]],
            y = .data[[pcy]],
            color = .data[[color_col]]
          )
        )
        
      } else {
        
        p <- ggplot(
          df,
          aes(
            x = .data[[pcx]],
            y = .data[[pcy]]
          )
        )
      }
      
      if (isTRUE(input$show_ellipse)) {
        
        if (isTRUE(use_color) && isTRUE(color_is_discrete)) {
          
          ellipse_df <- df %>%
            filter(!is.na(.data[[color_col]])) %>%
            group_by(.data[[color_col]]) %>%
            filter(n() >= 3) %>%
            ungroup()
          
          if (nrow(ellipse_df) > 0) {
            
            p <- p +
              stat_ellipse(
                data = ellipse_df,
                aes(
                  group = .data[[color_col]],
                  color = .data[[color_col]]
                ),
                level = input$ellipse_level,
                linewidth = 0.5,
                linetype = "dashed",
                show.legend = FALSE
              )
          }
          
        } else if (nrow(df) >= 3) {
          
          p <- p +
            stat_ellipse(
              level = input$ellipse_level,
              linewidth = 0.5,
              linetype = "dashed",
              color = "grey40"
            )
        }
      }
      
      if (isTRUE(use_color)) {
        
        p <- p +
          geom_point(
            size = input$point_size,
            alpha = 0.88
          )
        
      } else {
        
        p <- p +
          geom_point(
            size = input$point_size,
            alpha = 0.88,
            color = input$single_point_col
          )
      }
      
      if (isTRUE(use_color) && isTRUE(color_is_discrete)) {
        
        manual_cols <- pca_group_colors()
        
        if (!is.null(manual_cols)) {
          
          p <- p +
            scale_color_manual(
              values = manual_cols,
              na.value = "grey70"
            )
        }
      }
      
      if (isTRUE(input$show_sample_labels)) {
        
        p <- p +
          geom_text(
            aes(label = Sample),
            vjust = -0.8,
            size = input$label_size,
            check_overlap = TRUE,
            show.legend = FALSE
          )
      }
      
      p +
        geom_hline(
          yintercept = 0,
          linewidth = 0.3,
          color = "grey85"
        ) +
        geom_vline(
          xintercept = 0,
          linewidth = 0.3,
          color = "grey85"
        ) +
        theme_bw(base_size = input$font_size) +
        labs(
          title = input$plot_title,
          x = x_lab,
          y = y_lab,
          color = if (isTRUE(use_color)) color_col else NULL,
          shape = if (isTRUE(use_shape)) shape_col else NULL
        ) +
        theme(
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = input$font_size, color = "black"),
          axis.title = element_text(size = input$font_size, color = "black"),
          plot.title = element_text(size = input$font_size + 2, face = "bold", color = "black", hjust = 0.5),
          legend.title = element_text(color = "black", size = input$font_size),
          legend.text = element_text(color = "black", size = input$font_size)
        )
    })
    
    output$pca_plot <- renderPlot({
      
      req(create_pca_plot())
      
      print(create_pca_plot())
      
    }, res = 120)
    
    # =====================================================
    # 12) Scree Plot
    # =====================================================
    
    create_scree_plot <- reactive({
      
      req(pca_variance())
      
      df <- pca_variance()
      
      top_n <- min(20, nrow(df))
      df <- df[seq_len(top_n), , drop = FALSE]
      
      ggplot(
        df,
        aes(
          x = factor(PC, levels = PC),
          y = Variance_percent,
          group = 1
        )
      ) +
        geom_col(
          fill = "#2A78D6",
          alpha = 0.85
        ) +
        geom_line(
          aes(y = Cumulative_percent, color = "Cumulative Explained Variance"),
          linewidth = 0.7
        ) +
        geom_point(
          aes(y = Cumulative_percent, color = "Cumulative Explained Variance"),
          size = 3
        ) +
        scale_color_manual(
          values = c("Cumulative Explained Variance" = "grey30")
        ) +
        theme_bw(base_size = input$font_size) +
        labs(
          title = "Variance explained",
          x = NULL,
          y = "Variance explained (%)"
        ) +
        theme(
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = input$font_size, color = "black"),
          axis.text = element_text(size = input$font_size, color = "black"),
          axis.title = element_text(size = input$font_size, color = "black"),
          plot.title = element_text(face = "bold", color = "black", size = input$font_size + 2, hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = input$font_size, color = "black"),
          legend.margin = margin(t = -5)
        )
    })
    
    output$scree_plot <- renderPlot({
      
      req(create_scree_plot())
      
      print(create_scree_plot())
      
    }, res = 120)
    
    # =====================================================
    # 13) Sample Correlation Heatmap
    # =====================================================
    
    create_corr_heatmap <- reactive({
      
      req(processed_matrix())
      
      mat <- processed_matrix()
      
      corr <- suppressWarnings(
        cor(
          mat,
          method = input$corr_method,
          use = "pairwise.complete.obs"
        )
      )
      
      corr[!is.finite(corr)] <- NA_real_
      
      ord <- sample_order()
      ord <- ord[ord %in% colnames(corr)]
      
      corr <- corr[ord, ord, drop = FALSE]
      
      anno <- sample_annotation()
      color_col <- input$color_col %||% "__none__"
      
      has_group_annotation <- !is.null(anno) &&
        color_col != "__none__" &&
        color_col %in% colnames(anno)
      
      group_bar_df <- NULL
      group_cols <- NULL
      
      if (isTRUE(has_group_annotation)) {
        
        group_bar_df <- anno %>%
          mutate(
            Sample = as.character(Sample),
            Group = as.character(.data[[color_col]])
          ) %>%
          filter(
            Sample %in% ord,
            !is.na(Group),
            Group != ""
          ) %>%
          distinct(Sample, .keep_all = TRUE)
        
        group_bar_df <- tibble(Sample = ord) %>%
          left_join(
            group_bar_df %>%
              select(Sample, Group),
            by = "Sample"
          )
        
        group_bar_df$Sample_x <- factor(
          group_bar_df$Sample,
          levels = ord
        )
        
        group_bar_df$Sample_y <- factor(
          group_bar_df$Sample,
          levels = rev(ord)
        )
        
        group_levels <- pca_group_levels()
        group_cols <- pca_group_colors()
        
        if (!is.null(group_levels)) {
          group_bar_df$Group <- factor(
            group_bar_df$Group,
            levels = group_levels
          )
        }
      }
      
      plot_df <- as.data.frame(corr) %>%
        rownames_to_column("Sample_1") %>%
        pivot_longer(
          cols = -Sample_1,
          names_to = "Sample_2",
          values_to = "Correlation"
        )
      
      plot_df$Sample_1 <- factor(
        plot_df$Sample_1,
        levels = rev(ord)
      )
      
      plot_df$Sample_2 <- factor(
        plot_df$Sample_2,
        levels = ord
      )
      
      show_sample_names <- isTRUE(input$show_corr_sample_names)
      
      x_text_theme <- if (show_sample_names) {
        element_text(
          angle = 60,
          hjust = 1,
          vjust = 1,
          size = 7,
          color = "black"
        )
      } else {
        element_blank()
      }
      
      y_text_theme <- if (show_sample_names) {
        element_text(
          size = 7,
          color = "black"
        )
      } else {
        element_blank()
      }
      
      base_heatmap <- ggplot(
        plot_df,
        aes(
          x = Sample_2,
          y = Sample_1,
          fill = Correlation
        )
      ) +
        geom_raster() +
        scale_fill_gradient2(
          low = input$corr_col_low,
          mid = input$corr_col_mid,
          high = input$corr_col_high,
          midpoint = 0,
          limits = c(-1, 1),
          na.value = "grey90"
        ) +
        theme_bw(base_size = input$font_size) +
        theme(
          axis.title = element_blank(),
          axis.text.x = x_text_theme,
          axis.text.y = y_text_theme,
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_blank(),
          legend.title = element_text(size = input$font_size, color = "black"),
          legend.text = element_text(size = input$font_size, color = "black"),
          plot.margin = margin(2, 2, 2, 2)
        ) +
        labs(
          title = NULL,
          fill = "r"
        )
      
      corr_legend <- cowplot::get_legend(
        base_heatmap +
          theme(
            legend.position = "right",
            legend.box.margin = margin(0, 0, 0, 4)
          )
      )
      
      main_heatmap <- base_heatmap +
        theme(
          legend.position = "none"
        )
      
      if (!isTRUE(has_group_annotation) || is.null(group_bar_df) || is.null(group_cols)) {
        
        return(
          cowplot::plot_grid(
            main_heatmap,
            corr_legend,
            ncol = 2,
            rel_widths = c(1, 0.16),
            align = "h"
          )
        )
      }
      
      top_annotation <- ggplot(
        group_bar_df,
        aes(
          x = Sample_x,
          y = 1,
          fill = Group
        )
      ) +
        geom_tile() +
        scale_fill_manual(
          values = group_cols,
          na.value = "grey80"
        ) +
        theme_void() +
        theme(
          legend.position = "none",
          plot.margin = margin(0, 2, 0, 2)
        )
      
      left_annotation <- ggplot(
        group_bar_df,
        aes(
          x = 1,
          y = Sample_y,
          fill = Group
        )
      ) +
        geom_tile() +
        scale_fill_manual(
          values = group_cols,
          na.value = "grey80"
        ) +
        theme_void() +
        theme(
          legend.position = "none",
          plot.margin = margin(2, 0, 2, 0)
        )
      
      empty_corner <- ggplot() +
        theme_void()
      
      empty_legend_top <- ggplot() +
        theme_void()
      
      upper_row <- cowplot::plot_grid(
        empty_corner,
        top_annotation,
        empty_legend_top,
        ncol = 3,
        rel_widths = c(0.025, 1, 0.16),
        align = "h"
      )
      
      lower_row <- cowplot::plot_grid(
        left_annotation,
        main_heatmap,
        corr_legend,
        ncol = 3,
        rel_widths = c(0.025, 1, 0.16),
        align = "h"
      )
      
      cowplot::plot_grid(
        upper_row,
        lower_row,
        ncol = 1,
        rel_heights = c(0.025, 1),
        align = "v"
      )
    })
    
    output$sample_corr_heatmap <- renderPlot({
      
      req(create_corr_heatmap())
      
      print(create_corr_heatmap())
      
    }, res = 120)
    
    # =====================================================
    # 14) Tables
    # =====================================================
    
    output$scores_table <- renderDT({
      
      req(pca_scores())
      
      datatable(
        pca_scores(),
        rownames = FALSE,
        options = list(
          pageLength = 10,
          scrollX = TRUE
        )
      )
    })
    
    output$loadings_table <- renderDT({
      
      req(pca_loadings())
      
      datatable(
        pca_loadings(),
        rownames = FALSE,
        options = list(
          pageLength = 10,
          scrollX = TRUE
        )
      )
    })
    
    processed_matrix_df <- reactive({
      
      req(processed_matrix())
      
      as.data.frame(processed_matrix()) %>%
        rownames_to_column("Feature")
    })
    
    output$processed_matrix_table <- renderDT({
      
      req(processed_matrix_df())
      
      preview_df <- processed_matrix_df() %>%
        slice_head(n = 1000)
      
      datatable(
        preview_df,
        rownames = FALSE,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          deferRender = TRUE
        ),
        caption = htmltools::tags$caption(
          style = "caption-side: top; text-align: left;",
          paste0(
            "Preview showing first ",
            nrow(preview_df),
            " features. The downloaded processed matrix contains all features."
          )
        )
      )
    })
    
    # =====================================================
    # 15) Downloads
    # =====================================================
    
    output$download_pca_plot <- downloadHandler(
      
      filename = function() {
        paste0(
          "pca_plot_",
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
        
        open_plot_device(
          file = file,
          fmt = input$exp_fmt,
          width = width_in,
          height = height_in,
          dpi = input$exp_dpi
        )
        
        print(create_pca_plot())
        dev.off()
      }
    )
    
    output$download_scree_plot <- downloadHandler(

      filename = function() {
        paste0(
          "pca_scree_plot_",
          Sys.Date(),
          "_",
          input$exp_width_scree, "x", input$exp_height_scree,
          ".",
          input$exp_fmt
        )
      },

      content = function(file) {

        width_in <- input$exp_width_scree * 0.393701
        height_in <- input$exp_height_scree * 0.393701

        open_plot_device(
          file = file,
          fmt = input$exp_fmt,
          width = width_in,
          height = height_in,
          dpi = input$exp_dpi
        )

        print(create_scree_plot())
        dev.off()
      }
    )

    output$download_corr_plot <- downloadHandler(

      filename = function() {
        paste0(
          "pca_sample_correlation_heatmap_",
          Sys.Date(),
          "_",
          input$exp_width_corr, "x", input$exp_height_corr,
          ".",
          input$exp_fmt
        )
      },

      content = function(file) {

        width_in <- input$exp_width_corr * 0.393701
        height_in <- input$exp_height_corr * 0.393701

        open_plot_device(
          file = file,
          fmt = input$exp_fmt,
          width = width_in,
          height = height_in,
          dpi = input$exp_dpi
        )

        print(create_corr_heatmap())
        dev.off()
      }
    )
    
    output$download_scores <- downloadHandler(
      
      filename = function() {
        paste0(
          "pca_scores_",
          Sys.Date(),
          ".txt"
        )
      },
      
      content = function(file) {
        
        write.table(
          pca_scores(),
          file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )
      }
    )
    
    output$download_loadings <- downloadHandler(
      
      filename = function() {
        paste0(
          "pca_loadings_",
          Sys.Date(),
          ".txt"
        )
      },
      
      content = function(file) {
        
        write.table(
          pca_loadings(),
          file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )
      }
    )
    
    output$download_processed_matrix <- downloadHandler(
      
      filename = function() {
        paste0(
          "pca_processed_matrix_",
          Sys.Date(),
          ".txt"
        )
      },
      
      content = function(file) {
        
        write.table(
          processed_matrix_df(),
          file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )
      }
    )
    
    # =====================================================
    # 16) Return Reactive Objects
    # =====================================================
    
    return(
      list(
        raw_expr_matrix = raw_expr_matrix,
        processed_matrix = processed_matrix,
        sample_annotation = sample_annotation,
        pca_obj = pca_obj,
        pca_scores = pca_scores,
        pca_loadings = pca_loadings
      )
    )
  })
}