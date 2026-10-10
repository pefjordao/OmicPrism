# =========================================================
# EzOmicsViz — Expression Heatmap Module
# =========================================================

library(shiny)
library(shinyjs)
library(tidyverse)
library(DT)
library(svglite)
library(colourpicker)
library(readxl)

options(shiny.maxRequestSize = 500 * 1024^2)

# =========================================================
# Helper Functions
# =========================================================

hm_clean_column_names <- function(df) {
  
  colnames(df) <- colnames(df) %>%
    as.character() %>%
    stringr::str_trim()
  
  df
}

hm_smart_read_unchecked <- function(path, original_name = NULL) {
  
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
          hm_clean_column_names()
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
      hm_clean_column_names(),
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
        hm_clean_column_names(),
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
        hm_clean_column_names(),
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
        hm_clean_column_names(),
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

hm_smart_read <- function(...) {
  omics_validate_table(hm_smart_read_unchecked(...))
}

hm_read_feature_list <- function(path, original_name = NULL) {
  omics_read_id_list(path, original_name)
}

hm_open_plot_device <- function(file, fmt, width, height, dpi) {
  
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

hm_collapse_duplicate_features_highest_mean <- function(expr_df, feature_ids) {
  
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

hm_row_var <- function(mat) {
  apply(mat, 1, var, na.rm = TRUE)
}

hm_scale_rows <- function(mat) {
  
  scaled <- t(scale(t(mat)))
  scaled[!is.finite(scaled)] <- 0
  
  scaled
}

hm_get_cluster_order <- function(mat, margin = "rows", distance = "correlation", method = "complete") {
  
  margin <- match.arg(
    margin,
    choices = c("rows", "cols")
  )
  
  if (margin == "rows") {
    x <- mat
    labels <- rownames(mat)
  } else {
    x <- t(mat)
    labels <- colnames(mat)
  }
  
  if (nrow(x) < 2) {
    return(labels)
  }
  
  dist_obj <- tryCatch(
    {
      if (distance == "correlation") {
        
        corr <- suppressWarnings(
          cor(
            t(x),
            use = "pairwise.complete.obs",
            method = "pearson"
          )
        )
        
        corr[!is.finite(corr)] <- 0
        diag(corr) <- 1
        
        as.dist(1 - corr)
        
      } else {
        
        dist(
          x,
          method = "euclidean"
        )
      }
    },
    error = function(e) NULL
  )
  
  if (is.null(dist_obj)) {
    return(labels)
  }
  
  hc <- tryCatch(
    hclust(
      dist_obj,
      method = method
    ),
    error = function(e) NULL
  )
  
  if (is.null(hc)) {
    return(labels)
  }
  
  labels[hc$order]
}

# =========================================================
# UI
# =========================================================

heatmapModuleUI <- function(id) {
  
  ns <- NS(id)
  
  tagList(
    
    useShinyjs(),
    
    tags$head(
      tags$style(HTML("
        .heatmap-panel {
          background: #ffffff;
          border-radius: 8px;
          padding: 20px;
          border: 1px solid #eeeeee;
          margin-bottom: 20px;
        }

        .heatmap-help {
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
            "Upload normalized/log2 expression matrix (.txt, .csv, .tsv, .xls, .xlsx)",
            accept = c(".txt", ".csv", ".tsv", ".xls", ".xlsx")
          ),
          
          checkboxInput(
            ns("transpose"),
            "Matrix is transposed: samples in rows, features in columns",
            value = FALSE
          ),
          
          uiOutput(ns("id_column_ui")),
          
          p(
            class = "heatmap-help",
            "Input values must already be normalized and log2-transformed when appropriate. The module does not normalize, log-transform, or impute values. Duplicated feature identifiers are handled automatically by keeping the row with the highest mean abundance across samples."
          ),
          
          hr(),
          
          h4("3. Feature Selection"),
          
          selectInput(
            ns("feature_selection"),
            "Features to display:",
            choices = c(
              "Top variable features" = "top_variable",
              "Uploaded feature list" = "uploaded_list",
              "Manual feature selection" = "manual",
              "All features" = "all"
            ),
            selected = "top_variable"
          ),
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'top_variable'", ns("feature_selection")),
            
            numericInput(
              ns("top_variable_features"),
              "Number of top variable features:",
              value = 10,
              min = 2,
              step = 5
            )
          ),
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'uploaded_list'", ns("feature_selection")),
            
            fileInput(
              ns("feature_list_file"),
              "Upload feature list (.txt, .csv, .tsv, .xls, .xlsx)",
              accept = c(".txt", ".csv", ".tsv", ".xls", ".xlsx")
            )
          ),
          
          uiOutput(ns("manual_feature_ui")),
          
          numericInput(
            ns("max_features_display"),
            "Maximum features allowed in heatmap:",
            value = 500,
            min = 10,
            step = 10
          ),
          
          p(
            class = "heatmap-help",
            "Very large heatmaps are slow and usually unreadable. Use top variable features or an uploaded feature list for publication-style figures."
          ),
          
          hr(),
          
          h4("4. Sample Annotation"),
          
          fileInput(
            ns("anno_file"),
            "Upload sample annotation (.txt, .csv, .tsv, .xls, .xlsx)",
            accept = c(".txt", ".csv", ".tsv", ".xls", ".xlsx")
          ),
          
          uiOutput(ns("annotation_ui")),
          
          checkboxInput(
            ns("facet_by_group"),
            "Facet heatmap by group",
            value = TRUE
          ),
          
          checkboxInput(
            ns("order_samples_by_group"),
            "Order samples by group when not clustering columns",
            value = TRUE
          ),
          
          p(
            class = "heatmap-help",
            "Optional. Annotation is used to group samples visually and organize columns in the heatmap."
          ),
          
          hr(),
          
          h4("5. Clustering Options"),
          
          checkboxInput(
            ns("cluster_rows"),
            "Cluster features",
            value = TRUE
          ),
          
          checkboxInput(
            ns("cluster_cols"),
            "Cluster samples",
            value = FALSE
          ),
          
          selectInput(
            ns("distance_method"),
            "Clustering distance:",
            choices = c(
              "Correlation distance" = "correlation",
              "Euclidean distance" = "euclidean"
            ),
            selected = "correlation"
          ),
          
          selectInput(
            ns("hclust_method"),
            "Hierarchical clustering method:",
            choices = c(
              "Complete" = "complete",
              "Average" = "average",
              "Ward.D2" = "ward.D2",
              "Single" = "single"
            ),
            selected = "complete"
          ),
          
          hr(),
          
          h4("6. Plot Options"),

          textInput(
            ns("plot_title"),
            "Plot Title",
            value = "Expression heatmap"
          ),

          checkboxInput(
            ns("show_feature_names"),
            "Show feature names",
            value = TRUE
          ),
          
          checkboxInput(
            ns("show_sample_names"),
            "Show sample names",
            value = FALSE
          ),
          
          sliderInput(
            ns("font_size"),
            "Font size:",
            min = 6,
            max = 18,
            value = 12,
            step = 1
          ),
          
          sliderInput(
            ns("feature_label_size"),
            "Feature label size:",
            min = 4,
            max = 14,
            value = 8,
            step = 1
          ),
          
          sliderInput(
            ns("sample_label_size"),
            "Sample label size:",
            min = 4,
            max = 14,
            value = 7,
            step = 1
          ),
          
          checkboxInput(
            ns("tile_border"),
            "Show tile borders",
            value = FALSE
          ),
          
          hr(),
          
          h4("7. Colors"),
          
          colourInput(
            ns("col_low"),
            "Low",
            "#2A78D6"
          ),

          colourInput(
            ns("col_mid"),
            "Mid",
            "#F0EFEC"
          ),

          colourInput(
            ns("col_high"),
            "High",
            "#E34948"
          ),
          
          hr(),
          
          h4("8. Export"),
          
          fluidRow(
            column(
              6,
              numericInput(
                ns("exp_width"),
                "Width (cm)",
                value = 18,
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
            ns("download_heatmap"),
            "Download heatmap",
            class = "btn-success",
            style = "width:100%"
          ),
          
          br(),
          br(),
          
          downloadButton(
            ns("download_selected_matrix"),
            "Download selected z-score matrix",
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
          class = "heatmap-panel",
          h4("Input Summary"),
          verbatimTextOutput(ns("input_summary"))
        ),
        
        div(
          class = "heatmap-panel",
          h4("Expression Heatmap"),
          uiOutput(ns("heatmap_plot_ui"))
        ),
        
        div(
          class = "heatmap-panel",
          h4("Heatmap Matrix"),
          DTOutput(ns("selected_matrix_table"))
        )
      )
    )
  )
}

# =========================================================
# SERVER
# =========================================================

heatmapModuleServer <- function(id) {
  
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    # =====================================================
    # 1) Raw Data
    # =====================================================
    
    raw_data <- reactive({
      
      req(input$expr_file)
      
      hm_smart_read(
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
      
      label_txt <- if (isTRUE(input$transpose)) {
        "Sample ID column:"
      } else {
        "Gene/protein/feature ID column:"
      }
      
      selectInput(
        ns("id_col"),
        label_txt,
        choices = cols,
        selected = cols[1]
      )
    })
    
    # =====================================================
    # 3) Annotation UI
    # =====================================================
    
    output$annotation_ui <- renderUI({
      
      req(input$anno_file)
      
      anno <- hm_smart_read(
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
          ns("anno_group_col"),
          "Group / condition column:",
          choices = c("None" = "__none__", cols),
          selected = if (length(cols) >= 2) cols[2] else "__none__"
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
      
      if (isTRUE(input$transpose)) {
        
        sample_ids <- df[[input$id_col]] %>%
          as.character() %>%
          stringr::str_trim()
        
        expr_df <- df[, setdiff(colnames(df), input$id_col), drop = FALSE]
        
        expr_df <- omics_numeric_frame(expr_df)
        
        validate(
          need(
            ncol(expr_df) >= 2,
            "No numeric feature columns were found after removing the sample ID column."
          )
        )
        
        mat <- t(as.matrix(expr_df))
        colnames(mat) <- omics_validate_ids(sample_ids)
        feature_ids <- rownames(mat)
        
        feature_df <- as.data.frame(
          mat,
          check.names = FALSE
        )
        
        mat <- hm_collapse_duplicate_features_highest_mean(
          expr_df = feature_df,
          feature_ids = feature_ids
        )
        
      } else {
        
        feature_ids <- df[[input$id_col]]
        
        expr_df <- df[, setdiff(colnames(df), input$id_col), drop = FALSE]
        
        expr_df <- omics_numeric_frame(expr_df)
        
        validate(
          need(
            ncol(expr_df) >= 2,
            "No numeric sample columns were found after removing the feature ID column."
          )
        )
        
        mat <- hm_collapse_duplicate_features_highest_mean(
          expr_df = expr_df,
          feature_ids = feature_ids
        )
      }
      
      storage.mode(mat) <- "numeric"
      
      # ---------------------------------------------------
      # Remove rows with missing or non-finite values
      # ---------------------------------------------------
      
      missing_or_nonfinite_rows <- apply(mat, 1, function(x) {
        any(!is.finite(x))
      })
      
      removed_missing_rows <- sum(missing_or_nonfinite_rows)
      
      if (removed_missing_rows > 0) {
        mat <- mat[!missing_or_nonfinite_rows, , drop = FALSE]
      }
      
      attr(mat, "removed_missing_rows") <- removed_missing_rows
      
      validate(
        need(
          nrow(mat) >= 2,
          "Too few features remain after removing missing or non-finite values."
        ),
        need(
          ncol(mat) >= 2,
          "Heatmap requires at least 2 samples."
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
      
      if (isTRUE(input$transpose)) {
        ids <- setdiff(colnames(df), input$id_col)
      } else {
        ids <- df[[input$id_col]]
      }
      
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
      
      req(raw_expr_matrix())
      
      mat <- raw_expr_matrix()
      
      validate(
        need(
          nrow(mat) >= 2,
          "Too few features are available after duplicate handling."
        ),
        need(
          ncol(mat) >= 2,
          "Heatmap requires at least 2 samples."
        ),
        need(
          all(is.finite(mat)),
          "The processed matrix contains missing or non-finite values."
        )
      )
      
      mat
    })
    
    # =====================================================
    # 7) Manual Feature UI
    # =====================================================
    
    output$manual_feature_ui <- renderUI({
      
      req(processed_matrix())
      
      if (input$feature_selection != "manual") {
        return(NULL)
      }
      
      features <- rownames(processed_matrix())
      
      selectizeInput(
        ns("manual_features"),
        "Select features:",
        choices = features,
        selected = head(features, min(30, length(features))),
        multiple = TRUE,
        options = list(
          plugins = list("remove_button"),
          maxOptions = 5000
        )
      )
    })
    
    # =====================================================
    # 8) Uploaded Feature List
    # =====================================================
    
    uploaded_feature_list <- reactive({
      
      if (is.null(input$feature_list_file)) {
        return(character(0))
      }
      
      hm_read_feature_list(
        path = input$feature_list_file$datapath,
        original_name = input$feature_list_file$name
      )
    })
    
    # =====================================================
    # 9) Sample Annotation
    # =====================================================
    
    sample_annotation <- reactive({
      
      if (is.null(input$anno_file)) {
        return(NULL)
      }
      
      req(
        processed_matrix(),
        input$anno_sample_col,
        input$anno_group_col
      )
      
      anno <- hm_smart_read(
        path = input$anno_file$datapath,
        original_name = input$anno_file$name
      )
      
      validate(
        need(
          input$anno_sample_col %in% colnames(anno),
          "Selected sample column not found in annotation."
        )
      )
      
      samples <- colnames(processed_matrix())
      anno <- omics_annotation(anno, input$anno_sample_col, samples,
                               required = input$anno_group_col)

      
      if (input$anno_group_col == "__none__") {
        
        out <- tibble(
          Sample = samples,
          Group = "All samples"
        )
        
      } else {
        
        validate(
          need(
            input$anno_group_col %in% colnames(anno),
            "Selected group column not found in annotation."
          )
        )
        
        out <- anno %>%
          transmute(
            Sample = as.character(.data[[input$anno_sample_col]]),
            Group = as.character(.data[[input$anno_group_col]])
          ) %>%
          filter(
            !is.na(Sample),
            Sample != ""
          ) %>%
          distinct(
            Sample,
            .keep_all = TRUE
          )
        
        out <- tibble(Sample = samples) %>%
          left_join(
            out,
            by = "Sample"
          ) %>%
          mutate(
            Group = ifelse(
              is.na(Group) | Group == "",
              "Unannotated",
              Group
            )
          )
      }
      
      out$Sample <- factor(
        out$Sample,
        levels = samples
      )
      
      out$Group <- factor(
        out$Group,
        levels = unique(out$Group)
      )
      
      out
    })
    
    # =====================================================
    # 10) Selected Matrix
    # =====================================================
    
    selected_matrix_raw <- reactive({
      
      req(
        processed_matrix(),
        input$feature_selection,
        input$max_features_display
      )
      
      mat <- processed_matrix()
      
      if (input$feature_selection == "top_variable") {
        omics_number_parameter(input$top_variable_features, "Number of variable features", 2, .Machine$integer.max, TRUE)
        
        rv <- hm_row_var(mat)
        rv[!is.finite(rv)] <- NA_real_
        
        validate(
          need(
            any(!is.na(rv)),
            "No valid feature variance available."
          )
        )
        
        top_n <- min(
          input$top_variable_features,
          nrow(mat)
        )
        
        keep_idx <- order(
          rv,
          decreasing = TRUE,
          na.last = NA
        )[seq_len(top_n)]
        
        mat <- mat[keep_idx, , drop = FALSE]
      }
      
      if (input$feature_selection == "uploaded_list") {
        
        req(input$feature_list_file)
        
        features <- uploaded_feature_list()
        features_found <- features[features %in% rownames(mat)]
        
        validate(
          need(
            length(features_found) > 0,
            "None of the uploaded features were found in the processed matrix."
          )
        )
        
        mat <- mat[features_found, , drop = FALSE]
      }
      
      if (input$feature_selection == "manual") {
        
        req(input$manual_features)
        
        features <- input$manual_features
        features_found <- features[features %in% rownames(mat)]
        
        validate(
          need(
            length(features_found) > 0,
            "No manually selected features were found in the processed matrix."
          )
        )
        
        mat <- mat[features_found, , drop = FALSE]
      }
      
      if (input$feature_selection == "all") {
        mat <- mat
      }
      
      validate(
        need(
          nrow(mat) >= 2,
          "At least 2 features are required for the heatmap."
        ),
        need(
          nrow(mat) <= input$max_features_display,
          paste0(
            "Too many features selected (",
            nrow(mat),
            "). Increase the maximum allowed features or use top variable / uploaded feature list selection."
          )
        ),
        need(
          all(is.finite(mat)),
          "The selected heatmap matrix contains missing or non-finite values."
        )
      )
      
      mat
    })
    
    heatmap_matrix <- reactive({
      
      req(selected_matrix_raw())
      
      hm_scale_rows(selected_matrix_raw())
    })
    
    # =====================================================
    # 11) Row and Column Ordering
    # =====================================================
    
    row_order <- reactive({
      
      req(heatmap_matrix())
      
      mat <- heatmap_matrix()
      
      if (isTRUE(input$cluster_rows)) {
        
        hm_get_cluster_order(
          mat = mat,
          margin = "rows",
          distance = input$distance_method,
          method = input$hclust_method
        )
        
      } else {
        
        rownames(mat)
      }
    })
    
    col_order <- reactive({
      
      req(heatmap_matrix())
      
      mat <- heatmap_matrix()
      
      if (isTRUE(input$cluster_cols)) {
        
        hm_get_cluster_order(
          mat = mat,
          margin = "cols",
          distance = input$distance_method,
          method = input$hclust_method
        )
        
      } else {
        
        anno <- sample_annotation()
        
        if (
          !is.null(anno) &&
          isTRUE(input$order_samples_by_group) &&
          "Group" %in% colnames(anno)
        ) {
          
          anno %>%
            arrange(Group, Sample) %>%
            pull(Sample) %>%
            as.character()
          
        } else {
          
          colnames(mat)
        }
      }
    })
    
    # =====================================================
    # 12) Input Summary
    # =====================================================
    
    output$input_summary <- renderPrint({
      
      req(
        raw_expr_matrix(),
        processed_matrix(),
        selected_matrix_raw(),
        duplicate_summary()
      )
      
      raw_mat <- raw_expr_matrix()
      proc_mat <- processed_matrix()
      sel_mat <- selected_matrix_raw()
      dup <- duplicate_summary()
      anno <- sample_annotation()
      
      cat("Input matrix requirement: normalized/log2 matrix provided by user\n")
      cat("Duplicated feature handling: keep row with highest mean across samples\n")
      cat("Matrix transposed:", isTRUE(input$transpose), "\n")
      cat("Row z-score for heatmap visualization: applied automatically\n\n")
      
      cat("Raw features after duplicate handling:", nrow(raw_mat), "\n")
      cat(
        "Features removed due to missing/non-finite values:",
        attr(raw_mat, "removed_missing_rows") %||% 0,
        "\n"
      )
      cat("Processed features:", nrow(proc_mat), "\n")
      cat("Selected features for heatmap:", nrow(sel_mat), "\n")
      cat("Samples:", ncol(sel_mat), "\n\n")
      
      cat("Raw feature IDs:", dup$total_ids, "\n")
      cat("Unique raw feature IDs:", dup$unique_ids, "\n")
      cat("Duplicated rows:", dup$duplicated_rows, "\n")
      cat("Feature names with duplicates:", dup$duplicated_feature_names, "\n\n")
      
      cat("Feature selection:", input$feature_selection, "\n")
      cat("Cluster features:", input$cluster_rows, "\n")
      cat("Cluster samples:", input$cluster_cols, "\n\n")
      
      cat("Processed matrix range:\n")
      print(range(proc_mat, na.rm = TRUE))
      
      cat("\nHeatmap z-score matrix range:\n")
      print(range(heatmap_matrix(), na.rm = TRUE))
      
      cat("\n")
      
      if (is.null(anno)) {
        cat("Sample annotation: not provided\n")
      } else {
        cat("Sample annotation: provided\n")
        print(table(anno$Group))
      }
      
      if (input$feature_selection == "uploaded_list") {
        cat("\nUploaded feature list size:", length(uploaded_feature_list()), "\n")
        cat("Features found:", sum(uploaded_feature_list() %in% rownames(proc_mat)), "\n")
      }
    })
    
    # =====================================================
    # 13) Heatmap Plot
    # =====================================================
    
    create_heatmap_plot <- reactive({
      
      req(
        heatmap_matrix(),
        row_order(),
        col_order()
      )
      
      mat <- heatmap_matrix()
      
      rows <- row_order()
      cols <- col_order()
      
      rows <- rows[rows %in% rownames(mat)]
      cols <- cols[cols %in% colnames(mat)]
      
      mat <- mat[rows, cols, drop = FALSE]
      
      plot_df <- as.data.frame(mat) %>%
        rownames_to_column("Feature") %>%
        pivot_longer(
          cols = -Feature,
          names_to = "Sample",
          values_to = "Expression"
        )
      
      plot_df$Feature <- factor(
        plot_df$Feature,
        levels = rev(rows)
      )
      
      plot_df$Sample <- factor(
        plot_df$Sample,
        levels = cols
      )
      
      anno <- sample_annotation()
      
      if (!is.null(anno)) {
        
        plot_df <- plot_df %>%
          left_join(
            anno %>%
              mutate(Sample = as.character(Sample)),
            by = "Sample"
          )
        
        plot_df$Group <- factor(
          plot_df$Group,
          levels = unique(as.character(anno$Group))
        )
      }
      
      fill_scale <- scale_fill_gradient2(
        low = input$col_low,
        mid = input$col_mid,
        high = input$col_high,
        midpoint = 0,
        name = "Row z-score"
      )
      
      x_text_theme <- if (isTRUE(input$show_sample_names)) {
        element_text(
          angle = 60,
          hjust = 1,
          vjust = 1,
          size = input$sample_label_size,
          colour = "black"
        )
      } else {
        element_blank()
      }
      
      y_text_theme <- if (isTRUE(input$show_feature_names)) {
        element_text(
          size = input$feature_label_size,
          colour = "black"
        )
      } else {
        element_blank()
      }
      
      tile_color <- if (isTRUE(input$tile_border)) {
        "grey90"
      } else {
        NA
      }
      
      base_plot <- ggplot(
        plot_df,
        aes(
          x = Sample,
          y = Feature,
          fill = Expression
        )
      ) +
        geom_tile(
          color = tile_color,
          linewidth = 0.1
        ) +
        fill_scale +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme_bw(base_size = input$font_size) +
        theme(
          axis.title = element_blank(),
          axis.text.x = x_text_theme,
          axis.text.y = y_text_theme,
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = input$font_size + 2, face = "bold", color = "black", hjust = 0.5),
          legend.title = element_text(size = input$font_size, color = "black"),
          legend.text = element_text(size = input$font_size, color = "black"),
          strip.background = element_rect(fill = "grey95", color = "grey70"),
          strip.text = element_text(size = input$font_size, face = "bold", color = "black")
        ) +
        labs(
          title = input$plot_title
        )
      
      if (
        !is.null(anno) &&
        isTRUE(input$facet_by_group) &&
        "Group" %in% colnames(plot_df)
      ) {
        
        base_plot +
          facet_grid(
            . ~ Group,
            scales = "free_x",
            space = "free_x"
          )
        
      } else {
        
        base_plot
      }
    })
    
    output$heatmap_plot_ui <- renderUI({
      
      req(selected_matrix_raw())
      
      n_features <- nrow(selected_matrix_raw())
      plot_height <- max(500, min(1400, 18 * n_features + 180))
      
      plotOutput(
        ns("heatmap_plot"),
        height = paste0(plot_height, "px")
      )
    })
    
    output$heatmap_plot <- renderPlot({
      
      req(create_heatmap_plot())
      
      print(create_heatmap_plot())
      
    }, res = 120)
    
    # =====================================================
    # 14) Tables
    # =====================================================
    
    selected_matrix_df <- reactive({
      
      req(heatmap_matrix())
      
      as.data.frame(heatmap_matrix()) %>%
        rownames_to_column("Feature")
    })
    
    processed_matrix_df <- reactive({
      
      req(processed_matrix())
      
      as.data.frame(processed_matrix()) %>%
        rownames_to_column("Feature")
    })
    
    output$selected_matrix_table <- renderDT({
      
      req(selected_matrix_df())
      
      datatable(
        selected_matrix_df(),
        rownames = FALSE,
        options = list(
          pageLength = 10,
          scrollX = TRUE
        )
      )
    })
    
    # =====================================================
    # 15) Downloads
    # =====================================================
    
    output$download_heatmap <- downloadHandler(
      
      filename = function() {
        paste0(
          "expression_heatmap_",
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

        hm_open_plot_device(
          file = file,
          fmt = input$exp_fmt,
          width = width_in,
          height = height_in,
          dpi = input$exp_dpi
        )

        print(create_heatmap_plot())
        dev.off()
      }
    )

    output$download_selected_matrix <- downloadHandler(

      filename = function() {
        paste0("expression_heatmap_selected_zscore_matrix_", Sys.Date(), ".txt")
      },

      content = function(file) {
        write.table(
          selected_matrix_df(),
          file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )
      }
    )

    output$download_processed_matrix <- downloadHandler(

      filename = function() {
        paste0("expression_heatmap_processed_matrix_", Sys.Date(), ".txt")
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
        selected_matrix_raw = selected_matrix_raw,
        heatmap_matrix = heatmap_matrix,
        sample_annotation = sample_annotation
      )
    )
  })
}