# =========================================================
# EzOmicsViz — Volcano Plot Module
# =========================================================

library(shiny)
library(EnhancedVolcano)
library(tidyverse)
library(colourpicker)
library(DT)
library(svglite)
library(readxl)

# =========================================================
# Helper Functions
# =========================================================

guess_column <- function(cols, pattern) {
  
  matched_cols <- grep(
    pattern = pattern,
    x = cols,
    value = TRUE,
    ignore.case = TRUE
  )
  
  if (length(matched_cols) == 0) {
    cols[1]
  } else {
    matched_cols[1]
  }
}

sanitize_filename <- function(x) {
  
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  
  if (is.na(x) || x == "") {
    x <- "selected_column"
  }
  
  x
}

is_adjusted_p_column <- function(col_name) {
  
  if (is.null(col_name) || length(col_name) == 0 || is.na(col_name)) {
    return(FALSE)
  }
  
  grepl(
    pattern = "padj|adj|adjust|adjusted|fdr|qvalue|q.value|qval|q_val",
    x = col_name,
    ignore.case = TRUE
  )
}

p_value_axis_label <- function(col_name) {

  if (is_adjusted_p_column(col_name)) {
    expression(-Log[10]~"Adjusted"~P)
  } else {
    expression(-Log[10]~P)
  }
}

is_logfc_column <- function(col_name) {

  if (is.null(col_name) || length(col_name) == 0 || is.na(col_name)) {
    return(FALSE)
  }

  grepl(
    pattern = "log2foldchange|log2fc|logfc|avg_log2fc|^fc$",
    x = col_name,
    ignore.case = TRUE
  )
}

logfc_axis_label <- function(col_name) {

  if (is_logfc_column(col_name)) {
    expression(Log[2]~"fold change")
  } else {
    col_name
  }
}

p_value_export_label <- function(col_name) {
  
  if (is_adjusted_p_column(col_name)) {
    "adjusted_p_value"
  } else {
    "p_value"
  }
}

as_numeric_flexible <- omics_numeric

clean_column_names <- function(df) {
  
  colnames(df) <- colnames(df) %>%
    as.character() %>%
    stringr::str_trim()
  
  df
}

try_read_delim_safe <- function(path, delim, decimal_mark = ".") {
  
  out <- tryCatch(
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
        clean_column_names()
    },
    error = function(e) {
      NULL
    }
  )
  
  out
}

try_read_whitespace_safe <- function(path, decimal_mark = ".") {
  
  out <- tryCatch(
    {
      readr::read_table(
        file = path,
        show_col_types = FALSE,
        col_types = readr::cols(.default = readr::col_character()),
        trim_ws = TRUE,
        locale = readr::locale(decimal_mark = decimal_mark)
      ) %>%
        omics_delimited_frame() %>%
        clean_column_names()
    },
    error = function(e) {
      NULL
    }
  )
  
  out
}

try_read_excel_safe <- function(path) {
  
  out <- tryCatch(
    {
      readxl::read_excel(
        path = path,
        sheet = 1,
          .name_repair = "minimal"
      ) %>%
        as.data.frame() %>%
        clean_column_names()
    },
    error = function(e) {
      NULL
    }
  )
  
  out
}

score_read_table <- function(df) {
  
  if (is.null(df)) {
    return(-Inf)
  }
  
  if (!is.data.frame(df)) {
    return(-Inf)
  }
  
  if (ncol(df) <= 0) {
    return(-Inf)
  }
  
  non_empty_names <- sum(!is.na(colnames(df)) & colnames(df) != "")
  
  ncol(df) * 1000 + non_empty_names
}

read_uploaded_table_unchecked <- function(path, original_name = NULL) {
  
  ext_source <- original_name %||% path
  ext <- tolower(tools::file_ext(ext_source))
  
  if (ext %in% c("xls", "xlsx")) {
    
    validate(
      need(
        requireNamespace("readxl", quietly = TRUE),
        "Package 'readxl' is required to read Excel files."
      )
    )
    
    df <- try_read_excel_safe(path)
    
    validate(
      need(
        !is.null(df) && ncol(df) > 1,
        paste(
          "The uploaded Excel file could not be parsed into multiple columns.",
          "Use the first sheet with a header row and at least two columns."
        )
      )
    )
    
  } else {
    
    read_attempts <- list(
      comma_dot = try_read_delim_safe(
        path = path,
        delim = ",",
        decimal_mark = "."
      ),
      
      semicolon_dot = try_read_delim_safe(
        path = path,
        delim = ";",
        decimal_mark = "."
      ),
      
      semicolon_comma = try_read_delim_safe(
        path = path,
        delim = ";",
        decimal_mark = ","
      ),
      
      tab_dot = try_read_delim_safe(
        path = path,
        delim = "\t",
        decimal_mark = "."
      ),
      
      tab_comma = try_read_delim_safe(
        path = path,
        delim = "\t",
        decimal_mark = ","
      ),
      
      whitespace_dot = try_read_whitespace_safe(
        path = path,
        decimal_mark = "."
      ),
      
      whitespace_comma = try_read_whitespace_safe(
        path = path,
        decimal_mark = ","
      )
    )
    
    scores <- vapply(
      read_attempts,
      score_read_table,
      numeric(1)
    )
    
    best_name <- names(which.max(scores))
    df <- read_attempts[[best_name]]
    
    validate(
      need(
        !is.null(df) && ncol(df) > 1,
        paste(
          "The uploaded file could not be parsed into multiple columns.",
          "Use a file separated by comma, semicolon, tab, or spaces, with a header row."
        )
      )
    )
  }
  
  df <- clean_column_names(df)
  
  if (
    colnames(df)[1] == "X" ||
    is.na(colnames(df)[1]) ||
    colnames(df)[1] == ""
  ) {
    
    df <- df %>%
      rownames_to_column(var = "gene_id_internal")
    
  } else {
    
    df$gene_id_internal <- df[[1]]
  }
  
  df
}

read_uploaded_table <- function(...) {
  omics_validate_table(read_uploaded_table_unchecked(...))
}

# =========================================================
# UI
# =========================================================

volcanoModuleUI <- function(id) {
  
  ns <- NS(id)
  
  tagList(
    
    tags$head(
      tags$style(HTML("
        .volcano-panel {
          background: #ffffff;
          border-radius: 8px;
          padding: 20px;
          border: 1px solid #eeeeee;
          margin-bottom: 20px;
        }

        .volcano-panel-centered {
          background: #ffffff;
          border-radius: 8px;
          padding: 20px;
          border: 1px solid #eeeeee;
          margin-bottom: 20px;
          text-align: center;
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
            inputId = ns("upload_file"),
            label   = "Upload File (.txt, .csv, .tsv, .xls, .xlsx)",
            accept  = c(".txt", ".csv", ".tsv", ".xls", ".xlsx")
          ),
          
          uiOutput(ns("column_selector")),
          
          hr(),
          
          h4("2. Highlight Selection"),
          
          numericInput(
            inputId = ns("top_n"),
            label   = "Show Top N genes/proteins per direction:",
            value   = 5,
            min     = 1,
            step    = 1
          ),
          
          hr(),
          
          h4("3. Thresholds"),
          
          numericInput(
            inputId = ns("p_cutoff"),
            label   = "P-value / Adjusted P-value cutoff:",
            value   = 0.05,
            min     = 0,
            max     = 1,
            step    = 0.01
          ),
          
          numericInput(
            inputId = ns("fc_cutoff"),
            label   = "Log2FC cutoff:",
            value   = 1.0,
            min     = 0,
            step    = 0.1
          ),
          
          hr(),
          
          h4("4. Axis Control"),
          
          sliderInput(
            inputId = ns("range_x"),
            label   = "X-axis range (log2FC):",
            min     = -10,
            max     = 10,
            value   = c(-5, 5),
            step    = 0.5
          ),
          
          sliderInput(
            inputId = ns("range_y"),
            label   = "Y-axis range (-log10 P / Adjusted P):",
            min     = 0,
            max     = 100,
            value   = c(0, 20),
            step    = 1
          ),
          
          hr(),
          
          h4("5. Appearance"),
          
          textInput(
            inputId = ns("plot_title"),
            label   = "Plot Title",
            value   = "Volcano plot"
          ),
          
          sliderInput(
            inputId = ns("point_size"),
            label   = "Point Size",
            min     = 0.5,
            max     = 5,
            value   = 2.0,
            step    = 0.5
          ),
          
          sliderInput(
            inputId = ns("lab_size"),
            label   = "Gene Label Size",
            min     = 1,
            max     = 8,
            value   = 3.5,
            step    = 0.5
          ),

          sliderInput(
            inputId = ns("font_size"),
            label   = "Font size:",
            min     = 8,
            max     = 18,
            value   = 12,
            step    = 1
          ),

          selectInput(
            inputId = ns("legend_position"),
            label   = "Legend Position",
            choices = c(
              "Right"  = "right",
              "Bottom" = "bottom",
              "Top"    = "top",
              "Left"   = "left",
              "None"   = "none"
            ),
            selected = "right"
          ),
          
          hr(),
          
          h4("6. Colors"),
          
          fluidRow(
            column(
              width = 6,
              colourpicker::colourInput(
                inputId = ns("col_up"),
                label   = "Upregulated",
                value   = "#E34948"
              )
            ),
            
            column(
              width = 6,
              colourpicker::colourInput(
                inputId = ns("col_down"),
                label   = "Downregulated",
                value   = "#2A78D6"
              )
            )
          ),
          
          colourpicker::colourInput(
            inputId = ns("col_ns"),
            label   = "NS",
            value   = "#898781"
          ),
          
          hr(),
          
          h4("7. Export"),
          
          fluidRow(
            column(
              width = 6,
              numericInput(
                inputId = ns("exp_width"),
                label   = "Width (cm)",
                value   = 14,
                min     = 5,
                max     = 60
              )
            ),

            column(
              width = 6,
              numericInput(
                inputId = ns("exp_height"),
                label   = "Height (cm)",
                value   = 10,
                min     = 5,
                max     = 60
              )
            )
          ),

          fluidRow(
            column(
              width = 6,
              selectInput(
                inputId = ns("exp_fmt"),
                label   = "Format",
                choices = c("pdf", "png", "tiff", "svg"),
                selected = "svg"
              )
            ),

            column(
              width = 6,
              conditionalPanel(
                condition = sprintf("input['%s'] == 'png' || input['%s'] == 'tiff'", ns("exp_fmt"), ns("exp_fmt")),
                numericInput(
                  inputId = ns("exp_dpi"),
                  label   = "DPI",
                  value   = 300,
                  min     = 72,
                  max     = 1200,
                  step    = 50
                )
              )
            )
          ),
          
          downloadButton(
            outputId = ns("download_volcano"),
            label    = "Download Image",
            class    = "btn-success",
            style    = "width:100%;"
          ),
          
          br(),
          br(),
          
          downloadButton(
            outputId = ns("download_sig"),
            label    = "Download Significant Features",
            class    = "btn-info",
            style    = "width:100%;"
          )
        )
      ),
      
      column(
        width = 8,
        
        div(
          class = "volcano-panel-centered",
          plotOutput(
            outputId = ns("volcano_plot"),
            height   = "500px"
          )
        ),
        
        div(
          class = "volcano-panel",
          h4("Significant Results Table"),
          DTOutput(ns("sig_table"))
        )
      )
    )
  )
}

# =========================================================
# SERVER
# =========================================================

volcanoModuleServer <- function(id) {
  
  moduleServer(id, function(input, output, session) {
    
    # =====================================================
    # 1) Read Raw Data
    # =====================================================
    
    raw_data <- reactive({
      
      req(input$upload_file)
      
      read_uploaded_table(
        path = input$upload_file$datapath,
        original_name = input$upload_file$name
      )
    })
    
    # =====================================================
    # 2) Column Selectors
    # =====================================================
    
    output$column_selector <- renderUI({
      
      df <- raw_data()
      cols <- colnames(df)
      
      tagList(
        selectInput(
          inputId  = session$ns("col_lab"),
          label    = "Label Column:",
          choices  = cols,
          selected = guess_column(
            cols,
            "protein|gene|symbol|id|name"
          )
        ),
        
        selectInput(
          inputId  = session$ns("col_logfc"),
          label    = "Log2FC Column:",
          choices  = cols,
          selected = guess_column(
            cols,
            "log2FoldChange|logFC|log2FC|avg_log2FC|FC"
          )
        ),
        
        selectInput(
          inputId  = session$ns("col_p"),
          label    = "P-value / Adjusted P-value Column:",
          choices  = cols,
          selected = guess_column(
            cols,
            "adj.P.Val|padj|adj_p|adjusted|fdr|qvalue|P.Value|pvalue|p.value"
          )
        )
      )
    })
    
    # =====================================================
    # 3) Prepared Data
    # =====================================================
    
    prepared_data <- reactive({
      
      req(input$col_lab, input$col_logfc, input$col_p)
      
      df <- raw_data()
      
      validate(
        need(input$col_logfc %in% colnames(df), "Log2FC column not found."),
        need(input$col_p %in% colnames(df), "P-value column not found."),
        need(input$col_lab %in% colnames(df), "Label column not found."),
        need(
          length(unique(c(input$col_lab, input$col_logfc, input$col_p))) == 3,
          "Label, Log2FC, and P-value columns must all be different."
        )
      )
      
      df[[input$col_logfc]] <- as_numeric_flexible(df[[input$col_logfc]])
      df[[input$col_p]]     <- as_numeric_flexible(df[[input$col_p]])
      df[[input$col_lab]]   <- as.character(df[[input$col_lab]])
      omics_excluded_rows(sum(is.na(df[[input$col_logfc]]) | is.na(df[[input$col_p]])), "Volcano input")
      omics_probability(df[[input$col_p]])
      if (any(df[[input$col_p]] == 0, na.rm = TRUE)) {
        showNotification("Zero P-values use a finite numerical floor for drawing. Exported values remain zero.",
                         type = "message", duration = 10)
      }
      validate(need(all(is.na(df[[input$col_logfc]]) | is.finite(df[[input$col_logfc]])),
                    "Log2FC values must be finite."))
      
      df <- df %>%
        filter(
          !is.na(.data[[input$col_logfc]]),
          !is.na(.data[[input$col_p]])
        )
      
      validate(
        need(
          nrow(df) > 0,
          paste0(
            "No valid rows after reading numeric values. Check whether '",
            input$col_logfc,
            "' and '",
            input$col_p,
            "' are the correct columns, and whether the file uses comma, semicolon, tab, whitespace, or a valid Excel sheet."
          )
        )
      )
      
      df
    })
    
    # =====================================================
    # 4) Dynamic Axis Adjustment
    # =====================================================
    
    volcano_thresholds <- reactive({
      omics_number_parameter(input$p_cutoff, "P-value cutoff", 0, 1)
      omics_number_parameter(input$fc_cutoff, "Log2FC cutoff", 0, Inf)
      TRUE
    })

    observe({
      
      req(input$upload_file, input$col_logfc, input$col_p)
      
      df <- prepared_data()
      
      fc_val <- df[[input$col_logfc]]
      p_val  <- df[[input$col_p]]
      
      fc_min <- floor(min(fc_val, na.rm = TRUE)) - 1
      fc_max <- ceiling(max(fc_val, na.rm = TRUE)) + 1
      
      p_safe <- omics_plot_probabilities(p_val)
      y_max  <- ceiling(max(-log10(p_safe), na.rm = TRUE)) + 2
      
      updateSliderInput(
        session = session,
        inputId = "range_x",
        min     = fc_min,
        max     = fc_max,
        value   = c(fc_min, fc_max)
      )
      
      updateSliderInput(
        session = session,
        inputId = "range_y",
        min     = 0,
        max     = y_max,
        value   = c(0, y_max)
      )
    })
    
    # =====================================================
    # 5) Significant Data
    # =====================================================
    
    sig_data <- reactive({
      req(volcano_thresholds())
      
      df <- prepared_data()
      
      df %>%
        filter(
          .data[[input$col_p]] <= input$p_cutoff,
          abs(.data[[input$col_logfc]]) >= input$fc_cutoff
        )
    })
    
    sig_export_data <- reactive({
      
      sig_data() %>%
        mutate(
          selected_p_column = input$col_p,
          selected_p_metric = p_value_export_label(input$col_p),
          p_cutoff_used = input$p_cutoff,
          log2fc_cutoff_used = input$fc_cutoff,
          .before = 1
        )
    })
    
    # =====================================================
    # 6) Labels To Highlight
    # =====================================================
    
    selected_labels <- reactive({
      req(volcano_thresholds())
      omics_number_parameter(input$top_n, "Number of labels per direction", 1, .Machine$integer.max, TRUE)
      
      df <- prepared_data()
      
      df_sig <- df %>%
        filter(.data[[input$col_p]] <= input$p_cutoff)
      
      if (nrow(df_sig) == 0) {
        return(character())
      }
      
      top_down <- df_sig %>%
        arrange(.data[[input$col_logfc]]) %>%
        slice_head(n = input$top_n)
      
      top_up <- df_sig %>%
        arrange(desc(.data[[input$col_logfc]])) %>%
        slice_head(n = input$top_n)
      
      unique(c(
        top_down[[input$col_lab]],
        top_up[[input$col_lab]]
      ))
    })
    
    # =====================================================
    # 7) Custom Color Vector
    # =====================================================
    
    custom_colors <- reactive({
      req(volcano_thresholds())
      
      df <- prepared_data()
      
      regulation <- case_when(
        df[[input$col_logfc]] <= -input$fc_cutoff &
          df[[input$col_p]] <= input$p_cutoff ~ "Downregulated",
        
        df[[input$col_logfc]] >= input$fc_cutoff &
          df[[input$col_p]] <= input$p_cutoff ~ "Upregulated",
        
        TRUE ~ "NS"
      )
      
      keyvals <- case_when(
        regulation == "Downregulated" ~ input$col_down,
        regulation == "Upregulated"   ~ input$col_up,
        TRUE                          ~ input$col_ns
      )
      
      names(keyvals) <- regulation
      
      keyvals
    })
    
    # =====================================================
    # 8) Volcano Plot
    # =====================================================
    
    create_volcano <- reactive({
      req(volcano_thresholds())
      
      req(
        input$col_lab,
        input$col_logfc,
        input$col_p,
        input$range_x,
        input$range_y
      )
      
      df <- prepared_data()
      df[[input$col_p]] <- omics_plot_probabilities(df[[input$col_p]])
      
      legend_position <- input$legend_position %||% "right"
      
      p <- EnhancedVolcano(
        toptable = df,
        lab      = df[[input$col_lab]],
        x        = input$col_logfc,
        y        = input$col_p,
        
        selectLab = selected_labels(),
        
        title    = input$plot_title,
        subtitle = NULL,
        caption  = NULL,
        
        xlab = logfc_axis_label(input$col_logfc),
        ylab = p_value_axis_label(input$col_p),
        
        pCutoff  = input$p_cutoff,
        FCcutoff = input$fc_cutoff,
        
        xlim = input$range_x,
        ylim = input$range_y,
        
        pointSize = input$point_size,
        labSize   = input$lab_size,
        
        colCustom = custom_colors(),
        colAlpha  = 0.8,
        
        boxedLabels = TRUE,
        
        drawConnectors  = TRUE,
        widthConnectors = 0.4,
        colConnectors   = "black",
        
        gridlines.major = FALSE,
        gridlines.minor = FALSE,
        
        border       = "full",
        borderWidth  = 0.8,
        borderColour = "black",
        
        legendPosition = legend_position,
        legendLabSize  = input$font_size,
        legendIconSize = 4
      )
      
      p +
        guides(
          colour = guide_legend(
            title = NULL,
            override.aes = list(
              size = 4
            )
          )
        ) +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          
          axis.text = element_text(
            color = "black",
            size  = input$font_size
          ),

          axis.title = element_text(
            color = "black",
            size  = input$font_size
          ),

          axis.ticks = element_line(
            color = "black"
          ),

          plot.title = element_text(
            color = "black",
            face  = "bold",
            size  = input$font_size + 2,
            hjust = 0.5
          ),

          legend.position = legend_position,
          legend.title = element_blank(),
          legend.text = element_text(
            color = "black",
            size  = input$font_size
          ),
          legend.key.size = grid::unit(
            4,
            "mm"
          )
        )
    })
    
    # =====================================================
    # 9) Render Plot
    # =====================================================
    
    output$volcano_plot <- renderPlot(
      expr = {
        print(create_volcano())
      },
      height = 500,
      width = 500
    )
    
    # =====================================================
    # 10) Significant Results Table
    # =====================================================
    
    output$sig_table <- renderDT({
      
      datatable(
        sig_data(),
        rownames = FALSE,
        options = list(
          pageLength = 10,
          scrollX    = TRUE
        )
      )
    })
    
    # =====================================================
    # 11) Download Significant Results
    # =====================================================
    
    output$download_sig <- downloadHandler(
      
      filename = function() {
        
        paste0(
          "significant_features_",
          p_value_export_label(input$col_p),
          "_",
          sanitize_filename(input$col_p),
          "_",
          Sys.Date(),
          ".txt"
        )
      },
      
      content = function(file) {
        
        readr::write_delim(
          x     = sig_export_data(),
          file  = file,
          delim = "\t"
        )
      }
    )
    
    # =====================================================
    # 12) Download Volcano Plot
    # =====================================================
    
    output$download_volcano <- downloadHandler(
      
      filename = function() {
        paste0("volcano_plot_", Sys.Date(), "_", input$exp_width, "x", input$exp_height, ".", input$exp_fmt)
      },
      
      content = function(file) {
        
        width_in  <- input$exp_width * 0.393701
        height_in <- input$exp_height * 0.393701
        
        if (input$exp_fmt == "png") {
          
          png(
            filename = file,
            width    = width_in,
            height   = height_in,
            units    = "in",
            res      = input$exp_dpi
          )
          
        } else if (input$exp_fmt == "pdf") {
          
          pdf(
            file   = file,
            width  = width_in,
            height = height_in
          )
          
        } else if (input$exp_fmt == "tiff") {
          
          tiff(
            filename    = file,
            width       = width_in,
            height      = height_in,
            units       = "in",
            res         = input$exp_dpi,
            compression = "lzw"
          )
          
        } else if (input$exp_fmt == "svg") {
          
          svglite::svglite(
            file   = file,
            width  = width_in,
            height = height_in
          )
        }
        
        print(create_volcano())
        dev.off()
      }
    )
  })
}
