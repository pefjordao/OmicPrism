# =========================================================
# EzOmicsViz — Overlap Module
# =========================================================

library(shiny)
library(eulerr)
library(tidyverse)
library(shinyWidgets)
library(colourpicker)
library(svglite)
library(DT)
library(cowplot)
library(grid)
library(readxl)

# =========================================================
# Helper Functions
# =========================================================

open_venn_device <- function(file, fmt, width, height, dpi) {
  
  if (fmt == "png") {
    
    png(
      filename = file,
      width    = width,
      height   = height,
      units    = "in",
      res      = dpi
    )
    
  } else if (fmt == "pdf") {
    
    pdf(
      file   = file,
      width  = width,
      height = height
    )
    
  } else if (fmt == "tiff") {
    
    tiff(
      filename    = file,
      width       = width,
      height      = height,
      units       = "in",
      res         = dpi,
      compression = "lzw"
    )
    
  } else if (fmt == "svg") {
    
    svglite::svglite(
      file   = file,
      width  = width,
      height = height
    )
  }
}

format_intersection_label <- function(active_sets, all_sets) {
  
  active_sets <- active_sets[!is.na(active_sets) & active_sets != ""]
  
  if (length(active_sets) == 0) {
    return("None")
  }
  
  label <- paste(active_sets, collapse = " & ")
  
  if (length(active_sets) < length(all_sets)) {
    label <- paste(label, "only")
  }
  
  label
}

blend_colors <- function(colors) {
  if (length(colors) == 1) return(colors)
  rgb_vals <- col2rgb(colors)
  avg_rgb <- rowMeans(rgb_vals)
  rgb(avg_rgb[1], avg_rgb[2], avg_rgb[3], maxColorValue = 255)
}

default_group_colors <- function(n) {
  
  base_cols <- c(
    "#2A78D6",
    "#EB6834",
    "#1BAF7A",
    "#EDA100",
    "#E87BA4"
  )

  rep(base_cols, length.out = n)
}

# =========================================================
# UI
# =========================================================

vennModuleUI <- function(id) {
  
  ns <- NS(id)
  
  tagList(
    
    tags$head(
      tags$style(HTML("
        .set-block {
          background: #ffffff;
          border: 1px solid #dee2e6;
          border-radius: 8px;
          padding: 15px;
          margin-bottom: 15px;
          box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }

        .set-header {
          display: flex;
          gap: 10px;
          align-items: center;
          margin-bottom: 8px;
        }

        .plot-container {
          background: #ffffff;
          border-radius: 8px;
          padding: 20px;
          text-align: center;
          border: 1px solid #eeeeee;
          margin-bottom: 20px;
        }

        .table-container {
          background: #ffffff;
          border-radius: 8px;
          padding: 20px;
          border: 1px solid #eeeeee;
          margin-bottom: 20px;
        }

        .nav-tabs {
          margin-bottom: 15px;
        }

        .well {
          background-color: #ffffff;
          border: 1px solid #e3e3df;
          border-radius: 8px;
          box-shadow: none;
        }

        .block-help {
          color: #6b7280;
          font-size: 0.90rem;
          margin-bottom: 15px;
        }

        table.dataTable td.dt-intersection-cell {
          white-space: normal !important;
          word-break: break-word;
          min-width: 260px;
          max-width: 520px;
        }

        table.dataTable td.dt-center-cell {
          text-align: center;
          vertical-align: middle;
        }
      "))
    ),
    
    fluidRow(
      
      column(
        width = 4,
        
        wellPanel(
          h4("1. Diagram Controls"),
          
          sliderInput(
            inputId = ns("n_sets"),
            label   = "Number of Groups:",
            min     = 2,
            max     = 5,
            value   = 2,
            step    = 1
          ),
          
          hr(),
          
          radioGroupButtons(
            inputId  = ns("plot_type"),
            label    = "Main Visualization",
            choices  = c("Venn", "Euler", "UpSet"),
            selected = "Venn",
            justified = TRUE
          ),
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'UpSet'", ns("plot_type")),
            
            h5("UpSet Options"),
            
            numericInput(
              inputId = ns("top_intersections"),
              label   = "Top intersections:",
              value   = 20,
              min     = 1,
              step    = 1
            ),
            
            numericInput(
              inputId = ns("min_intersection_size"),
              label   = "Minimum intersection size:",
              value   = 1,
              min     = 1,
              step    = 1
            ),
            
            selectInput(
              inputId = ns("upset_order"),
              label   = "Order intersections by:",
              choices = c(
                "Intersection size" = "count",
                "Intersection degree" = "degree",
                "Intersection name" = "name"
              ),
              selected = "count"
            ),
            
            h5("UpSet Colors"),
            
            colourInput(
              inputId = ns("upset_bar_color"),
              label   = "Intersection bar color:",
              value   = "#2A78D6"
            ),
            
            colourInput(
              inputId = ns("upset_active_dot_color"),
              label   = "Active dot color:",
              value   = "#222222"
            ),
            
            colourInput(
              inputId = ns("upset_inactive_dot_color"),
              label   = "Inactive dot color:",
              value   = "#d9d9d9"
            ),
            
            colourInput(
              inputId = ns("upset_line_color"),
              label   = "Connection line color:",
              value   = "#333333"
            )
          ),
          
          h5("Appearance"),

          textInput(
            inputId = ns("plot_title"),
            label   = "Plot Title (Main Plot)",
            value   = ""
          ),

          textInput(
            inputId = ns("size_plot_title"),
            label   = "Plot Title (Set Size Barplot)",
            value   = "Set size summary"
          ),

          sliderInput(
            inputId = ns("font_size"),
            label   = "Font size:",
            min     = 8,
            max     = 18,
            value   = 12,
            step    = 1
          ),
          
          sliderInput(
            inputId = ns("opacity"),
            label   = "Opacity",
            min     = 0.1,
            max     = 1,
            value   = 0.5,
            step    = 0.1
          )
        ),
        
        wellPanel(
          h4("2. Data Input"),
          
          p(
            class = "block-help",
            "Enter one gene, protein, or ID per line, or upload a file for each group (.txt, .csv, .tsv, .xls, .xlsx)."
          ),
          
          uiOutput(ns("individual_sets_ui"))
        ),
        
        wellPanel(
          h4("3. Export"),

          h5("Main Plot"),

          fluidRow(
            column(
              width = 6,
              numericInput(
                inputId = ns("exp_width"),
                label   = "Width (cm)",
                value   = 15,
                min     = 5,
                max     = 60,
                step    = 1
              )
            ),

            column(
              width = 6,
              numericInput(
                inputId = ns("exp_height"),
                label   = "Height (cm)",
                value   = 11,
                min     = 5,
                max     = 60,
                step    = 1
              )
            )
          ),

          h5("Set Size Barplot"),

          fluidRow(
            column(
              width = 6,
              numericInput(
                inputId = ns("exp_width_size"),
                label   = "Width (cm)",
                value   = 14,
                min     = 5,
                max     = 60,
                step    = 1
              )
            ),

            column(
              width = 6,
              numericInput(
                inputId = ns("exp_height_size"),
                label   = "Height (cm)",
                value   = 10,
                min     = 5,
                max     = 60,
                step    = 1
              )
            )
          ),

          fluidRow(
            column(
              width = 6,
              selectInput(
                inputId = ns("export_format"),
                label   = "Format",
                choices = c("pdf", "png", "tiff", "svg"),
                selected = "svg"
              )
            ),

            column(
              width = 6,
              conditionalPanel(
                condition = sprintf("input['%s'] == 'png' || input['%s'] == 'tiff'", ns("export_format"), ns("export_format")),
                numericInput(
                  inputId = ns("export_dpi"),
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
            outputId = ns("download_plot"),
            label    = "Download Main Plot",
            class    = "btn-success",
            style    = "width:100%;"
          ),
          
          br(),
          br(),
          
          downloadButton(
            outputId = ns("download_size_plot"),
            label    = "Download Set Size Plot",
            class    = "btn-success",
            style    = "width:100%;"
          ),
          
          br(),
          br(),
          
          downloadButton(
            outputId = ns("download_data"),
            label    = "Download All Intersections",
            class    = "btn-info",
            style    = "width:100%;"
          )
        )
      ),
      
      column(
        width = 8,
        
        div(
          class = "plot-container",
          tabsetPanel(
            tabPanel(
              title = "Main Plot",
              br(),
              plotOutput(
                outputId = ns("main_plot"),
                inline = TRUE
              )
            ),
            
            tabPanel(
              title = "Set Size Barplot",
              br(),
              plotOutput(
                outputId = ns("set_size_plot"),
                height = "520px"
              )
            )
          )
        ),
        
        div(
          class = "table-container",
          h4("Intersection Summary"),
          DTOutput(ns("intersection_summary_table"))
        ),
        
        div(
          class = "table-container",
          h4("Intersections Explorer"),
          uiOutput(ns("tabs_ui"))
        )
      )
    )
  )
}

# =========================================================
# SERVER
# =========================================================

vennModuleServer <- function(id) {
  
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    # =====================================================
    # 1) Preserve Dynamic Input State
    # =====================================================
    
    stored_values <- reactiveValues()
    
    observe({
      
      n <- input$n_sets
      req(n)
      
      for (i in seq_len(n)) {
        
        name_id   <- paste0("name_", i)
        txt_id    <- paste0("txt_", i)
        col_id    <- paste0("col_", i)
        source_id <- paste0("source_", i)
        
        if (!is.null(input[[name_id]])) {
          stored_values[[name_id]] <- input[[name_id]]
        }
        
        if (!is.null(input[[txt_id]])) {
          stored_values[[txt_id]] <- input[[txt_id]]
        }
        
        if (!is.null(input[[col_id]])) {
          stored_values[[col_id]] <- input[[col_id]]
        }
        
        if (!is.null(input[[source_id]])) {
          stored_values[[source_id]] <- input[[source_id]]
        }
      }
    })
    
    # =====================================================
    # 2) Dynamic Data Input Blocks
    # =====================================================
    
    output$individual_sets_ui <- renderUI({
      
      n <- req(input$n_sets)
      default_cols <- default_group_colors(n)
      
      tagList(
        lapply(seq_len(n), function(i) {
          
          div(
            class = "set-block",
            
            div(
              class = "set-header",
              
              textInput(
                inputId = ns(paste0("name_", i)),
                label   = NULL,
                value   = isolate(stored_values[[paste0("name_", i)]]) %||% paste("Group", LETTERS[i]),
                width   = "65%"
              ),
              
              colourInput(
                inputId = ns(paste0("col_", i)),
                label   = NULL,
                value   = isolate(stored_values[[paste0("col_", i)]]) %||% default_cols[i],
                width   = "30%"
              )
            ),
            
            radioButtons(
              inputId  = ns(paste0("source_", i)),
              label    = NULL,
              choices  = c("Text" = "txt", "File" = "file"),
              selected = isolate(stored_values[[paste0("source_", i)]]) %||% "txt",
              inline   = TRUE
            ),
            
            conditionalPanel(
              condition = sprintf("input['%s'] == 'txt'", ns(paste0("source_", i))),
              
              textAreaInput(
                inputId     = ns(paste0("txt_", i)),
                label       = NULL,
                value       = isolate(stored_values[[paste0("txt_", i)]]) %||% "",
                placeholder = "Enter one ID per line...",
                rows        = 3
              )
            ),
            
            conditionalPanel(
              condition = sprintf("input['%s'] == 'file'", ns(paste0("source_", i))),
              
              fileInput(
                inputId = ns(paste0("file_", i)),
                label   = NULL,
                accept  = c(".txt", ".csv", ".tsv", ".xls", ".xlsx")
              )
            )
          )
        })
      )
    })
    
    # =====================================================
    # 3) Read Uploaded Files
    # =====================================================
    
    read_group_file <- function(file_path, original_name = NULL) {
      
      ext_source <- original_name %||% file_path
      ext <- tolower(tools::file_ext(ext_source))
      
      if (ext %in% c("xlsx", "xls")) {
        
        validate(
          need(
            requireNamespace("readxl", quietly = TRUE),
            "Package 'readxl' is required to read Excel files."
          )
        )
        
        x <- tryCatch(
          {
            readxl::read_excel(
              path = file_path,
              sheet = 1,
              col_names = FALSE,
              col_types = "text"
            )
          },
          error = function(e) NULL
        )
        
        validate(
          need(
            !is.null(x) && ncol(x) >= 1,
            "Invalid Excel file. The first sheet must contain at least one column."
          )
        )
        
        values <- x[[1]]
        
      } else if (ext == "csv") {
        
        x <- tryCatch(
          {
            readr::read_csv(
              file = file_path,
              col_names = FALSE,
              col_types = readr::cols(.default = readr::col_character()),
              show_col_types = FALSE
            )
          },
          error = function(e) NULL
        )
        
        validate(
          need(
            !is.null(x) && ncol(x) >= 1,
            "Invalid CSV file. The file must contain at least one column."
          )
        )
        
        values <- x[[1]]
        
      } else if (ext == "tsv") {
        
        x <- tryCatch(
          {
            readr::read_tsv(
              file = file_path,
              col_names = FALSE,
              col_types = readr::cols(.default = readr::col_character()),
              show_col_types = FALSE
            )
          },
          error = function(e) NULL
        )
        
        validate(
          need(
            !is.null(x) && ncol(x) >= 1,
            "Invalid TSV file. The file must contain at least one column."
          )
        )
        
        values <- x[[1]]
        
      } else {
        
        values <- tryCatch(
          {
            readLines(file_path, warn = FALSE)
          },
          error = function(e) character()
        )
      }
      
      values <- values %>%
        as.character() %>%
        stringr::str_trim()
      
      # Only conventional explicit headers are removed; arbitrary first IDs survive.
      headers <- c("gene", "genes", "gene_id", "gene_name", "gene_symbol", "symbol",
                   "id", "feature", "features", "feature_id", "protein", "protein_id")
      first <- which(!is.na(values) & values != "")[1]
      if (!is.na(first) && tolower(values[first]) %in% headers) values <- values[-first]
      unique(values[!is.na(values) & values != ""])
    }
    
    # =====================================================
    # 4) Group Configuration
    # =====================================================
    
    group_config <- reactive({
      
      n <- req(input$n_sets)
      
      raw_names <- sapply(seq_len(n), function(i) {
        
        nm <- input[[paste0("name_", i)]] %||% paste("Group", LETTERS[i])
        nm <- stringr::str_trim(nm)
        
        if (nm == "") {
          nm <- paste("Group", LETTERS[i])
        }
        
        nm
      })
      
      set_names <- make.unique(raw_names, sep = "_")
      
      cols <- sapply(seq_len(n), function(i) {
        input[[paste0("col_", i)]] %||% default_group_colors(n)[i]
      })
      
      tibble(
        Index   = seq_len(n),
        RawName = raw_names,
        SetName = set_names,
        Color   = cols
      )
    })
    
    # =====================================================
    # 5) Get Sets
    # =====================================================
    
    get_data <- reactive({
      
      cfg <- group_config()
      sets <- list()
      
      for (row_i in seq_len(nrow(cfg))) {
        
        i <- cfg$Index[row_i]
        nm <- cfg$SetName[row_i]
        
        src <- input[[paste0("source_", i)]] %||% "txt"
        
        values <- if (src == "txt") {
          
          txt <- input[[paste0("txt_", i)]] %||% ""
          
          txt %>%
            strsplit(split = "\\n") %>%
            unlist() %>%
            stringr::str_trim()
          
        } else {
          
          f <- input[[paste0("file_", i)]]
          
          if (!is.null(f)) {
            read_group_file(
              file_path = f$datapath,
              original_name = f$name
            )
          } else {
            character()
          }
        }
        
        values <- values %>%
          as.character() %>%
          stringr::str_trim()
        
        values <- unique(values[!is.na(values) & values != ""])
        
        if (length(values) > 0) {
          sets[[nm]] <- values
        }
      }
      
      validate(
        need(
          length(sets) >= 2,
          "Please provide data for at least 2 groups."
        )
      )
      
      sets
    })
    
    current_set_colors <- reactive({
      
      sets <- get_data()
      cfg <- group_config()
      
      cols <- cfg$Color
      names(cols) <- cfg$SetName
      
      cols[names(sets)]
    })
    
    # =====================================================
    # 6) Draw Venn / Euler Object
    # =====================================================
    
    draw_venn_obj <- function() {
      
      raw_sets <- get_data()
      raw_cols <- current_set_colors()

      names(raw_sets) <- gsub("&", "and", names(raw_sets))
      names(raw_cols) <- gsub("&", "and", names(raw_cols))
      
      if (input$plot_type == "Euler") {
        
        set_sigs <- sapply(raw_sets, function(x) paste(sort(x), collapse = "||"))
        unique_sigs <- unique(set_sigs)
        
        sets <- list()
        cols <- c()
        
        for (sig in unique_sigs) {
          idx <- which(set_sigs == sig)

          combined_name <- paste(names(raw_sets)[idx], collapse = "\n")
          
          sets[[combined_name]] <- raw_sets[[idx[1]]]

          mixed_color <- blend_colors(raw_cols[names(raw_sets)[idx]])
          names(mixed_color) <- combined_name 
          
          cols <- c(cols, mixed_color)
        }
        
        fit <- eulerr::euler(sets)
        
      } else {
  
        sets <- raw_sets
        cols <- raw_cols
        
        fit <- eulerr::venn(sets)
        
      }
      
      grid::grid.newpage()
      
      venn_grob <- plot(
        fit,
        fills = list(
          fill  = cols,
          alpha = input$opacity
        ),
        labels = list(
          font = 2,
          cex  = input$font_size / 12
        ),
        quantities = list(
          cex = (input$font_size / 12) * 0.8,
          col = "black"
        ),
        edges = list(
          col = "black",
          lwd = 1
        )
      )

      if (inherits(venn_grob, "grob") || inherits(venn_grob, "gTree")) {
        grid::grid.draw(venn_grob)
      }

      if (nzchar(input$plot_title %||% "")) {
        grid::grid.text(
          input$plot_title,
          x = 0.5,
          y = 0.97,
          gp = grid::gpar(fontsize = input$font_size + 2, fontface = "bold")
        )
      }

      invisible(venn_grob)
    }
    
    # =====================================================
    # 7) Intersection Membership Matrix
    # =====================================================
    
    intersection_matrix <- reactive({
      
      sets <- get_data()
      
      all_elements <- unique(unlist(sets))
      
      membership <- as.data.frame(
        lapply(sets, function(s) all_elements %in% s),
        check.names = FALSE
      )
      
      rownames(membership) <- all_elements
      
      membership
    })
    
    # =====================================================
    # 8) Long Intersection Table
    # =====================================================
    
    intersection_table <- reactive({
      
      df <- intersection_matrix()
      set_names <- colnames(df)
      
      df_out <- df %>%
        rownames_to_column("ID") %>%
        mutate(
          Intersection = apply(
            df,
            1,
            function(x) {
              active_sets <- set_names[x]
              format_intersection_label(
                active_sets = active_sets,
                all_sets = set_names
              )
            }
          )
        ) %>%
        relocate(ID, Intersection)
      
      df_out
    })
    
    # =====================================================
    # 9) Intersection Summary
    # =====================================================
    
    intersection_summary <- reactive({
      
      df <- intersection_table()
      set_names <- colnames(intersection_matrix())
      
      out <- df %>%
        group_by(Intersection) %>%
        summarise(
          Count = n(),
          across(
            all_of(set_names),
            ~ any(.x),
            .names = "{.col}"
          ),
          .groups = "drop"
        )
      
      out$Degree <- rowSums(
        out[, set_names, drop = FALSE]
      )
      
      out %>%
        arrange(desc(Count), desc(Degree), Intersection) %>%
        relocate(Intersection, Count, Degree)
    })
    
    intersection_profiles <- reactive({
      
      intersection_summary() %>%
        mutate(
          TabID = paste0("intersection_table_", row_number())
        )
    })
    
    set_size_summary <- reactive({
      
      sets <- get_data()
      
      tibble(
        Group = names(sets),
        Size = lengths(sets)
      ) %>%
        arrange(desc(Size))
    })
    
    # =====================================================
    # 10) UpSet Plot
    # =====================================================
    
    upset_intersections <- reactive({
      
      df <- intersection_summary()
      
      df <- df %>%
        filter(Count >= input$min_intersection_size)
      
      if (input$upset_order == "count") {
        df <- df %>%
          arrange(desc(Count), desc(Degree), Intersection)
      } else if (input$upset_order == "degree") {
        df <- df %>%
          arrange(desc(Degree), desc(Count), Intersection)
      } else {
        df <- df %>%
          arrange(Intersection)
      }
      
      df %>%
        slice_head(n = input$top_intersections) %>%
        mutate(
          IntersectionPlot = stringr::str_wrap(
            Intersection,
            width = 18
          ),
          IntersectionPlot = factor(
            IntersectionPlot,
            levels = unique(IntersectionPlot)
          )
        )
    })
    
    create_upset_plot <- reactive({
      
      df <- upset_intersections()
      set_names <- colnames(intersection_matrix())
      
      upset_bar_color <- input$upset_bar_color %||% "#2A78D6"
      upset_active_dot_color <- input$upset_active_dot_color %||% "#222222"
      upset_inactive_dot_color <- input$upset_inactive_dot_color %||% "#d9d9d9"
      upset_line_color <- input$upset_line_color %||% "#333333"
      
      upset_font_factor <- (input$font_size %||% 12) / 12
      upset_base_size <- input$font_size %||% 12
      dot_size_inactive <- 2.64 * upset_font_factor
      dot_size_active <- 3.36 * upset_font_factor
      line_width <- 0.42 * upset_font_factor
      
      validate(
        need(
          nrow(df) > 0,
          "No intersections pass the current UpSet filters."
        )
      )
      
      matrix_df <- df %>%
        select(Intersection, IntersectionPlot, all_of(set_names)) %>%
        pivot_longer(
          cols = all_of(set_names),
          names_to = "Set",
          values_to = "Active"
        ) %>%
        mutate(
          Set = factor(Set, levels = rev(set_names))
        )
      
      segment_df <- matrix_df %>%
        filter(Active) %>%
        group_by(IntersectionPlot) %>%
        summarise(
          y_min = Set[which.min(as.numeric(Set))],
          y_max = Set[which.max(as.numeric(Set))],
          .groups = "drop"
        )
      
      top_bar <- ggplot(
        df,
        aes(
          x = IntersectionPlot,
          y = Count
        )
      ) +
        geom_col(
          fill = upset_bar_color,
          width = 0.72
        ) +
        theme_bw(base_size = upset_base_size) +
        labs(
          x = NULL,
          y = "Intersection size"
        ) +
        theme(
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = upset_base_size, color = "black"),
          axis.title.y = element_text(size = upset_base_size, color = "black"),
          plot.margin = margin(5, 10, 0, 10)
        )
      
      matrix_plot <- ggplot(
        matrix_df,
        aes(
          x = IntersectionPlot,
          y = Set
        )
      ) +
        geom_point(
          data = matrix_df %>% filter(!Active),
          color = upset_inactive_dot_color,
          size = dot_size_inactive
        ) +
        geom_segment(
          data = segment_df,
          aes(
            x = IntersectionPlot,
            xend = IntersectionPlot,
            y = y_min,
            yend = y_max
          ),
          inherit.aes = FALSE,
          linewidth = line_width,
          color = upset_line_color
        ) +
        geom_point(
          data = matrix_df %>% filter(Active),
          color = upset_active_dot_color,
          size = dot_size_active
        ) +
        theme_bw(base_size = upset_base_size) +
        labs(
          x = "Intersection",
          y = NULL
        ) +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(
            angle = 35,
            hjust = 1,
            size = upset_base_size,
            color = "black"
          ),
          axis.text.y = element_text(size = upset_base_size, color = "black"),
          axis.title.x = element_text(size = upset_base_size, color = "black"),
          plot.margin = margin(0, 10, 5, 10)
        )

      upset_grid <- cowplot::plot_grid(
        top_bar,
        matrix_plot,
        ncol = 1,
        align = "v",
        axis = "lr",
        rel_heights = c(0.55, 0.45)
      )

      if (nzchar(input$plot_title %||% "")) {
        cowplot::plot_grid(
          cowplot::ggdraw() +
            cowplot::draw_label(
              input$plot_title,
              fontface = "bold",
              size = upset_base_size + 2
            ),
          upset_grid,
          ncol = 1,
          rel_heights = c(0.08, 0.92)
        )
      } else {
        upset_grid
      }
    })
    
    # =====================================================
    # 11) Set Size Barplot
    # =====================================================
    
    create_set_size_plot <- reactive({
      
      df <- set_size_summary()
      cols <- current_set_colors()
      
      df <- df %>%
        mutate(
          Group = factor(
            Group,
            levels = rev(Group)
          )
        )
      
      ggplot(
        df,
        aes(
          x = Size,
          y = Group,
          fill = Group
        )
      ) +
        geom_col(
          width = 0.70,
          color = "grey25",
          linewidth = 0.35,
          alpha = input$opacity
        ) +
        scale_fill_manual(
          values = cols
        ) +
        theme_bw(base_size = input$font_size %||% 12) +
        labs(
          title = input$size_plot_title,
          x = "Number of unique IDs",
          y = NULL
        ) +
        theme(
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.text = element_text(size = input$font_size, color = "black"),
          axis.title = element_text(size = input$font_size, color = "black"),
          plot.title = element_text(size = input$font_size + 2, face = "bold", color = "black", hjust = 0.5),
          legend.position = "none"
        )
    })
    
    # =====================================================
    # 12) Render Main Plot
    # =====================================================
    
    output$main_plot <- renderPlot(
      
      expr = {
        
        if (input$plot_type == "UpSet") {
          print(create_upset_plot())
        } else {
          draw_venn_obj()
        }
      },
      
      height = function() {
        if (input$plot_type == "UpSet") {
          650
        } else {
          600
        }
      },
      
      width = function() {
        if (input$plot_type == "UpSet") {
          850
        } else {
          600
        }
      },
      
      res = 120
    )
    
    output$set_size_plot <- renderPlot({
      
      req(create_set_size_plot())
      
      print(create_set_size_plot())
      
    }, res = 120)
    
    # =====================================================
    # 13) Intersection Summary Table
    # =====================================================
    
    output$intersection_summary_table <- renderDT({
      
      datatable(
        intersection_summary(),
        rownames = FALSE,
        class = "cell-border stripe hover compact",
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          autoWidth = FALSE,
          columnDefs = list(
            list(
              targets = 0,
              width = "360px",
              className = "dt-intersection-cell"
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
    # 14) Intersections Explorer UI
    # =====================================================
    
    output$tabs_ui <- renderUI({
      
      profiles <- intersection_profiles()
      
      tab_panels <- lapply(seq_len(nrow(profiles)), function(i) {
        
        profile_name <- profiles$Intersection[i]
        profile_n    <- profiles$Count[i]
        output_id    <- profiles$TabID[i]
        
        tabPanel(
          title = paste0(profile_name, " (", profile_n, ")"),
          br(),
          DTOutput(ns(output_id))
        )
      })
      
      do.call(tabsetPanel, tab_panels)
    })
    
    # =====================================================
    # 15) Render Dynamic Intersection Tables
    # =====================================================
    
    observe({
      
      df <- intersection_table()
      profiles <- intersection_profiles()
      
      for (i in seq_len(nrow(profiles))) {
        
        local({
          
          profile_local <- profiles$Intersection[i]
          output_id <- profiles$TabID[i]
          
          output[[output_id]] <- renderDT({
            
            df %>%
              filter(Intersection == profile_local) %>%
              select(ID, Intersection, everything()) %>%
              datatable(
                rownames = FALSE,
                class = "cell-border stripe hover compact",
                options = list(
                  pageLength = 10,
                  scrollX = TRUE,
                  autoWidth = FALSE,
                  columnDefs = list(
                    list(
                      targets = 1,
                      width = "360px",
                      className = "dt-intersection-cell"
                    ),
                    list(
                      targets = "_all",
                      className = "dt-center-cell"
                    )
                  )
                )
              )
          })
        })
      }
    })
    
    # =====================================================
    # 16) Download Main Plot
    # =====================================================
    
    output$download_plot <- downloadHandler(
      
      filename = function() {
        plot_label <- stringr::str_to_lower(input$plot_type)
        paste0(plot_label, "_plot_", Sys.Date(), "_", input$exp_width, "x", input$exp_height, ".", input$export_format)
      },
      
      content = function(file) {
        
        width_in  <- input$exp_width * 0.393701
        height_in <- input$exp_height * 0.393701
        
        open_venn_device(
          file   = file,
          fmt    = input$export_format,
          width  = width_in,
          height = height_in,
          dpi    = input$export_dpi
        )
        
        if (input$plot_type == "UpSet") {
          print(create_upset_plot())
        } else {
          draw_venn_obj()
        }
        
        dev.off()
      }
    )
    
    # =====================================================
    # 16.1) Download Set Size Plot
    # =====================================================
    
    output$download_size_plot <- downloadHandler(
      
      filename = function() {
        paste0("set_size_plot_", Sys.Date(), "_", input$exp_width_size, "x", input$exp_height_size, ".", input$export_format)
      },

      content = function(file) {

        width_in  <- input$exp_width_size * 0.393701
        height_in <- input$exp_height_size * 0.393701

        open_venn_device(
          file   = file,
          fmt    = input$export_format,
          width  = width_in,
          height = height_in,
          dpi    = input$export_dpi
        )

        print(create_set_size_plot())
        
        dev.off()
      }
    )
    
    # =====================================================
    # 17) Download All Intersections
    # =====================================================
    
    output$download_data <- downloadHandler(
      
      filename = function() {
        paste0("venn_intersections_", Sys.Date(), ".txt")
      },
      
      content = function(file) {
        
        readr::write_delim(
          x     = intersection_table(),
          file  = file,
          delim = "\t"
        )
      }
    )
  })
}