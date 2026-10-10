# OmicPrism home page and bundled example downloads.

example_files <- list(
  list(id = "example_01", file = "01_Overlap_ListA.xlsx", label = "List A \u00b7 XLSX"),
  list(id = "example_02", file = "02_Overlap_ListB.txt", label = "List B \u00b7 TXT"),
  list(id = "example_03", file = "03_Overlap_ListC.csv", label = "List C \u00b7 CSV"),
  list(id = "example_04", file = "04_PCA_PXD048675_Neutrophil_Abundance.xlsx", label = "Expression matrix"),
  list(id = "example_05", file = "05_PCA_PXD048675_Neutrophil_SampleAnnotation.txt", label = "Sample annotation"),
  list(id = "example_06", file = "06_Volcano_PXD048675_Neutrophil_Limma_CAPvsControl.txt", label = "Differential results"),
  list(id = "example_07", file = "07_Heatmap_Leite2021_Sepsis_Abundance.txt", label = "Expression matrix"),
  list(id = "example_08", file = "08_Heatmap_Leite2021_Sepsis_SampleAnnotation.txt", label = "Sample annotation"),
  list(id = "example_09", file = "09_Mfuzz_Hato2019_Sepsis_Abundance.txt", label = "Expression matrix"),
  list(id = "example_10", file = "10_ORA_Hato2019_Sepsis_MfuzzClusterCoreProteins.txt", label = "Cluster feature list"),
  list(id = "example_11", file = "11_ORA_Hato2019_Sepsis_Background.txt", label = "Background universe"),
  list(id = "example_12", file = "12_ssGSEA_PXD057739_Bcell_Abundance.txt", label = "Expression matrix"),
  list(id = "example_13", file = "13_ssGSEA_PXD057739_Bcell_SampleAnnotation.txt", label = "Sample annotation"),
  list(id = "example_14", file = "14_SingleGSEA_Leite2021_Sepsis_Limma.txt", label = "Ranked results table"),
  list(id = "example_15", file = "15_ComparativeGSEA_PXD075181_RAW264.7_Limma_LPS15min_vs_Control.txt", label = "LPS 15 min contrast"),
  list(id = "example_16", file = "16_ComparativeGSEA_PXD075181_RAW264.7_Limma_LPS6h_vs_Control.txt", label = "LPS 6 h contrast")
)

home_modules <- list(
  list(id = "go_overlap", icon = "circle-half-stroke", title = "Overlap", category = "SET ANALYSIS",
       text = "Compare gene or protein lists with Venn, Euler and UpSet diagrams.",
       input = "Gene / protein lists", examples = 1:3),
  list(id = "go_pca", icon = "chart-simple", title = "PCA", category = "SAMPLE QC",
       text = "Inspect sample structure, group separation and correlations.",
       input = "Expression matrix + sample annotation", examples = 4:5),
  list(id = "go_volcano", icon = "mountain", title = "Volcano plot", category = "DIFFERENTIAL",
       text = "Explore differential expression using effect size and significance.",
       input = "Differential results table", examples = 6L),
  list(id = "go_heatmap", icon = "table-cells", title = "Heatmap", category = "EXPRESSION",
       text = "View expression patterns across samples with clustering and annotation.",
       input = "Expression matrix + optional annotation", examples = 7:8),
  list(id = "go_mfuzz", icon = "chart-line", title = "Time-series Mfuzz", category = "TIME SERIES",
       text = "Find soft clusters of features with similar temporal profiles.",
       input = "Ordered expression matrix", examples = 9L),
  list(id = "go_ora", icon = "list-check", title = "ORA", category = "ENRICHMENT",
       text = "Identify enriched pathways in a feature list against a background universe.",
       input = "Feature list + background", examples = 10:11),
  list(id = "go_ssgsea", icon = "chart-column", title = "ssGSEA", category = "PATHWAY SCORES",
       text = "Score pathway activity in individual samples and compare groups.",
       input = "Expression matrix + sample annotation", examples = 12:13),
  list(id = "go_gsea_single", icon = "chart-line", title = "Single GSEA", category = "ENRICHMENT",
       text = "Test pathway enrichment along a ranked list from one contrast.",
       input = "Ranked results table", examples = 14L),
  list(id = "go_gsea_comp", icon = "layer-group", title = "Comparative GSEA", category = "CONTRASTS",
       text = "Compare pathway enrichment between two ranked contrasts.",
       input = "Two ranked results tables", examples = 15:16)
)

tool_card <- function(module, index) {
  tags$article(
    class = "tool-card",
    tags$div(class = "tool-heading",
      tags$span(class = "tool-icon", icon(module$icon)),
      tags$h3(class = "tool-title", module$title),
      tags$span(class = "tool-category", title = sprintf("Module %02d", index), module$category)
    ),
    tags$p(class = "tool-text", module$text),
    tags$div(class = "tool-meta",
      tags$span(class = "input-label", "Input"),
      tags$span(module$input)
    ),
    tags$div(class = "tool-actions",
      actionLink(module$id, tagList("Open module", icon("arrow-right")),
                 class = "module-open")
    ),
    tags$div(class = "tool-examples",
      tags$p(class = "examples-label", "Download examples"),
      tags$div(class = "example-links",
        lapply(module$examples, function(i) {
          example <- example_files[[i]]
          downloadLink(example$id, tagList(icon("download"), example$label),
                       class = "example-link", title = example$file)
        })
      )
    )
  )
}

home_page_ui <- function() {
  tabPanel(
    title = "Home", value = "home", icon = icon("house"),
    tags$main(class = "home-workspace",
      tags$section(class = "home-intro",
        tags$div(class = "intro-copy",
          tags$p(class = "home-kicker", "Bioinformatics workbench"),
          tags$h1("OmicPrism"),
          tags$p(class = "home-description",
                 "Expression analysis, temporal clustering and pathway enrichment."),
          tags$p(class = "home-context", "Proteomics \u00b7 Transcriptomics")
        ),
        tags$img(src = "omicprism-icon.png", class = "home-logo",
                 alt = "OmicPrism prism emblem", width = 128, height = 148)
      ),
      tags$section(class = "module-directory",
        tags$div(class = "directory-heading",
          tags$div(tags$h2("Analysis modules"),
                    tags$p("Choose a module or download its example files to try the workflow.")),
          tags$p(class = "directory-count", "9 modules \u00b7 16 example files")
        ),
        tags$div(class = "module-grid", lapply(seq_along(home_modules), function(i) {
          tool_card(home_modules[[i]], i)
        }))
      )
    )
  )
}

register_example_downloads <- function(output) {
  for (example in example_files) {
    local({
      entry <- example
      source_path <- file.path("case", entry$file)
      output[[entry$id]] <- downloadHandler(
        filename = function() entry$file,
        contentType = switch(tolower(tools::file_ext(entry$file)),
          xlsx = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
          csv = "text/csv", "text/plain"),
        content = function(file) {
          if (!file.exists(source_path)) stop("The example file is currently unavailable.", call. = FALSE)
          if (!file.copy(source_path, file, overwrite = TRUE)) {
            stop("The example file could not be downloaded. Please try again.", call. = FALSE)
          }
        }
      )
    })
  }
}
