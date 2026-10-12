# =========================================================
# MODULE: Shared Helpers (EN)
# =========================================================

# =========================================================
# Helper Functions
# =========================================================

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)) || identical(x, "")) {
    y
  } else {
    x
  }
}

# Strict per-value conversion: decimal comma/dot and scientific notation.
# Do not guess grouping separators or extract numbers from malformed text.
omics_numeric <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  text <- trimws(as.character(x))
  missing <- is.na(text) | toupper(text) %in% c("", "NA", "N/A", "NULL", "NAN")
  special <- !missing & text %in% c("Inf", "+Inf", "-Inf", "Infinity", "+Infinity", "-Infinity")
  valid <- missing | special | grepl(
    "^[+-]?([0-9]+([.,][0-9]*)?|[.,][0-9]+)([eE][+-]?[0-9]+)?$", text
  )
  shiny::validate(shiny::need(all(valid), paste0(
    "Invalid numeric value(s): ", paste(head(unique(text[!valid]), 5), collapse = ", "),
    ". Use decimal comma or decimal dot, without thousands separators or extra text."
  )))
  text[missing] <- NA_character_
  text <- sub("Infinity", "Inf", text, fixed = TRUE)
  suppressWarnings(as.numeric(chartr(",", ".", text)))
}

omics_validate_ids <- function(ids, label = "Sample IDs") {
  ids <- trimws(as.character(ids))
  shiny::validate(
    shiny::need(all(!is.na(ids) & nzchar(ids)), paste0(label, " must not be empty.")),
    shiny::need(!anyDuplicated(ids), paste0(label, " must be unique. Duplicates: ",
      paste(head(unique(ids[duplicated(ids)]), 5), collapse = ", ")))
  )
  ids
}

omics_validate_table <- function(df) {
  force(df)
  colnames(df) <- omics_validate_ids(colnames(df), "Column names")
  df
}

omics_delimited_frame <- function(x) {
  problems <- readr::problems(x)
  shiny::validate(shiny::need(nrow(problems) == 0,
    "Inconsistent table rows. Check the delimiter and the number of fields in each row."))
  as.data.frame(x)
}

omics_excluded_rows <- function(n, context) {
  if (n > 0 && !is.null(shiny::getDefaultReactiveDomain())) {
    shiny::showNotification(paste0(context, ": ", n,
      " row(s) excluded because IDs or numeric values are missing."),
      type = "warning", duration = 10)
  }
}

omics_clean_id_list <- function(values) {
  values <- trimws(as.character(values))
  headers <- c("gene", "genes", "gene_id", "gene_name", "gene_symbol", "symbol",
               "id", "feature", "features", "feature_id", "protein", "protein_id")
  first <- which(!is.na(values) & nzchar(values))[1]
  if (!is.na(first) && tolower(values[first]) %in% headers) values <- values[-first]
  unique(values[!is.na(values) & nzchar(values)])
}

omics_read_id_list <- function(path, original_name = NULL) {
  ext <- tolower(tools::file_ext(original_name %||% path))
  if (ext %in% c("xlsx", "xls")) {
    x <- readxl::read_excel(path, sheet = 1, col_names = FALSE, col_types = "text")
    values <- x[[1]]
  } else if (ext %in% c("csv", "tsv")) {
    x <- readr::read_delim(path, delim = if (ext == "csv") "," else "\t",
      col_names = FALSE, col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE)
    x <- omics_delimited_frame(x)
    values <- x[[1]]
  } else {
    values <- sub("[,;\t].*$", "", readLines(path, warn = FALSE))
  }
  omics_clean_id_list(values)
}

omics_numeric_frame <- function(df) {
  df <- omics_validate_table(df)
  df[] <- lapply(df, omics_numeric)
  empty <- !vapply(df, function(x) any(is.finite(x)), logical(1))
  shiny::validate(shiny::need(!any(empty), paste0(
    "Columns without finite numeric values: ", paste(names(df)[empty], collapse = ", "),
    ". Correct these columns; samples cannot be silently discarded."
  )))
  df
}

omics_annotation <- function(anno, key, samples, required = character()) {
  anno <- omics_validate_table(anno)
  shiny::validate(shiny::need(key %in% names(anno), "Sample ID column not found in annotation."))
  ids <- trimws(as.character(anno[[key]]))
  shiny::validate(shiny::need(all(!is.na(ids) & nzchar(ids)), "Annotation contains empty sample IDs."))
  anno[[key]] <- ids
  duplicates <- unique(ids[duplicated(ids)])
  conflicts <- duplicates[vapply(duplicates, function(id) {
    nrow(unique(anno[ids == id, , drop = FALSE])) > 1
  }, logical(1))]
  shiny::validate(shiny::need(length(conflicts) == 0, paste0(
    "Conflicting annotation rows for sample(s): ", paste(head(conflicts, 10), collapse = ", "), "."
  )))
  anno <- anno[!duplicated(ids), , drop = FALSE]
  missing <- setdiff(samples, anno[[key]])
  shiny::validate(shiny::need(length(missing) == 0, paste0(
    "Missing annotation for ", length(missing), " sample(s): ",
    paste(head(missing, 10), collapse = ", "), ". Check sample IDs."
  )))
  required <- setdiff(required, c("__none__", ""))
  for (column in required) {
    shiny::validate(shiny::need(column %in% names(anno), paste0("Annotation column not found: ", column)))
    values <- trimws(as.character(anno[[column]][match(samples, anno[[key]])]))
    shiny::validate(shiny::need(all(!is.na(values) & nzchar(values)),
      paste0("Missing annotation values in column: ", column)))
  }
  anno
}

omics_number_parameter <- function(x, label, min = -Inf, max = Inf, integer = FALSE) {
  valid <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x >= min && x <= max && (!integer || x == floor(x))
  shiny::validate(shiny::need(valid, paste0(label, " must be ",
    if (integer) "an integer" else "a finite number", " between ", min, " and ", max, ".")))
  invisible(x)
}

omics_set_sizes <- function(min_size, max_size) {
  omics_number_parameter(min_size, "Minimum set size", 1, .Machine$integer.max, TRUE)
  omics_number_parameter(max_size, "Maximum set size", min_size, .Machine$integer.max, TRUE)
}

omics_probability <- function(x, label = "P-values") {
  valid <- is.na(x) | (is.finite(x) & x >= 0 & x <= 1)
  shiny::validate(shiny::need(all(valid), paste0(label, " must be finite values between 0 and 1.")))
  invisible(x)
}

# Zero may represent numerical underflow. Use a finite floor only for drawing;
# statistical filtering and exported values continue to use the original zero.
omics_plot_probabilities <- function(x) {
  positive <- x[is.finite(x) & x > 0]
  floor <- .Machine$double.xmin * .Machine$double.eps
  replacement <- if (length(positive)) max(min(positive) / 10, floor) else floor
  x[!is.na(x) & x == 0] <- replacement
  x
}

# A result is usable only with the data/analysis parameters that produced it.
# Appearance and display filters remain reactive without rerunning the analysis.
omics_event_reactive <- function(eventExpr, valueExpr, state, label,
                                 env = parent.frame()) {
  event <- substitute(eventExpr)
  value <- substitute(valueExpr)
  force(state)
  force(label)
  cached <- shiny::eventReactive(event, {
    stamp <- state()
    answer <- eval(value, envir = new.env(parent = env))
    list(state = stamp, result = answer)
  }, event.quoted = TRUE, event.env = env)
  shiny::reactive({
    out <- cached()
    shiny::validate(shiny::need(identical(out$state, state()), paste0(
      "Data or analysis parameters changed. Run ", label, " again to update results."
    )))
    out$result
  })
}

# =========================================================
# MSigDB Category / Subcategory Choices
# =========================================================

.omics_msigdb_catalog_cache <- new.env(parent = emptyenv())

omics_msigdb_catalog <- function() {
  if (!exists("HS", envir = .omics_msigdb_catalog_cache, inherits = FALSE)) {
    shiny::validate(shiny::need(
      "db_species" %in% names(formals(msigdbr::msigdbr_collections)),
      "Update msigdbr to a version supporting the MSigDB collection catalog."
    ))
    catalog <- tryCatch(msigdbr::msigdbr_collections(db_species = "HS"), error = identity)
    shiny::validate(shiny::need(!inherits(catalog, "error"), paste0(
      "Cannot load the MSigDB collection catalog. Check the network or msigdbr cache. ",
      if (inherits(catalog, "error")) conditionMessage(catalog) else ""
    )))
    shiny::validate(shiny::need(nrow(catalog) > 0, "No collections available for this MSigDB database."))
    .omics_msigdb_catalog_cache[["HS"]] <- catalog
  }
  .omics_msigdb_catalog_cache[["HS"]]
}

omics_msigdb_species_choices <- function() {
  sort(unique(msigdbr::msigdbr_species()$species_name))
}

msigdb_category_choices <- function() {
  catalog <- omics_msigdb_catalog()
  codes <- unique(catalog$gs_collection)
  labels <- catalog$gs_collection_name[match(codes, catalog$gs_collection)]
  # Some catalog rows label the subcollection rather than the main collection.
  main_labels <- c(H = "Hallmark gene sets", C1 = "Positional gene sets",
    C2 = "Curated gene sets", C3 = "Regulatory gene sets", C4 = "Computational gene sets",
    C5 = "Ontology gene sets", C6 = "Oncogenic gene sets", C7 = "Immunologic gene sets",
    C8 = "Cell type signature gene sets", C9 = "Computational perturbation signatures")
  known <- codes %in% names(main_labels)
  labels[known] <- unname(main_labels[codes[known]])
  stats::setNames(codes, paste0(codes, ": ", labels))
}

msigdb_subcategory_choices <- function(category) {
  catalog <- omics_msigdb_catalog()
  codes <- unique(catalog$gs_subcollection[catalog$gs_collection == category])
  codes <- codes[!is.na(codes) & nzchar(codes)]
  c("All subcollections" = "", stats::setNames(codes, codes))
}

omics_msigdb_ui <- function(ns) {
  shiny::tagList(
    shiny::selectInput(ns("species"), "Species:",
      choices = omics_msigdb_species_choices(), selected = "Homo sapiens"),
    shiny::selectInput(ns("category"), "Main category:", choices = NULL),
    shiny::selectInput(ns("subcategory"), "Subcategory:", choices = NULL),
    shiny::uiOutput(ns("msigdb_note"))
  )
}

# A single species menu uses the HS reference collections for all supported
# organisms, as in the original app. Uploaded identifiers are never converted.
omics_msigdb_server <- function(input, session, source_input) {
  source <- shiny::reactive(input[[source_input]] %||% "msigdb")
  shiny::observeEvent(source(), {
    shiny::req(source() == "msigdb")
    choices <- omics_msigdb_species_choices()
    selected <- shiny::isolate(input$species)
    if (is.null(selected) || !selected %in% choices) {
      selected <- "Homo sapiens"
    }
    shiny::updateSelectInput(session, "species", choices = choices, selected = selected)
    categories <- msigdb_category_choices()
    selected <- shiny::isolate(input$category)
    if (is.null(selected) || !selected %in% categories) {
      selected <- if ("H" %in% categories) "H" else unname(categories[1])
    }
    shiny::updateSelectInput(session, "category", choices = categories, selected = selected)
  }, ignoreInit = FALSE)
  shiny::observeEvent(list(source(), input$category), {
    shiny::req(source() == "msigdb", input$category)
    shiny::req(input$category %in% msigdb_category_choices())
    choices <- msigdb_subcategory_choices(input$category)
    selected <- shiny::isolate(input$subcategory) %||% ""
    if (!selected %in% choices) selected <- ""
    shiny::updateSelectInput(session, "subcategory", choices = choices, selected = selected)
  }, ignoreInit = FALSE)
  session$output$msigdb_note <- shiny::renderUI({
    shiny::req(source() == "msigdb")
    # Loading the catalog here also gives a visible validation error if unavailable.
    omics_msigdb_catalog()
    shiny::p(class = "help-block",
      "Select the species of your data. Uploaded IDs are preserved; only exact matches are used for gene set membership. Unmatched IDs do not block the remaining matches.")
  })
  shiny::reactive({
    if (source() != "msigdb") return(NULL)
    subcategory <- input$subcategory %||% ""
    if (nzchar(subcategory)) subcategory else NULL
  })
}

omics_get_msigdb <- function(species, category, subcategory = NULL) {
  shiny::req(species, category)
  shiny::validate(
    shiny::need(species %in% omics_msigdb_species_choices(),
      "Species is not supported by msigdbr."),
    shiny::need(category %in% msigdb_category_choices(),
      "Collection is not available in the selected MSigDB database."),
    shiny::need(is.null(subcategory) || subcategory %in% msigdb_subcategory_choices(category),
      "Subcollection is not available in the selected MSigDB collection.")
  )
  args <- list(species = species, db_species = "HS", collection = category)
  if (!is.null(subcategory) && nzchar(subcategory)) args$subcollection <- subcategory
  msig <- tryCatch(do.call(msigdbr::msigdbr, args), error = identity)
  shiny::validate(shiny::need(!inherits(msig, "error"), paste0(
    "Cannot load the selected MSigDB gene sets. ",
    if (inherits(msig, "error")) conditionMessage(msig) else ""
  )))
  shiny::validate(shiny::need(nrow(msig) > 0 && all(c("gs_name", "gene_symbol") %in% names(msig)),
    "The selected MSigDB collection has no usable gene symbols."))
  msig
}

# GMT: name, description (may be blank), then one or more exact identifiers.
# Duplicate set names are ambiguous and must never be silently dropped/merged.
omics_read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- lines[nzchar(trimws(lines))]
  shiny::validate(shiny::need(length(lines) > 0, "The GMT file is empty."))
  # Preserve an original empty trailing field so malformed rows are rejected.
  fields <- strsplit(paste0(lines, "\t"), "\t", fixed = TRUE)
  valid <- vapply(fields, function(x) length(x) >= 3L && nzchar(x[1]) &&
    all(nzchar(x[-c(1, 2)])) && all(x[-2] == trimws(x[-2])), logical(1))
  shiny::validate(shiny::need(all(valid), paste0(
    "Invalid GMT row(s): ", paste(head(which(!valid), 5), collapse = ", "),
    ". Use tab-separated set name, description, and exact IDs without surrounding whitespace."
  )))
  set_names <- vapply(fields, `[[`, character(1), 1)
  shiny::validate(shiny::need(!anyDuplicated(set_names),
    "GMT gene set names must be unique; duplicate names are not merged."))
  stats::setNames(lapply(fields, function(x) unique(x[-c(1, 2)])), set_names)
}

omics_validate_gene_set_ids <- function(ids, gene_sets, label = "Input") {
  ids <- unique(as.character(ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  reference <- unique(as.character(unlist(gene_sets, use.names = FALSE)))
  reference <- reference[!is.na(reference) & nzchar(reference)]
  unmatched <- setdiff(ids, reference)
  case_only <- unmatched[toupper(unmatched) %in% toupper(reference)]
  n_exact <- length(intersect(ids, reference))
  shiny::validate(shiny::need(n_exact > 0, paste0(label,
    " has no exact identifier overlap with the selected gene sets. Check species and identifier format. ",
    if (length(case_only)) paste0(length(case_only), " ID(s) match only if capitalization is ignored. ") else "",
    "No gene/protein IDs are converted.")))
  session <- shiny::getDefaultReactiveDomain()
  if (length(case_only) && !is.null(session)) {
    shiny::showNotification(paste0(label, ": ", length(case_only),
      " ID(s) differ in capitalization and will not match the gene sets: ",
      paste(head(case_only, 5), collapse = ", "), ". Continuing with ", n_exact,
      " exact matches; uploaded IDs are unchanged."), type = "warning", duration = 10,
      id = session$ns(paste0("omics-case-", gsub("[^[:alnum:]]", "-", label))))
  }
  invisible(n_exact)
}
