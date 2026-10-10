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

msigdb_category_choices <- function() {
  c(
    "H: hallmark gene sets" = "H",
    "C1: positional gene sets" = "C1",
    "C2: curated gene sets" = "C2",
    "C3: regulatory gene sets" = "C3",
    "C4: computational gene sets" = "C4",
    "C5: ontology gene sets" = "C5",
    "C6: oncogenic gene sets" = "C6",
    "C7: immunologic gene sets" = "C7",
    "C8: cell type signature gene sets" = "C8"
  )
}

msigdb_subcategory_choices <- function(category) {

  switch(
    category,

    "H" = c("None" = ""),
    "C1" = c("None" = ""),

    "C2" = c(
      "CGP: chemical and genetic perturbations" = "CGP",
      "CP:BIOCARTA: BioCarta gene sets" = "CP:BIOCARTA",
      "CP:KEGG_MEDICUS: KEGG Medicus gene sets" = "CP:KEGG_MEDICUS",
      "CP:PID: PID gene sets" = "CP:PID",
      "CP:REACTOME: Reactome gene sets" = "CP:REACTOME",
      "CP:WIKIPATHWAYS: WikiPathways gene sets" = "CP:WIKIPATHWAYS",
      "CP:KEGG_LEGACY: KEGG Legacy gene sets" = "CP:KEGG_LEGACY"
    ),

    "C3" = c(
      "MIR:MIRDB: miRDB microRNA targets" = "MIR:MIRDB",
      "MIR:MIR_LEGACY: legacy microRNA targets" = "MIR:MIR_LEGACY",
      "TFT:GTRD: GTRD transcription factor targets" = "TFT:GTRD",
      "TFT:TFT_LEGACY: legacy transcription factor targets" = "TFT:TFT_LEGACY"
    ),

    "C4" = c(
      "3CA: Curated Cancer Cell Atlas gene sets" = "3CA",
      "CGN: cancer gene neighborhoods" = "CGN",
      "CM: cancer modules" = "CM"
    ),

    "C5" = c(
      "GO:BP: GO biological process" = "GO:BP",
      "GO:CC: GO cellular component" = "GO:CC",
      "GO:MF: GO molecular function" = "GO:MF",
      "HPO: Human Phenotype Ontology" = "HPO"
    ),

    "C6" = c("None" = ""),

    "C7" = c(
      "IMMUNESIGDB: ImmuneSigDB gene sets" = "IMMUNESIGDB",
      "VAX: vaccine response gene sets" = "VAX"
    ),

    "C8" = c("None" = "")
  )
}
