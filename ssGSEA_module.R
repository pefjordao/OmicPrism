# =========================================================
# EzOmicsViz — ssGSEA Module
# =========================================================

library(shiny)
library(shinyjs)
library(tidyverse)
library(DT)
library(GSVA)
library(msigdbr)
library(svglite)
library(colourpicker)
library(rstatix)
library(ggpubr)
library(readxl)

# =========================================================
# Helper Functions
# =========================================================

ssgsea_clean_column_names <- function(df) {
  
  colnames(df) <- colnames(df) %>%
    as.character() %>%
    stringr::str_trim()
  
  df
}

ssgsea_try_read_delim <- function(path, delim, decimal_mark = ".") {
  
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
        ssgsea_clean_column_names()
    },
    error = function(e) NULL
  )
}

ssgsea_try_read_whitespace <- function(path, decimal_mark = ".") {
  
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
        ssgsea_clean_column_names()
    },
    error = function(e) NULL
  )
}

ssgsea_try_read_excel <- function(path) {
  
  tryCatch(
    {
      readxl::read_excel(
        path = path,
        sheet = 1,
          .name_repair = "minimal"
      ) %>%
        as.data.frame() %>%
        ssgsea_clean_column_names()
    },
    error = function(e) NULL
  )
}

ssgsea_score_table <- function(df) {
  
  if (is.null(df) || !is.data.frame(df) || ncol(df) <= 0) {
    return(-Inf)
  }
  
  non_empty_names <- sum(!is.na(colnames(df)) & colnames(df) != "")
  
  ncol(df) * 1000 + non_empty_names
}

ssgsea_smart_read_unchecked <- function(path, original_name = NULL, min_cols = 2) {
  
  ext_source <- original_name %||% path
  ext <- tolower(tools::file_ext(ext_source))
  
  if (ext %in% c("xls", "xlsx")) {
    
    validate(
      need(
        requireNamespace("readxl", quietly = TRUE),
        "Package 'readxl' is required to read Excel files."
      )
    )
    
    df <- ssgsea_try_read_excel(path)
    
    validate(
      need(
        !is.null(df) && ncol(df) >= min_cols,
        "Invalid Excel file. The first sheet must contain a valid table."
      )
    )
    
    return(as.data.frame(df))
  }
  
  read_attempts <- list(
    comma_dot = ssgsea_try_read_delim(
      path = path,
      delim = ",",
      decimal_mark = "."
    ),
    
    semicolon_dot = ssgsea_try_read_delim(
      path = path,
      delim = ";",
      decimal_mark = "."
    ),
    
    semicolon_comma = ssgsea_try_read_delim(
      path = path,
      delim = ";",
      decimal_mark = ","
    ),
    
    tab_dot = ssgsea_try_read_delim(
      path = path,
      delim = "\t",
      decimal_mark = "."
    ),
    
    tab_comma = ssgsea_try_read_delim(
      path = path,
      delim = "\t",
      decimal_mark = ","
    ),
    
    whitespace_dot = ssgsea_try_read_whitespace(
      path = path,
      decimal_mark = "."
    ),
    
    whitespace_comma = ssgsea_try_read_whitespace(
      path = path,
      decimal_mark = ","
    )
  )
  
  scores <- vapply(
    read_attempts,
    ssgsea_score_table,
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

ssgsea_smart_read <- function(...) {
  omics_validate_table(ssgsea_smart_read_unchecked(...))
}

open_ssgsea_plot_device <- function(file, fmt, width, height, dpi) {
  
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

format_p_value <- function(p) {

  ifelse(
    is.na(p),
    "NA",
    format.pval(p, digits = 3, eps = 0.001)
  )
}

format_p_label <- function(p, prefix) {

  formatted <- trimws(format_p_value(p))
  ifelse(startsWith(formatted, "<"),
    paste0(prefix, " ", formatted), paste0(prefix, " = ", formatted))
}


strip_ansi_codes <- function(x) {
  
  x <- as.character(x)
  
  # Remove ANSI color/style escape sequences such as ESC[38;5;232m.
  x <- gsub(
    pattern = "\033\\[[0-9;]*[[:alpha:]]",
    replacement = "",
    x = x,
    perl = TRUE
  )
  
  # Remove remaining non-printing control characters, preserving line breaks and tabs.
  x <- gsub(
    pattern = "[\001-\010\013\014\016-\037\177]",
    replacement = "",
    x = x,
    perl = TRUE
  )
  
  x
}


p_to_stars <- function(p) {
  
  case_when(
    is.na(p)  ~ "NA",
    p <= 1e-4 ~ "****",
    p <= 1e-3 ~ "***",
    p <= 1e-2 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ "ns"
  )
}


make_gene_set_overlap_diagnostics <- function(feature_ids, gene_sets) {
  
  feature_ids <- unique(as.character(feature_ids))
  feature_ids <- feature_ids[!is.na(feature_ids) & feature_ids != ""]
  
  gene_set_ids <- unique(as.character(unlist(gene_sets, use.names = FALSE)))
  gene_set_ids <- gene_set_ids[!is.na(gene_set_ids) & gene_set_ids != ""]
  
  exact_overlap <- length(intersect(feature_ids, gene_set_ids))
  
  case_insensitive_overlap <- length(
    intersect(
      stringr::str_to_upper(feature_ids),
      stringr::str_to_upper(gene_set_ids)
    )
  )
  
  list(
    exact_overlap = exact_overlap,
    case_insensitive_overlap = case_insensitive_overlap,
    n_features = length(feature_ids),
    n_gene_set_ids = length(gene_set_ids),
    example_features = head(feature_ids, 8),
    example_gene_set_ids = head(gene_set_ids, 8)
  )
}

format_gene_set_overlap_error <- function(diag, species = NULL) {
  
  species_txt <- species %||% "the selected species"
  
  base_msg <- paste0(
    "No exact overlap was found between the expression feature IDs and the selected gene sets. ",
    "The analysis was not run. ",
    "Selected species: ", species_txt, "."
  )
  
  if (!is.null(diag$case_insensitive_overlap) && diag$case_insensitive_overlap > 0) {
    
    base_msg <- paste0(
      base_msg,
      " A case-insensitive overlap was detected (n = ",
      diag$case_insensitive_overlap,
      "), suggesting an identifier capitalization mismatch. ",
      "Check whether the expression matrix uses the correct gene-symbol format for the selected species. ",
      "For example, human symbols are often uppercase, whereas mouse symbols commonly use title case."
    )
    
  } else {
    
    base_msg <- paste0(
      base_msg,
      " This usually indicates a species mismatch, an identifier type mismatch, or a gene-symbol formatting problem."
    )
  }
  
  paste0(
    base_msg,
    "

Expression feature examples: ",
    paste(diag$example_features, collapse = ", "),
    "
Gene-set identifier examples: ",
    paste(diag$example_gene_set_ids, collapse = ", ")
  )
}

get_default_box_fill_palette <- function() {
  c(
    "#AECBEF",
    "#F7C5B1",
    "#A8E0CC",
    "#F8DB9E",
    "#F6CCDC",
    "#9ECF9E",
    "#BAB4DD",
    "#F4B9B9"
  )
}

get_default_point_palette <- function() {
  c(
    "#2262AF",
    "#C0552A",
    "#168F64",
    "#C28400",
    "#BE6486",
    "#006B00",
    "#3C2F88",
    "#BA3B3B"
  )
}

parse_group_order <- function(group_order_text, detected_groups) {
  
  detected_groups <- as.character(detected_groups)
  detected_groups <- detected_groups[!is.na(detected_groups) & detected_groups != ""]
  
  if (length(detected_groups) == 0) {
    return(character(0))
  }
  
  if (is.null(group_order_text) || length(group_order_text) == 0 || is.na(group_order_text) || group_order_text == "") {
    return(detected_groups)
  }
  
  requested_groups <- unlist(
    strsplit(
      group_order_text,
      split = "[,;|]"
    )
  )
  
  requested_groups <- requested_groups %>%
    as.character() %>%
    stringr::str_trim()
  
  requested_groups <- requested_groups[!is.na(requested_groups) & requested_groups != ""]
  requested_groups <- requested_groups[requested_groups %in% detected_groups]
  requested_groups <- unique(requested_groups)
  
  missing_groups <- setdiff(detected_groups, requested_groups)
  
  c(requested_groups, missing_groups)
}

get_group_color_values <- function(input, detected_groups, prefix, default_palette) {
  
  detected_groups <- as.character(detected_groups)
  detected_groups <- detected_groups[!is.na(detected_groups) & detected_groups != ""]
  
  if (length(detected_groups) == 0) {
    return(character(0))
  }
  
  values <- vapply(
    seq_along(detected_groups),
    function(i) {
      input[[paste0(prefix, i)]] %||% rep(default_palette, length.out = length(detected_groups))[i]
    },
    FUN.VALUE = character(1)
  )
  
  stats::setNames(values, detected_groups)
}

apply_expression_transform <- function(mat, method, pseudocount = 1) {
  
  method <- method %||% "none"
  
  if (method == "none") {
    return(mat)
  }
  
  min_value <- suppressWarnings(min(mat, na.rm = TRUE))
  
  if (method == "log2_x") {
    
    validate(
      need(
        is.finite(min_value) && min_value > 0,
        "log2(x) requires all expression values to be > 0. Use log2(x + 1) if your matrix contains zeros."
      )
    )
    
    return(log2(mat))
  }
  
  if (method == "log2_x_plus_1") {
    
    validate(
      need(
        is.finite(min_value) && min_value >= 0,
        "log2(x + 1) requires non-negative values. If your data are already log-transformed or centered, choose 'No transformation'."
      )
    )
    
    return(log2(mat + 1))
  }
  
  if (method == "log2_x_plus_pseudocount") {
    
    validate(
      need(
        is.finite(pseudocount) && pseudocount > 0,
        "The pseudocount must be a positive numeric value."
      )
    )
    
    validate(
      need(
        is.finite(min_value) && min_value + pseudocount > 0,
        "log2(x + pseudocount) requires x + pseudocount to be > 0 for all values."
      )
    )
    
    return(log2(mat + pseudocount))
  }
  
  mat
}

collapse_duplicate_features <- function(expr_df, feature_ids, method) {
  
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
  
  if (method == "mean") {
    
    tmp <- tmp %>%
      group_by(.FeatureID) %>%
      summarise(
        across(
          all_of(sample_cols),
          ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)
        ),
        .groups = "drop"
      )
    
    mat <- as.matrix(tmp[, sample_cols, drop = FALSE])
    rownames(mat) <- tmp$.FeatureID
    
    return(mat)
  }
  
  if (method == "highest_mean") {
    
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
    
    return(mat)
  }
  
  if (method == "highest_variance") {
    
    tmp_mat <- as.matrix(tmp[, sample_cols, drop = FALSE])
    row_var <- apply(tmp_mat, 1, var, na.rm = TRUE)
    row_var[!is.finite(row_var)] <- -Inf
    
    tmp$.row_var <- row_var
    
    tmp <- tmp %>%
      group_by(.FeatureID) %>%
      slice_max(
        order_by = .row_var,
        n = 1,
        with_ties = FALSE
      ) %>%
      ungroup() %>%
      select(-.row_var)
    
    mat <- as.matrix(tmp[, sample_cols, drop = FALSE])
    rownames(mat) <- tmp$.FeatureID
    
    return(mat)
  }
  
  if (method == "first") {
    
    tmp <- tmp %>%
      distinct(.FeatureID, .keep_all = TRUE)
    
    mat <- as.matrix(tmp[, sample_cols, drop = FALSE])
    rownames(mat) <- tmp$.FeatureID
    
    return(mat)
  }
  
  tmp <- tmp %>%
    group_by(.FeatureID) %>%
    summarise(
      across(
        all_of(sample_cols),
        ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  
  mat <- as.matrix(tmp[, sample_cols, drop = FALSE])
  rownames(mat) <- tmp$.FeatureID
  
  mat
}

get_msigdb_collection <- function(species, category, subcategory = NULL) {
  omics_get_msigdb(species, category, subcategory)
}

add_y_positions <- function(pairwise_df, plot_df) {
  
  if (is.null(pairwise_df) || nrow(pairwise_df) == 0) {
    return(pairwise_df)
  }
  
  y_min <- min(plot_df$Score, na.rm = TRUE)
  y_max <- max(plot_df$Score, na.rm = TRUE)
  y_span <- y_max - y_min
  
  if (!is.finite(y_span) || y_span == 0) {
    y_span <- 1
  }
  
  pairwise_df %>%
    arrange(p_plot) %>%
    mutate(
      y.position = y_max + row_number() * 0.14 * y_span
    )
}

run_pathway_group_statistics <- function(df, p_adjust_method = "BH", show_ns = FALSE) {
  
  df <- df %>%
    filter(
      !is.na(Group),
      !is.na(Score),
      is.finite(Score)
    ) %>%
    mutate(
      Group = droplevels(factor(Group))
    )
  
  n_groups <- nlevels(df$Group)
  
  normality_table <- df %>%
    group_by(Group) %>%
    summarise(
      n = n(),
      mean = mean(Score, na.rm = TRUE),
      median = median(Score, na.rm = TRUE),
      sd = sd(Score, na.rm = TRUE),
      shapiro_p = ifelse(
        n() >= 3 && is.finite(sd(Score, na.rm = TRUE)) && sd(Score, na.rm = TRUE) > 0,
        tryCatch(
          stats::shapiro.test(Score)$p.value,
          error = function(e) NA_real_
        ),
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      normality_pass = !is.na(shapiro_p) & shapiro_p > 0.05
    )
  
  all_groups_testable_for_normality <- all(normality_table$n >= 3)
  all_normal <- all_groups_testable_for_normality &&
    all(!is.na(normality_table$shapiro_p)) &&
    all(normality_table$shapiro_p > 0.05)
  
  variance_table <- tryCatch(
    {
      rstatix::levene_test(df, Score ~ Group) %>%
        as_tibble()
    },
    error = function(e) {
      tibble(
        df1 = NA_real_,
        df2 = NA_real_,
        statistic = NA_real_,
        p = NA_real_
      )
    }
  )
  
  levene_p <- variance_table$p[1]
  variance_homogeneous <- !is.na(levene_p) && levene_p > 0.05
  
  if (n_groups < 2) {
    
    return(
      list(
        n_groups = n_groups,
        normality_table = normality_table,
        variance_table = variance_table,
        omnibus_table = tibble(
          test = "No test",
          statistic = NA_real_,
          p = NA_real_
        ),
        pairwise_table = tibble(),
        plot_pairwise_table = tibble(),
        omnibus_label = "Only one group available: no statistical test performed",
        assumption_caption = "At least two groups are required for group comparison."
      )
    )
  }
  
  if (n_groups == 2) {
    
    if (isTRUE(all_normal)) {
      
      test_name <- if (isTRUE(variance_homogeneous)) {
        "Student t-test"
      } else {
        "Welch t-test"
      }
      
      tt <- stats::t.test(
        Score ~ Group,
        data = df,
        var.equal = isTRUE(variance_homogeneous)
      )
      
      group_levels <- levels(df$Group)
      
      omnibus_table <- tibble(
        test = test_name,
        statistic = unname(tt$statistic),
        p = tt$p.value
      )
      
      pairwise_table <- tibble(
        group1 = group_levels[1],
        group2 = group_levels[2],
        test = test_name,
        p = tt$p.value,
        p.adj = tt$p.value,
        p.signif = p_to_stars(tt$p.value),
        p.label = format_p_label(tt$p.value, "P"),
        p_plot = tt$p.value
      )
      
    } else {
      
      wt <- stats::wilcox.test(
        Score ~ Group,
        data = df,
        exact = FALSE
      )
      
      group_levels <- levels(df$Group)
      
      omnibus_table <- tibble(
        test = "Mann-Whitney/Wilcoxon rank-sum test",
        statistic = unname(wt$statistic),
        p = wt$p.value
      )
      
      pairwise_table <- tibble(
        group1 = group_levels[1],
        group2 = group_levels[2],
        test = "Mann-Whitney/Wilcoxon rank-sum test",
        p = wt$p.value,
        p.adj = wt$p.value,
        p.signif = p_to_stars(wt$p.value),
        p.label = format_p_label(wt$p.value, "P"),
        p_plot = wt$p.value
      )
    }
    
  } else {
    
    if (isTRUE(all_normal) && isTRUE(variance_homogeneous)) {
      
      fit <- stats::aov(Score ~ Group, data = df)
      sm <- summary(fit)[[1]]
      
      omnibus_table <- tibble(
        test = "One-way ANOVA",
        statistic = sm$`F value`[1],
        p = sm$`Pr(>F)`[1]
      )
      
      pairwise_table <- tryCatch(
        {
          rstatix::tukey_hsd(df, Score ~ Group) %>%
            as_tibble() %>%
            mutate(
              test = "Tukey post hoc",
              p = NA_real_,
              adjustment = "Tukey; pairwise comparisons within this pathway",
              p.signif = p_to_stars(p.adj),
              p.label = format_p_label(p.adj, "Adj. P"),
              p_plot = p.adj
            ) %>%
            select(group1, group2, test, p, p.adj, p.signif, p.label, p_plot, everything())
        },
        error = function(e) validate(need(FALSE,
          paste0("Tukey post hoc failed: ", conditionMessage(e))))
      )
      
    } else if (isTRUE(all_normal) && !isTRUE(variance_homogeneous)) {
      
      welch <- tryCatch(
        {
          rstatix::welch_anova_test(df, Score ~ Group) %>%
            as_tibble()
        },
        error = function(e) tibble(statistic = NA_real_, p = NA_real_)
      )
      
      omnibus_table <- tibble(
        test = "Welch ANOVA",
        statistic = welch$statistic[1],
        p = welch$p[1]
      )
      
      pairwise_table <- tryCatch(
        {
          rstatix::games_howell_test(df, Score ~ Group) %>%
            as_tibble() %>%
            mutate(
              test = "Games-Howell post hoc",
              p = NA_real_,
              adjustment = "Games-Howell; pairwise comparisons within this pathway",
              p.signif = p_to_stars(p.adj),
              p.label = format_p_label(p.adj, "Adj. P"),
              p_plot = p.adj
            ) %>%
            select(group1, group2, test, p, p.adj, p.signif, p.label, p_plot, everything())
        },
        error = function(e) validate(need(FALSE,
          paste0("Games-Howell post hoc failed: ", conditionMessage(e))))
      )
      
    } else {
      
      kw <- stats::kruskal.test(Score ~ Group, data = df)
      
      omnibus_table <- tibble(
        test = "Kruskal-Wallis test",
        statistic = unname(kw$statistic),
        p = kw$p.value
      )
      
      pairwise_table <- tryCatch(
        {
          rstatix::dunn_test(
            df,
            Score ~ Group,
            p.adjust.method = p_adjust_method
          ) %>%
            as_tibble() %>%
            mutate(
              test = paste0("Dunn post hoc; p.adjust = ", p_adjust_method),
              adjustment = paste0(p_adjust_method, "; pairwise comparisons within this pathway"),
              p.signif = p_to_stars(p.adj),
              p.label = format_p_label(p.adj, "Adj. P"),
              p_plot = p.adj
            ) %>%
            select(group1, group2, test, p, p.adj, p.signif, p.label, p_plot, everything())
        },
        error = function(e) validate(need(FALSE,
          paste0("Dunn post hoc failed: ", conditionMessage(e))))
      )
    }
  }
  
  plot_pairwise_table <- pairwise_table
  
  if (!isTRUE(show_ns) && nrow(plot_pairwise_table) > 0) {
    plot_pairwise_table <- plot_pairwise_table %>%
      filter(!is.na(p_plot), p_plot <= 0.05)
  }
  
  plot_pairwise_table <- add_y_positions(
    pairwise_df = plot_pairwise_table,
    plot_df = df
  )
  
  normality_txt <- normality_table %>%
    mutate(
      shapiro_txt = ifelse(
        is.na(shapiro_p),
        "NA",
        formatC(shapiro_p, digits = 3, format = "f")
      )
    ) %>%
    transmute(txt = paste0("Shapiro P (", Group, ") = ", shapiro_txt)) %>%
    pull(txt) %>%
    paste(collapse = "\n")
  
  variance_txt <- format_p_label(levene_p, "Levene P")
  
  assumption_caption <- paste0(
    "Normality by group:\n",
    normality_txt,
    "\nHomogeneity of variances: ",
    variance_txt
  )
  
  omnibus_label <- paste0(
    omnibus_table$test[1],
    ": ",
    format_p_label(omnibus_table$p[1], "P")
  )
  
  list(
    n_groups = n_groups,
    normality_table = normality_table,
    variance_table = variance_table,
    omnibus_table = omnibus_table,
    pairwise_table = pairwise_table,
    plot_pairwise_table = plot_pairwise_table,
    omnibus_label = omnibus_label,
    assumption_caption = assumption_caption
  )
}


run_pathway_omnibus_only <- function(df) {
  
  df <- df %>%
    filter(
      !is.na(Group),
      !is.na(Score),
      is.finite(Score)
    ) %>%
    mutate(
      Group = droplevels(factor(Group))
    )
  
  n_groups <- nlevels(df$Group)
  
  if (n_groups < 2) {
    return(
      tibble(
        test = "No test",
        statistic = NA_real_,
        p = NA_real_,
        n_groups = n_groups
      )
    )
  }
  
  normality_table <- df %>%
    group_by(Group) %>%
    summarise(
      n = n(),
      sd = sd(Score, na.rm = TRUE),
      shapiro_p = ifelse(
        n() >= 3 && is.finite(sd(Score, na.rm = TRUE)) && sd(Score, na.rm = TRUE) > 0,
        tryCatch(
          stats::shapiro.test(Score)$p.value,
          error = function(e) NA_real_
        ),
        NA_real_
      ),
      .groups = "drop"
    )
  
  all_groups_testable_for_normality <- all(normality_table$n >= 3)
  all_normal <- all_groups_testable_for_normality &&
    all(!is.na(normality_table$shapiro_p)) &&
    all(normality_table$shapiro_p > 0.05)
  
  variance_table <- tryCatch(
    {
      rstatix::levene_test(df, Score ~ Group) %>%
        as_tibble()
    },
    error = function(e) {
      tibble(p = NA_real_)
    }
  )
  
  levene_p <- variance_table$p[1]
  variance_homogeneous <- !is.na(levene_p) && levene_p > 0.05
  
  if (n_groups == 2) {
    
    if (isTRUE(all_normal)) {
      
      test_name <- if (isTRUE(variance_homogeneous)) {
        "Student t-test"
      } else {
        "Welch t-test"
      }
      
      tt <- tryCatch(
        stats::t.test(
          Score ~ Group,
          data = df,
          var.equal = isTRUE(variance_homogeneous)
        ),
        error = function(e) NULL
      )
      
      return(
        tibble(
          test = test_name,
          statistic = if (!is.null(tt)) unname(tt$statistic) else NA_real_,
          p = if (!is.null(tt)) tt$p.value else NA_real_,
          n_groups = n_groups
        )
      )
      
    } else {
      
      wt <- tryCatch(
        stats::wilcox.test(
          Score ~ Group,
          data = df,
          exact = FALSE
        ),
        error = function(e) NULL
      )
      
      return(
        tibble(
          test = "Mann-Whitney/Wilcoxon rank-sum test",
          statistic = if (!is.null(wt)) unname(wt$statistic) else NA_real_,
          p = if (!is.null(wt)) wt$p.value else NA_real_,
          n_groups = n_groups
        )
      )
    }
  }
  
  if (isTRUE(all_normal) && isTRUE(variance_homogeneous)) {
    
    fit <- tryCatch(
      stats::aov(Score ~ Group, data = df),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      return(
        tibble(
          test = "One-way ANOVA",
          statistic = NA_real_,
          p = NA_real_,
          n_groups = n_groups
        )
      )
    }
    
    sm <- summary(fit)[[1]]
    
    return(
      tibble(
        test = "One-way ANOVA",
        statistic = sm$`F value`[1],
        p = sm$`Pr(>F)`[1],
        n_groups = n_groups
      )
    )
  }
  
  if (isTRUE(all_normal) && !isTRUE(variance_homogeneous)) {
    
    welch <- tryCatch(
      {
        rstatix::welch_anova_test(df, Score ~ Group) %>%
          as_tibble()
      },
      error = function(e) tibble(statistic = NA_real_, p = NA_real_)
    )
    
    return(
      tibble(
        test = "Welch ANOVA",
        statistic = welch$statistic[1],
        p = welch$p[1],
        n_groups = n_groups
      )
    )
  }
  
  kw <- tryCatch(
    stats::kruskal.test(Score ~ Group, data = df),
    error = function(e) NULL
  )
  
  tibble(
    test = "Kruskal-Wallis test",
    statistic = if (!is.null(kw)) unname(kw$statistic) else NA_real_,
    p = if (!is.null(kw)) kw$p.value else NA_real_,
    n_groups = n_groups
  )
}


# =========================================================
# Pathway-wide FDR is calculated before display filtering or Top N selection.
# Failed/undefined tests retain NA and are not treated as valid hypotheses.
adjust_pathway_pvalues <- function(pvalues) {
  adjusted <- rep(NA_real_, length(pvalues))
  valid <- is.finite(pvalues) & pvalues >= 0 & pvalues <= 1
  adjusted[valid] <- stats::p.adjust(pvalues[valid], method = "BH")
  adjusted
}

filter_pathway_statistics <- function(stats_df, source = "adjusted", cutoff = 0.05) {
  validate(need(source %in% c("adjusted", "nominal"), "Invalid heatmap P-value source."))
  omics_number_parameter(cutoff, "Heatmap P-value/FDR cutoff", 0, 1)
  metric <- if (source == "adjusted") "padj" else "pvalue"
  stats_df %>%
    filter(!is.na(.data[[metric]]), .data[[metric]] <= cutoff) %>%
    arrange(.data[[metric]], desc(group_mean_variance), desc(score_range))
}

# UI
# =========================================================

ssgseaModuleUI <- function(id) {
  
  ns <- NS(id)
  
  tagList(
    
    useShinyjs(),
    
    tags$head(
      tags$style(HTML("
        .ssgsea-panel {
          background: #ffffff;
          border-radius: 8px;
          padding: 20px;
          border: 1px solid #eeeeee;
          margin-bottom: 20px;
        }

        .ssgsea-help {
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

        .selectize-control.single .selectize-input > div {
          max-width: calc(100% - 28px);
          overflow: hidden;
          text-overflow: ellipsis;
          white-space: nowrap;
        }
      "))
    ),
    
    fluidRow(
      
      column(
        width = 4,
        
        wellPanel(
          
          # -------------------------------------------------
          # 1) Data Input
          # -------------------------------------------------
          
          h4("1. Data Input"),
          
          fileInput(
            ns("upload"),
            "Upload expression matrix (.txt, .csv, .tsv, .xls, .xlsx)",
            accept = c(".txt", ".csv", ".tsv", ".xls", ".xlsx")
          ),
          
          uiOutput(ns("column_selector")),
          
          p(
            class = "ssgsea-help",
            "Input should contain one feature identifier column and numeric sample columns. Duplicated feature identifiers are handled internally by retaining the row with the highest mean across samples."
          ),
          
          hr(),
          
          # -------------------------------------------------
          # 2) Sample Annotation
          # -------------------------------------------------
          
          h4("2. Sample Annotation"),
          
          fileInput(
            ns("anno_upload"),
            "Upload sample annotation (.txt, .csv, .tsv, .xls, .xlsx)",
            accept = c(".txt", ".csv", ".tsv", ".xls", ".xlsx")
          ),
          
          uiOutput(ns("anno_column_selector")),
          
          p(
            class = "ssgsea-help",
            "Optional. If provided, samples will be grouped in the heatmap and used for the selected pathway boxplot."
          ),
          
          hr(),
          
          # -------------------------------------------------
          # 4) Gene Set Source
          # -------------------------------------------------
          
          h4("3. Gene Set Source"),
          
          radioButtons(
            ns("geneset_source"),
            "Gene set source:",
            choices = c(
              "MSigDB" = "msigdb",
              "Custom GMT file" = "custom"
            ),
            selected = "msigdb"
          ),
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'msigdb'", ns("geneset_source")),
            
            omics_msigdb_ui(ns)
          ),
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'custom'", ns("geneset_source")),
            
            fileInput(
              ns("gmt_upload"),
              "Upload custom gene sets (.gmt)",
              accept = c(".gmt", ".txt")
            ),
            p(class = "help-block", "Input IDs must match the GMT exactly. No gene/protein IDs are converted.")
          ),
          
          hr(),
          
          # -------------------------------------------------
          # 5) ssGSEA Parameters
          # -------------------------------------------------
          
          h4("4. ssGSEA Parameters"),
          
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
            ns("alpha"),
            "Alpha:",
            value = 0.25,
            min = 0,
            step = 0.05
          ),

          p(
            class = "ssgsea-help",
            "Alpha is the ssGSEA weighting exponent (Barbie et al. 2009 'tau'). Higher values emphasize genes at the extremes of each sample's expression ranking more heavily; 0 gives an unweighted (rank-only) score. Default: 0.25."
          ),

          checkboxInput(
            ns("normalize_scores"),
            "Normalize ssGSEA scores",
            value = TRUE
          ),
          
          actionButton(
            ns("run_ssgsea"),
            "Run ssGSEA",
            class = "btn-primary",
            style = "width:100%; font-weight:bold;"
          ),
          
          hr(),
          
          # -------------------------------------------------
          # 6) Plot and Statistics Options
          # -------------------------------------------------
          
          h4("5. Plot and Statistics Options"),
          
          p(
            class = "ssgsea-help",
            "Heatmap filtering uses the global comparison of ssGSEA scores between groups. BH FDR is calculated across all tested pathways before display selection; nominal P-values are available for exploration."
          ),
          
          sliderInput(
            ns("top_n"),
            "Number of pathways to display:",
            min = 5,
            max = 100,
            value = 30,
            step = 5
          ),
          
          radioButtons(
            ns("heatmap_p_source"),
            "Use for heatmap pathway filtering:",
            choices = c("Adjusted P-value / FDR (BH)" = "adjusted", "Nominal P-value" = "nominal"),
            selected = "adjusted"
          ),

          numericInput(
            ns("heatmap_pvalue_cutoff"),
            "Heatmap P-value / FDR cutoff:",
            value = 0.05,
            min = 0,
            max = 1,
            step = 0.01
          ),
          
          textInput(
            ns("plot_title"),
            "Plot Title",
            value = "ssGSEA heatmap"
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
            ns("heatmap_wrap_width"),
            "Pathway label wrap width:",
            min = 20,
            max = 80,
            value = 35,
            step = 5
          ),
          
          checkboxInput(
            ns("scale_rows"),
            "Z-score pathways by row in heatmap",
            value = TRUE
          ),
          
          checkboxInput(
            ns("show_heatmap_sample_names"),
            "Show sample names in heatmap",
            value = FALSE
          ),
          
          uiOutput(ns("posthoc_adjust_ui")),
          hr(),
          
          # -------------------------------------------------
          # 7) Colors
          # -------------------------------------------------
          
          h4("6. Colors"),
          
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

          colourInput(
            ns("profile_col"),
            "Selected pathway profile",
            "#2A78D6"
          ),
          
          hr(),
          
          h5("Boxplot group order"),
          
          uiOutput(ns("group_order_ui")),
          
          p(
            class = "ssgsea-help",
            "To change the boxplot order, type the group names separated by commas. Missing groups are automatically added at the end."
          ),
          
          hr(),
          
          h5("Boxplot colors"),
          
          uiOutput(ns("boxplot_color_ui")),
          
          hr(),
          
          h5("Boxplot line settings"),
          
          sliderInput(
            ns("box_linewidth"),
            "Boxplot line width:",
            min = 0.10,
            max = 1.00,
            value = 0.25,
            step = 0.05
          ),
          
          sliderInput(
            ns("point_stroke"),
            "Point outline width:",
            min = 0.00,
            max = 1.00,
            value = 0.20,
            step = 0.05
          ),
          
          sliderInput(
            ns("bracket_linewidth"),
            "Significance bracket line width:",
            min = 0.10,
            max = 1.00,
            value = 0.20,
            step = 0.05
          ),
          
          hr(),
          
          # -------------------------------------------------
          # 8) Export
          # -------------------------------------------------
          
          h4("7. Export"),

          h5("Heatmap"),

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
                value = 18,
                min = 5,
                max = 60
              )
            )
          ),

          h5("Selected pathway boxplot"),

          fluidRow(
            column(
              6,
              numericInput(
                ns("exp_width_box"),
                "Width (cm)",
                value = 14,
                min = 5,
                max = 60
              )
            ),

            column(
              6,
              numericInput(
                ns("exp_height_box"),
                "Height (cm)",
                value = 10,
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
            ns("download_boxplot"),
            "Download selected pathway boxplot",
            class = "btn-success",
            style = "width:100%"
          ),
          
          br(),
          br(),
          
          downloadButton(
            ns("download_table"),
            "Download ssGSEA matrix",
            class = "btn-info",
            style = "width:100%"
          ),
          
          br(),
          br(),
          
          downloadButton(
            ns("download_processed_matrix"),
            "Download processed expression matrix",
            class = "btn-secondary",
            style = "width:100%"
          ),
          
          br(),
          br(),
          
          downloadButton(
            ns("download_pathway_statistics"),
            "Download pathway P-values and FDR",
            class = "btn-secondary",
            style = "width:100%"
          )
        )
      ),
      
      column(
        width = 8,
        
        div(
          class = "ssgsea-panel",
          h4("Input Summary"),
          verbatimTextOutput(ns("input_summary"))
        ),
        
        div(
          class = "ssgsea-panel",
          h4("ssGSEA Summary"),
          verbatimTextOutput(ns("ssgsea_summary"))
        ),
        
        div(
          class = "ssgsea-panel",
          h4("ssGSEA Heatmap"),
          plotOutput(ns("ssgsea_heatmap"), height = "700px")
        ),
        
        div(
          class = "ssgsea-panel",
          h4("Selected Pathway Boxplot / Profile"),
          uiOutput(ns("pathway_selector_ui")),
          br(),
          plotOutput(ns("pathway_profile"), height = "420px")
        ),
        
        div(
          class = "ssgsea-panel",
          h4("Selected Pathway Statistical Summary"),
          verbatimTextOutput(ns("pathway_stats_summary")),
        ),
        
        div(
          class = "ssgsea-panel",
          h4("ssGSEA Score Table"),
          DTOutput(ns("ssgsea_table"))
        )
      )
    )
  )
}

# =========================================================
# SERVER
# =========================================================

ssgseaModuleServer <- function(id) {
  
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    # =====================================================
    # 1) Raw Data
    # =====================================================
    
    raw_data <- reactive({
      
      req(input$upload)
      
      ssgsea_smart_read(
        path = input$upload$datapath,
        original_name = input$upload$name,
        min_cols = 2
      )
    })
    
    # =====================================================
    # 2) Column Selector
    # =====================================================
    
    output$column_selector <- renderUI({
      
      req(raw_data())
      
      cols <- colnames(raw_data())
      
      selectInput(
        ns("col_id"),
        "Gene/protein/feature ID column:",
        choices = cols,
        selected = cols[1]
      )
    })
    
    # =====================================================
    # 2) Sample Annotation Selectors
    # =====================================================
    
    output$anno_column_selector <- renderUI({
      
      req(input$anno_upload)
      
      anno_df <- ssgsea_smart_read(
        path = input$anno_upload$datapath,
        original_name = input$anno_upload$name,
        min_cols = 2
      )
      
      cols <- colnames(anno_df)
      
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
          choices = cols,
          selected = if (length(cols) >= 2) cols[2] else cols[1]
        )
      )
    })
    
    # =====================================================
    # 3) MSigDB Subcategory Mapping
    # =====================================================
    
    get_safe_subcategory <- omics_msigdb_server(input, session, "geneset_source")
    
    # =====================================================
    # 4) Expression Matrix Processing
    # =====================================================
    
    expr_matrix <- reactive({
      
      req(
        raw_data(),
        input$col_id
      )
      
      df <- raw_data()
      
      validate(
        need(
          input$col_id %in% colnames(df),
          "Selected feature ID column not found."
        )
      )
      
      feature_ids <- df[[input$col_id]]
      
      expr_df <- df[, setdiff(colnames(df), input$col_id), drop = FALSE]
      
      expr_df <- omics_numeric_frame(expr_df)
      
      validate(
        need(
          ncol(expr_df) >= 1,
          "No numeric sample columns were found."
        )
      )
      
      mat <- collapse_duplicate_features(
        expr_df = expr_df,
        feature_ids = feature_ids,
        method = "highest_mean"
      )
      
      storage.mode(mat) <- "numeric"
      
      # ---------------------------------------------------
      # Remove missing / non-finite rows
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
          nrow(mat) > 1,
          "Too few genes/features remain after removing missing or non-finite values."
        ),
        need(
          ncol(mat) > 0,
          "No valid samples remain after processing."
        )
      )
      
      mat
    })
    
    # =====================================================
    # 5) Sample Annotation Processing
    # =====================================================
    
    sample_annotation <- reactive({
      
      if (is.null(input$anno_upload)) {
        return(NULL)
      }
      
      req(
        expr_matrix(),
        input$anno_sample_col,
        input$anno_group_col
      )
      
      anno_df <- ssgsea_smart_read(
        path = input$anno_upload$datapath,
        original_name = input$anno_upload$name,
        min_cols = 2
      )
      
      validate(
        need(
          input$anno_sample_col %in% colnames(anno_df),
          "Selected sample column not found in annotation."
        ),
        need(
          input$anno_group_col %in% colnames(anno_df),
          "Selected group column not found in annotation."
        )
      )
      
      anno_df <- omics_annotation(anno_df, input$anno_sample_col,
        colnames(expr_matrix()), required = input$anno_group_col)

      anno_df <- anno_df %>%
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
      
      validate(
        need(
          nrow(anno_df) > 0,
          "Sample annotation file is empty after filtering."
        )
      )
      
      expr_samples <- colnames(expr_matrix())
      
      anno_df <- tibble(Sample = expr_samples) %>%
        left_join(
          anno_df,
          by = "Sample"
        ) %>%
        mutate(
          Group = ifelse(
            is.na(Group) | Group == "",
            "Unannotated",
            Group
          )
        )
      
      anno_df$Sample <- factor(
        anno_df$Sample,
        levels = expr_samples
      )
      
      anno_df %>%
        arrange(Group, Sample) %>%
        mutate(
          Group = factor(
            Group,
            levels = unique(Group)
          )
        )
    })
    
    detected_groups <- reactive({
      
      anno <- sample_annotation()
      
      if (is.null(anno) || !"Group" %in% colnames(anno)) {
        return(character(0))
      }
      
      groups <- as.character(anno$Group)
      groups <- groups[!is.na(groups) & groups != ""]
      
      unique(groups)
    })
    
    group_order_levels <- reactive({
      
      parse_group_order(
        group_order_text = input$group_order_text,
        detected_groups = detected_groups()
      )
    })
    
    sample_order <- reactive({
      
      anno <- sample_annotation()
      
      if (is.null(anno)) {
        return(colnames(expr_matrix()))
      }
      
      if ("Group" %in% colnames(anno)) {
        
        group_levels <- group_order_levels()
        
        anno <- anno %>%
          mutate(
            Group = factor(
              as.character(Group),
              levels = group_levels
            )
          ) %>%
          arrange(
            Group,
            Sample
          )
        
        return(as.character(anno$Sample))
      }
      
      as.character(anno$Sample)
    })
    
    output$group_order_ui <- renderUI({
      
      groups <- detected_groups()
      
      if (length(groups) == 0) {
        return(
          p(
            class = "ssgsea-help",
            "Upload a sample annotation file and select a group column to control heatmap and boxplot order."
          )
        )
      }
      
      textInput(
        ns("group_order_text"),
        "Group order for heatmap and boxplot:",
        value = paste(groups, collapse = ", ")
      )
    })
    
    output$posthoc_adjust_ui <- renderUI({
      
      groups <- detected_groups()
      
      if (length(groups) <= 2) {
        return(NULL)
      }
      
      selectInput(
        ns("posthoc_adjust_method"),
        "Dunn post hoc P-value adjustment:",
        choices = c(
          "BH",
          "bonferroni",
          "holm",
          "hochberg",
          "hommel",
          "BY"
        ),
        selected = "BH"
      )
    })
    
    output$boxplot_color_ui <- renderUI({
      
      groups <- detected_groups()
      
      if (length(groups) == 0) {
        return(
          p(
            class = "ssgsea-help",
            "Upload a sample annotation file and select a group column to customize boxplot colors."
          )
        )
      }
      
      box_defaults <- rep(
        get_default_box_fill_palette(),
        length.out = length(groups)
      )
      
      point_defaults <- rep(
        get_default_point_palette(),
        length.out = length(groups)
      )
      
      tagList(
        lapply(
          seq_along(groups),
          function(i) {
            fluidRow(
              column(
                6,
                colourInput(
                  ns(paste0("box_fill_", i)),
                  paste0("Box fill — ", groups[i]),
                  box_defaults[i]
                )
              ),
              
              column(
                6,
                colourInput(
                  ns(paste0("point_col_", i)),
                  paste0("Points — ", groups[i]),
                  point_defaults[i]
                )
              )
            )
          }
        )
      )
    })
    
    # =====================================================
    # 6) Gene Set Loading
    # =====================================================
    
    gene_sets <- reactive({
      
      req(input$geneset_source)
      
      if (input$geneset_source == "msigdb") {
        
        req(
          input$species,
          input$category
        )
        
        msig <- get_msigdb_collection(
          species = input$species,
          category = input$category,
          subcategory = get_safe_subcategory()
        )
        
        validate(
          need(
            "gene_symbol" %in% colnames(msig),
            "MSigDB output does not contain the expected gene_symbol column."
          )
        )
        
        msig <- msig %>%
          filter(
            !is.na(gene_symbol),
            gene_symbol != ""
          ) %>%
          distinct(
            gs_name,
            gene_symbol
          )
        
        gsets <- split(
          msig$gene_symbol,
          msig$gs_name
        )
        
        gsets <- lapply(gsets, unique)
        gsets <- gsets[lengths(gsets) > 0]
        
        validate(
          need(
            length(gsets) > 0,
            "No valid MSigDB gene sets remained after filtering."
          )
        )
        
        gsets
        
      } else {
        
        req(input$gmt_upload)
        
        gsets <- omics_read_gmt(input$gmt_upload$datapath)
        
        gsets <- lapply(gsets, unique)
        gsets <- gsets[lengths(gsets) > 0]
        
        validate(
          need(
            length(gsets) > 0,
            "No valid gene sets found in the uploaded GMT file."
          )
        )
        
        gsets
      }
    })
    
    overlap_diagnostics <- reactive({
      
      req(
        expr_matrix(),
        gene_sets()
      )
      
      make_gene_set_overlap_diagnostics(
        feature_ids = rownames(expr_matrix()),
        gene_sets = gene_sets()
      )
    })
    
    overlap_n <- reactive({
      
      overlap_diagnostics()$exact_overlap
    })
    
    # =====================================================
    # 7) Input Summary
    # =====================================================
    
    output$input_summary <- renderPrint({
      
      req(
        expr_matrix(),
        gene_sets(),
        overlap_n()
      )
      
      mat <- expr_matrix()
      gsets <- gene_sets()
      anno <- sample_annotation()
      
      cat("Features:", nrow(mat), "\n")
      cat("Samples:", ncol(mat), "\n")
      cat("Gene sets loaded:", length(gsets), "\n")
      diag <- overlap_diagnostics()
      
      cat("Feature overlap with gene sets:", diag$exact_overlap, "\n")
      cat("Case-insensitive overlap with gene sets:", diag$case_insensitive_overlap, "\n")
      
      if (diag$exact_overlap == 0) {
        cat("\nWARNING:\n")
        cat(format_gene_set_overlap_error(diag, species = input$species), "\n")
        cat("\n")
      }
      
      cat("Duplicated feature handling: keep row with highest mean across samples\n")
      cat("Expression transformation: input matrix must already be log2-transformed\n")
      cat(
        "Features removed due to missing/non-finite values:",
        attr(mat, "removed_missing_rows") %||% 0,
        "\n"
      )
      cat("Expression range:\n")
      print(range(mat, na.rm = TRUE))
      cat("\n")
      
      if (input$geneset_source == "msigdb") {
        cat("Gene set source: MSigDB\n")

        cat("Species:", input$species, "\n")
        cat("Category:", input$category, "\n")
        cat("Subcategory:", get_safe_subcategory() %||% "None", "\n")
      } else {
        cat("Gene set source: Custom GMT\n")
        cat("GMT file:", input$gmt_upload$name, "\n")
      }
      
      cat("\n")
      
      if (is.null(anno)) {
        cat("Sample annotation: not provided\n")
      } else {
        cat("Sample annotation: provided\n")
        cat("Annotated samples:", nrow(anno), "\n")
        cat("Groups:", nlevels(anno$Group), "\n")
        print(table(anno$Group))
      }
    })
    
    # =====================================================
    # 8) Run ssGSEA
    # =====================================================
    
    ssgsea_scores_state <- reactive({
      list(matrix = expr_matrix(), source = input$geneset_source, gmt = input$gmt_upload, species = input$species, category = input$category, subcategory = get_safe_subcategory(),
        min_size = input$min_size, max_size = input$max_size,
        alpha = input$alpha, normalize = input$normalize_scores)
    })

    ssgsea_scores <- omics_event_reactive(input$run_ssgsea, {
      omics_set_sizes(input$min_size, input$max_size)
      omics_number_parameter(input$alpha, "ssGSEA alpha", 0, Inf)

      
      req(
        expr_matrix(),
        gene_sets()
      )
      
      mat <- expr_matrix()
      gsets <- gene_sets()
      omics_validate_gene_set_ids(rownames(mat), gsets, "Expression matrix")
      ov <- overlap_n()
      
      validate(
        need(
          input$min_size >= 1,
          "Minimum set size must be >= 1."
        ),
        need(
          input$max_size >= input$min_size,
          "Maximum set size must be >= minimum set size."
        ),
        need(
          ov > 0,
          format_gene_set_overlap_error(
            overlap_diagnostics(),
            species = input$species
          )
        )
      )
      
      if (ov <= 0) {
        showModal(
          modalDialog(
            title = "No overlap between expression IDs and gene sets",
            div(
              style = "white-space: pre-wrap;",
              format_gene_set_overlap_error(
                overlap_diagnostics(),
                species = input$species
              )
            ),
            easyClose = TRUE,
            footer = modalButton("OK")
          )
        )
        
        validate(
          need(
            FALSE,
            format_gene_set_overlap_error(
              overlap_diagnostics(),
              species = input$species
            )
          )
        )
      }
      
      withProgress(message = "Running ssGSEA...", value = 0.5, {
        
        es <- tryCatch(
          {
            param <- GSVA::ssgseaParam(
              exprData = mat,
              geneSets = gsets,
              alpha = input$alpha,
              normalize = isTRUE(input$normalize_scores),
              minSize = input$min_size,
              maxSize = input$max_size,
              checkNA = "auto",
              use = "na.rm"
            )
            
            GSVA::gsva(
              param,
              verbose = FALSE
            )
          },
          error = function(e) {
            e
          }
        )
      })
      
      if (inherits(es, "error")) {
        
        diag <- overlap_diagnostics()
        raw_error <- strip_ansi_codes(conditionMessage(es))
        
        if (grepl("No gene set left after mapping and filtering", raw_error, ignore.case = TRUE)) {
          
          error_msg <- paste0(
            "No gene set remained after ID mapping and set-size filtering.\n\n",
            "Check species/ID format or change Min/Max set size.\n\n",
            "Exact overlap: ",
            diag$exact_overlap,
            " | Case-insensitive overlap: ",
            diag$case_insensitive_overlap
          )
          
        } else {
          
          error_msg <- paste0(
            "ssGSEA failed.\n\n",
            raw_error,
            "\n\nExact overlap: ",
            diag$exact_overlap,
            " | Case-insensitive overlap: ",
            diag$case_insensitive_overlap
          )
        }
        
        showModal(
          modalDialog(
            title = "ssGSEA failed",
            div(
              style = "white-space: pre-wrap;",
              error_msg
            ),
            easyClose = TRUE,
            footer = modalButton("OK")
          )
        )
        
        validate(
          need(FALSE, error_msg)
        )
      }
      
      as.matrix(es)
    }, state = ssgsea_scores_state, label = "ssGSEA")
    
    output$ssgsea_summary <- renderPrint({
      
      req(ssgsea_scores())
      
      es <- ssgsea_scores()
      
      cat("Pathways scored:", nrow(es), "\n")
      cat("Samples scored:", ncol(es), "\n")
      cat("Heatmap pathway selection:", "significant", "\n")

      if (is.null(sample_annotation())) {
        cat("Pathways displayed in heatmap: unavailable (no sample annotation provided)\n")
      } else {
        cat("Pathways displayed in heatmap:", nrow(heatmap_matrix()), "\n")
      }

      cat("Heatmap filter:", input$heatmap_p_source %||% "adjusted", "\n")
      cat("Heatmap P-value/FDR cutoff:", input$heatmap_pvalue_cutoff %||% 0.05, "\n")

      cat("Score range:\n")
      print(range(es, na.rm = TRUE))
    })
    
    # =====================================================
    # 9) Heatmap Pathway Selection
    # =====================================================
    
    heatmap_pathway_stats <- reactive({
      
      req(
        ssgsea_scores(),
        sample_annotation()
      )
      
      es <- ssgsea_scores()
      anno <- sample_annotation()
      
      validate(
        need(
          "Group" %in% colnames(anno),
          "Sample annotation with a group column is required to select statistically different pathways."
        )
      )
      
      long_df <- as.data.frame(es) %>%
        rownames_to_column("Pathway") %>%
        pivot_longer(
          cols = -Pathway,
          names_to = "Sample",
          values_to = "Score"
        ) %>%
        left_join(
          anno %>%
            mutate(Sample = as.character(Sample)) %>%
            select(Sample, Group),
          by = "Sample"
        ) %>%
        filter(
          !is.na(Group),
          !is.na(Score),
          is.finite(Score)
        )
      
      validate(
        need(
          nrow(long_df) > 0,
          "No valid ssGSEA scores and group annotations are available for pathway statistics."
        ),
        need(
          dplyr::n_distinct(long_df$Group) >= 2,
          "At least two sample groups are required to select statistically different pathways."
        )
      )
      
      split_df <- split(long_df, long_df$Pathway)
      
      stat_df <- purrr::map_dfr(
        names(split_df),
        function(pathway_name) {
          
          df <- split_df[[pathway_name]] %>%
            mutate(Group = factor(Group)) %>%
            filter(!is.na(Group))
          
          omnibus <- tryCatch(
            run_pathway_omnibus_only(df),
            error = function(e) {
              tibble(
                test = "Error",
                statistic = NA_real_,
                p = NA_real_,
                n_groups = dplyr::n_distinct(df$Group)
              )
            }
          )
          
          group_means <- tapply(
            df$Score,
            df$Group,
            mean,
            na.rm = TRUE
          )
          
          group_mean_variance <- if (length(group_means) >= 2) {
            stats::var(group_means, na.rm = TRUE)
          } else {
            NA_real_
          }
          
          score_range <- diff(range(df$Score, na.rm = TRUE))
          score_variance <- stats::var(df$Score, na.rm = TRUE)
          
          tibble(
            Pathway = pathway_name,
            test = omnibus$test[1],
            statistic = omnibus$statistic[1],
            pvalue = omnibus$p[1],
            n_groups = omnibus$n_groups[1],
            group_mean_variance = group_mean_variance,
            score_range = score_range,
            score_variance = score_variance
          )
        }
      ) %>%
        mutate(
          padj = adjust_pathway_pvalues(pvalue)
        ) %>%
        arrange(
          pvalue,
          desc(group_mean_variance),
          desc(score_range)
        )
      
      stat_df
    })
    
    filtered_pathway_stats <- reactive({
      filter_pathway_statistics(heatmap_pathway_stats(),
        source = input$heatmap_p_source %||% "adjusted",
        cutoff = input$heatmap_pvalue_cutoff %||% 0.05)
    })

    heatmap_matrix <- reactive({
      omics_number_parameter(input$top_n, "Number of displayed pathways", 1, .Machine$integer.max, TRUE)
      omics_number_parameter(input$heatmap_pvalue_cutoff %||% 0.05, "Heatmap P-value cutoff", 0, 1)
      
      req(
        ssgsea_scores(),
        sample_annotation()
      )
      
      es <- ssgsea_scores()
      top_n <- min(
        input$top_n,
        nrow(es)
      )
      
      selected_pathways <- filtered_pathway_stats() %>%
        slice_head(n = top_n) %>%
        pull(Pathway)
      
      validate(
        need(
          length(selected_pathways) > 0,
          "No pathways passed the selected P-value/FDR cutoff for the heatmap. Review the cutoff or filtering metric."
        )
      )
      
      mat <- es[selected_pathways, , drop = FALSE]
      mat <- mat[, sample_order(), drop = FALSE]
      
      if (isTRUE(input$scale_rows)) {
        
        mat <- t(scale(t(mat)))
        mat[!is.finite(mat)] <- 0
      }
      
      mat
    })
    
    # =====================================================
    # 10) Heatmap Plot
    # =====================================================
    
    create_heatmap_plot <- reactive({
      
      req(heatmap_matrix())
      
      mat <- heatmap_matrix()
      anno <- sample_annotation()
      
      plot_df <- as.data.frame(mat) %>%
        rownames_to_column("Pathway") %>%
        pivot_longer(
          cols = -Pathway,
          names_to = "Sample",
          values_to = "Score"
        )
      
      pathway_levels <- rownames(mat)
      
      plot_df$Pathway <- factor(
        plot_df$Pathway,
        levels = rev(pathway_levels)
      )
      
      plot_df$Sample <- factor(
        plot_df$Sample,
        levels = colnames(mat)
      )
      
      
      pathway_label_map <- stats::setNames(
        pathway_levels,
        pathway_levels
      )
      
      score_rng <- range(
        plot_df$Score,
        na.rm = TRUE
      )
      
      fill_scale <- if (score_rng[1] < 0 && score_rng[2] > 0) {
        
        scale_fill_gradient2(
          low = input$col_low,
          mid = input$col_mid,
          high = input$col_high,
          midpoint = 0
        )
        
      } else {
        
        scale_fill_gradientn(
          colours = c(
            input$col_low,
            input$col_mid,
            input$col_high
          )
        )
      }
      
      x_text_theme <- if (isTRUE(input$show_heatmap_sample_names)) {
        element_text(
          angle = 60,
          hjust = 1,
          vjust = 1,
          colour = "black"
        )
      } else {
        element_blank()
      }
      
      x_tick_theme <- if (isTRUE(input$show_heatmap_sample_names)) {
        element_line()
      } else {
        element_blank()
      }
      
      if (!is.null(anno)) {
        
        anno_levels <- anno %>%
          distinct(
            Sample,
            Group
          ) %>%
          mutate(
            Sample = factor(
              Sample,
              levels = colnames(mat)
            )
          )
        
        plot_df <- plot_df %>%
          left_join(
            anno_levels %>%
              mutate(Sample = as.character(Sample)),
            by = "Sample"
          )
        
        group_levels <- group_order_levels()
        
        if (length(group_levels) > 0 && "Group" %in% colnames(plot_df)) {
          plot_df <- plot_df %>%
            mutate(
              Group = factor(
                as.character(Group),
                levels = group_levels
              )
            )
        }
        
        ggplot(
          plot_df,
          aes(
            x = Sample,
            y = Pathway,
            fill = Score
          )
        ) +
          geom_tile() +
          facet_grid(
            . ~ Group,
            scales = "free_x",
            space = "free_x"
          ) +
          fill_scale +
          scale_y_discrete(
            labels = function(x) {
              labels <- pathway_label_map[as.character(x)]
              labels[is.na(labels)] <- as.character(x)[is.na(labels)]
              labels %>%
                stringr::str_replace_all("_", " ") %>%
                stringr::str_wrap(width = input$heatmap_wrap_width %||% 45)
            }
          ) +
          theme_bw(base_size = input$font_size) +
          theme(
            axis.text.x = x_text_theme,
            axis.text.y = element_text(size = input$font_size, color = "black"),
            axis.ticks.x = x_tick_theme,
            axis.title = element_blank(),
            plot.title = element_text(size = input$font_size + 2, face = "bold", color = "black", hjust = 0.5),
            panel.grid = element_blank(),
            strip.background = element_rect(fill = "grey95"),
            strip.text = element_text(size = input$font_size, face = "bold", color = "black"),
            legend.title = element_text(size = input$font_size, color = "black"),
            legend.text = element_text(size = input$font_size, color = "black")
          ) +
          labs(
            title = input$plot_title,
            fill = "Score"
          )
        
      } else {
        
        ggplot(
          plot_df,
          aes(
            x = Sample,
            y = Pathway,
            fill = Score
          )
        ) +
          geom_tile() +
          fill_scale +
          scale_y_discrete(
            labels = function(x) {
              labels <- pathway_label_map[as.character(x)]
              labels[is.na(labels)] <- as.character(x)[is.na(labels)]
              labels %>%
                stringr::str_replace_all("_", " ") %>%
                stringr::str_wrap(width = input$heatmap_wrap_width %||% 45)
            }
          ) +
          theme_bw(base_size = input$font_size) +
          theme(
            axis.text.x = x_text_theme,
            axis.text.y = element_text(size = input$font_size, color = "black"),
            axis.ticks.x = x_tick_theme,
            axis.title = element_blank(),
            plot.title = element_text(size = input$font_size + 2, face = "bold", color = "black", hjust = 0.5),
            panel.grid = element_blank(),
            legend.title = element_text(size = input$font_size, color = "black"),
            legend.text = element_text(size = input$font_size, color = "black")
          ) +
          labs(
            title = input$plot_title,
            fill = "Score"
          )
      }
    })
    
    output$ssgsea_heatmap <- renderPlot({
      
      req(create_heatmap_plot())
      
      print(create_heatmap_plot())
      
    }, res = 120)
    
    # =====================================================
    # 11) Pathway Selector
    # =====================================================
    
    observeEvent(ssgsea_scores(), {
      
      choices <- rownames(ssgsea_scores())
      
      updateSelectInput(
        session = session,
        inputId = "selected_pathway",
        choices = choices,
        selected = choices[1]
      )
    })
    
    output$pathway_selector_ui <- renderUI({
      
      req(ssgsea_scores())
      
      selectInput(
        ns("selected_pathway"),
        "Pathway:",
        choices = rownames(ssgsea_scores()),
        selected = rownames(ssgsea_scores())[1]
      )
    })
    
    # =====================================================
    # 12) Selected Pathway Data and Statistics
    # =====================================================
    
    selected_pathway_df <- reactive({
      
      req(
        ssgsea_scores(),
        input$selected_pathway
      )
      
      es <- ssgsea_scores()
      anno <- sample_annotation()
      
      validate(
        need(
          input$selected_pathway %in% rownames(es),
          "Selected pathway not found."
        )
      )
      
      df <- tibble(
        Sample = colnames(es),
        Score = as.numeric(es[input$selected_pathway, ])
      )
      
      if (!is.null(anno)) {
        
        df <- df %>%
          left_join(
            anno %>%
              mutate(Sample = as.character(Sample)),
            by = "Sample"
          )
        
        group_levels <- group_order_levels()
        
        if (length(group_levels) == 0) {
          group_levels <- unique(as.character(df$Group))
        }
        
        df <- df %>%
          mutate(
            Group = factor(
              as.character(Group),
              levels = group_levels
            )
          ) %>%
          arrange(
            Group,
            Sample
          )
      }
      
      df
    })
    
    selected_pathway_stats <- reactive({
      
      req(selected_pathway_df())
      
      df <- selected_pathway_df()
      
      validate(
        need(
          "Group" %in% colnames(df),
          "Sample annotation is required for group comparison."
        )
      )
      
      result <- run_pathway_group_statistics(
        df = df,
        p_adjust_method = input$posthoc_adjust_method %||% "BH",
        show_ns = FALSE
      )
      if (result$n_groups < 2 || nrow(result$omnibus_table) == 0) return(result)

      global <- selected_pathway_global_stats()
      validate(need(nrow(global) == 1, "Global pathway statistics are unavailable."))
      n_tested <- sum(is.finite(heatmap_pathway_stats()$pvalue))
      result$omnibus_table <- result$omnibus_table %>%
        mutate(p.adj = global$padj[1], adjustment = "BH across all tested pathways",
          n_pathways_tested = n_tested)
      result$omnibus_label <- paste0(global$test[1], ": ",
        format_p_label(global$pvalue[1], "Nominal P"), " | ",
        format_p_label(global$padj[1], "FDR (BH)"))

      # With two groups there is one contrast per pathway: its discovery FDR
      # belongs to the family of all pathway tests, not a one-P-value adjustment.
      if (result$n_groups == 2 && nrow(result$pairwise_table) > 0) {
        result$pairwise_table <- result$pairwise_table %>%
          mutate(p.adj = global$padj[1], p_plot = p.adj,
            p.signif = p_to_stars(p.adj), p.label = format_p_label(p.adj, "FDR (BH)"),
            adjustment = "BH across all tested pathways")
        result$plot_pairwise_table <- add_y_positions(
          result$pairwise_table %>% filter(!is.na(p_plot), p_plot <= 0.05), df)
        scope <- "Bracket significance: BH FDR across all tested pathways"
      } else {
        scope <- paste0("Pairwise significance: adjustment within this pathway (",
          paste(unique(result$pairwise_table$test), collapse = "; "),
          "). Global FDR: BH across all tested pathways")
      }
      result$assumption_caption <- paste(result$assumption_caption, scope, sep = "\n")
      result
    })
    
    
    selected_pathway_global_stats <- reactive({
      
      req(
        ssgsea_scores(),
        input$selected_pathway,
        sample_annotation()
      )
      
      if (dplyr::n_distinct(sample_annotation()$Group) < 2) return(tibble())
      stats_df <- heatmap_pathway_stats()
      
      stats_df %>%
        filter(Pathway == input$selected_pathway) %>%
        slice_head(n = 1)
    })
    
    # =====================================================
    # 13) Selected Pathway Boxplot / Profile
    # =====================================================
    
    create_pathway_profile_plot <- reactive({
      
      req(selected_pathway_df())
      
      df <- selected_pathway_df()
      
      if ("Group" %in% colnames(df)) {
        
        stats_res <- selected_pathway_stats()
        
        y_rng <- range(df$Score, na.rm = TRUE)
        y_span <- diff(y_rng)
        
        if (!is.finite(y_span) || y_span == 0) {
          y_span <- 1
        }
        
        pairwise_plot <- stats_res$plot_pairwise_table
        
        y_upper <- if (!is.null(pairwise_plot) && nrow(pairwise_plot) > 0) {
          max(pairwise_plot$y.position, na.rm = TRUE) + 0.12 * y_span
        } else {
          y_rng[2] + 0.18 * y_span
        }
        
        detected_group_levels <- detected_groups()
        group_levels <- group_order_levels()
        
        if (length(group_levels) == 0) {
          group_levels <- unique(as.character(df$Group))
        }
        
        group_levels <- group_levels[group_levels %in% unique(as.character(df$Group))]
        
        if (length(group_levels) == 0) {
          group_levels <- unique(as.character(df$Group))
        }
        
        df <- df %>%
          mutate(
            Group = factor(
              as.character(Group),
              levels = group_levels
            )
          )
        
        fill_map <- get_group_color_values(
          input = input,
          detected_groups = detected_group_levels,
          prefix = "box_fill_",
          default_palette = get_default_box_fill_palette()
        )
        
        point_map <- get_group_color_values(
          input = input,
          detected_groups = detected_group_levels,
          prefix = "point_col_",
          default_palette = get_default_point_palette()
        )
        
        fill_vals <- fill_map[group_levels]
        point_vals <- point_map[group_levels]
        
        names(fill_vals) <- group_levels
        names(point_vals) <- group_levels
        
        global_stats <- selected_pathway_global_stats()
        
        subtitle_text <- stats_res$omnibus_label
        
        if (!is.null(global_stats) && nrow(global_stats) > 0) {
          subtitle_text <- paste0(global_stats$test[1], ": ",
            format_p_label(global_stats$pvalue[1], "Nominal P"), " | ",
            format_p_label(global_stats$padj[1], "FDR (BH)"))
        }
        
        p <- ggplot(
          df,
          aes(
            x = Group,
            y = Score
          )
        ) +
          geom_boxplot(
            aes(fill = Group),
            width = 0.58,
            outlier.shape = NA,
            alpha = 0.90,
            color = "grey30",
            linewidth = input$box_linewidth
          ) +
          geom_jitter(
            aes(
              color = Group,
              fill = Group
            ),
            width = 0.10,
            height = 0,
            size = 2.25,
            alpha = 0.90,
            shape = 21,
            stroke = input$point_stroke
          ) +
          scale_fill_manual(
            values = fill_vals
          ) +
          scale_color_manual(
            values = point_vals
          ) +
          coord_cartesian(
            ylim = c(y_rng[1] - 0.06 * y_span, y_upper),
            clip = "off"
          ) +
          theme_bw(base_size = input$font_size) +
          theme(
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(
              color = "grey92",
              linewidth = 0.35
            ),
            axis.text.x = element_text(
              angle = 30,
              hjust = 1,
              size = input$font_size,
              color = "black"
            ),
            axis.text.y = element_text(
              size = input$font_size,
              color = "black"
            ),
            axis.title = element_text(
              size = input$font_size,
              color = "black"
            ),
            plot.title = element_text(
              size = input$font_size + 2,
              face = "bold",
              color = "black",
              hjust = 0.5
            ),
            plot.subtitle = element_text(
              size = input$font_size,
              color = "black",
              hjust = 0.5
            ),
            plot.caption = element_text(
              hjust = 0,
              size = input$font_size,
              color = "black"
            ),
            legend.position = "none",
            panel.border = element_rect(
              color = "black",
              linewidth = 0.35
            ),
            plot.margin = margin(
              t = 10,
              r = 28,
              b = 10,
              l = 10
            )
          ) +
          labs(
            title = input$selected_pathway,
            subtitle = subtitle_text,
            x = NULL,
            y = "ssGSEA score",
            caption = stats_res$assumption_caption
          )
        
        if (!is.null(pairwise_plot) && nrow(pairwise_plot) > 0) {
          
          p <- p +
            ggpubr::stat_pvalue_manual(
              pairwise_plot,
              label = "p.signif",
              xmin = "group1",
              xmax = "group2",
              y.position = "y.position",
              tip.length = 0,
              bracket.size = input$bracket_linewidth,
              size = 3.8
            )
        }
        
        p
        
      } else {
        
        df <- df %>%
          mutate(
            Sample = factor(
              Sample,
              levels = colnames(ssgsea_scores())
            )
          )
        
        ggplot(
          df,
          aes(
            x = Sample,
            y = Score,
            group = 1
          )
        ) +
          geom_line(
            color = input$profile_col,
            linewidth = 0.8
          ) +
          geom_point(
            color = input$profile_col,
            size = 2.5
          ) +
          theme_bw(base_size = input$font_size) +
          theme(
            axis.text.x = element_text(
              angle = 60,
              hjust = 1,
              vjust = 1,
              size = input$font_size,
              color = "black"
            ),
            axis.text.y = element_text(
              size = input$font_size,
              color = "black"
            ),
            axis.title = element_text(
              size = input$font_size,
              color = "black"
            ),
            plot.title = element_text(
              size = input$font_size + 2,
              face = "bold",
              color = "black",
              hjust = 0.5
            )
          ) +
          labs(
            title = input$selected_pathway,
            x = NULL,
            y = "ssGSEA score"
          )
      }
    })
    
    output$pathway_profile <- renderPlot({
      
      req(create_pathway_profile_plot())
      
      print(create_pathway_profile_plot())
      
    }, res = 120)
    
    # =====================================================
    # 14) Selected Pathway Statistical Outputs
    # =====================================================
    
    output$pathway_stats_summary <- renderPrint({
      
      req(selected_pathway_df())
      
      df <- selected_pathway_df()
      
      if (!"Group" %in% colnames(df)) {
        cat("Sample annotation is required for group comparison.\n")
        return()
      }
      
      stats_res <- selected_pathway_stats()
      
      cat("Selected pathway:", input$selected_pathway, "\n\n")
      cat("Number of groups:", stats_res$n_groups, "\n")
      cat("Omnibus test:\n")
      print(stats_res$omnibus_table)
      cat("\nGlobal adjustment: BH across all tested pathways, before heatmap filtering or Top N selection.\n")
      cat("\nPairwise comparisons (adjustment scope shown explicitly):\n")
      print(stats_res$pairwise_table)
      
      cat("\nNormality by group:\n")
      print(stats_res$normality_table)
      
      cat("\nHomogeneity of variances:\n")
      print(stats_res$variance_table)
      
      cat("\nDecision used in the plot:\n")
      cat(stats_res$omnibus_label, "\n")
      cat(stats_res$assumption_caption, "\n")
    })
    
    # =====================================================
    # 15) Tables
    # =====================================================
    
    score_table_df <- reactive({
      
      req(ssgsea_scores())
      
      as.data.frame(ssgsea_scores()) %>%
        rownames_to_column("Pathway")
    })
    
    processed_matrix_df <- reactive({
      
      req(expr_matrix())
      
      as.data.frame(expr_matrix()) %>%
        rownames_to_column("Feature")
    })
    
    pathway_statistics_table <- reactive({
      
      req(heatmap_pathway_stats())
      
      source <- input$heatmap_p_source %||% "adjusted"
      metric <- if (source == "adjusted") "padj" else "pvalue"
      filtered_pathway_stats() %>%
        transmute(
          Pathway = Pathway,
          Test = test,
          Statistic = statistic,
          `Nominal P-value` = pvalue,
          `FDR (BH across pathways)` = padj,
          `Heatmap filter metric` = source,
          `Heatmap cutoff` = input$heatmap_pvalue_cutoff %||% 0.05,
          `Pathways tested for BH` = sum(is.finite(heatmap_pathway_stats()$pvalue)),
          Significance = p_to_stars(.data[[metric]])
        )
    })
    
    output$ssgsea_table <- renderDT({
      
      req(score_table_df())
      
      datatable(
        score_table_df(),
        options = list(
          pageLength = 10,
          scrollX = TRUE
        ),
        rownames = FALSE
      )
    })
    
    # =====================================================
    # 16) Downloads
    # =====================================================
    
    output$download_table <- downloadHandler(
      
      filename = function() {
        paste0(
          "ssgsea_scores_",
          Sys.Date(),
          ".txt"
        )
      },
      
      content = function(file) {
        
        write.table(
          score_table_df(),
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
          "ssgsea_processed_expression_matrix_",
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
    
    output$download_pathway_statistics <- downloadHandler(
      
      filename = function() {
        paste0(
          "ssgsea_significant_pathway_pvalues_",
          Sys.Date(),
          ".txt"
        )
      },
      
      content = function(file) {
        
        write.table(
          pathway_statistics_table(),
          file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )
      }
    )
    
    output$download_heatmap <- downloadHandler(
      
      filename = function() {
        paste0(
          "ssgsea_heatmap_",
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
        
        open_ssgsea_plot_device(
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
    
    output$download_boxplot <- downloadHandler(

      filename = function() {
        paste0(
          "ssgsea_selected_pathway_boxplot_",
          Sys.Date(),
          "_",
          input$exp_width_box, "x", input$exp_height_box,
          ".",
          input$exp_fmt
        )
      },

      content = function(file) {

        width_in <- input$exp_width_box * 0.393701
        height_in <- input$exp_height_box * 0.393701

        open_ssgsea_plot_device(
          file = file,
          fmt = input$exp_fmt,
          width = width_in,
          height = height_in,
          dpi = input$exp_dpi
        )

        print(create_pathway_profile_plot())
        dev.off()
      }
    )
    
    # =====================================================
    # 17) Return Reactive Objects
    # =====================================================
    
    return(
      list(
        expr_matrix = expr_matrix,
        sample_annotation = sample_annotation,
        gene_sets = gene_sets,
        ssgsea_scores = ssgsea_scores,
        selected_pathway_df = selected_pathway_df,
        selected_pathway_stats = selected_pathway_stats,
        selected_pathway_global_stats = selected_pathway_global_stats,
        heatmap_pathway_stats = heatmap_pathway_stats
      )
    )
  })
}
