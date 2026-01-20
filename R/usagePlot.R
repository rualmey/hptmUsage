#' Generate Usage Plots
#'
#' @description
#' This function generates line plots showing the usage of histone variants,
#' precursors, or PTMs across samples. It visualizes the relevant usage
#' normalization factors, corrected usages, and model estimates.
#'
#' @param object A `QFeatures` object containing the histone data, including
#'   models from [msqrob2::msqrob()] etc.
#' @param assays A named `list()` of assays (`character()` or `integer()`) to
#'   plot (hPTM and/or variant group level), where names should correspond to
#'   the assay level, i.e., "ptm" or "variant". For example `list(ptm = "ptm",
#'   variant = "variant"),`.
#' @param model A named `character()` of formula(e) that were fit using
#'   [msqrob2::msqrob()]. Used to reconstruct model estimates. For example,
#'   `c(meansReference = "~ group", regressionTreat = "~ cell_line + treatment")`.
#' @param significant_only Logical. Whether to generate lineplots for all
#'   features or only those that weres significant in **any** contrasts of the
#'   supplied model(s).
#' @param show_legend Logical. Whether to show the legend.
#' @returns A list of `ggplot` objects, nested by assay level and assay name.
#' @export
usagePlot <- function(
  object,
  assays,
  model,
  significant_only = TRUE,
  show_legend = TRUE
) {
  stopifnot("'object' must be a QFeatures object." = is(object, "QFeatures"))
  stopifnot("'assays' must only contain names 'ptm' and/or 'variant'" = all(names(assays) %in% c("ptm", "variant")))
  purrr::imap(assays, \(x, idx) list(QFeatures:::.normIndex(object, x)) |> purrr::set_names(idx))

  # retrieve feature data = summarized and estimated PTM or variant usage
  features <- purrr::imap(assays, \(x, idx) .get_features(x, idx, object, model, significant_only))
  summarized_usage <- purrr::imap(features, \(x, idx) .get_summarized_usage(x, idx, object))
  estimated_usage <- purrr::imap(features, \(x, idx) .get_estimated_usage(x, idx, object, model))

  # retrieve precursor data = abundance and usage of corresponding precursors
  precursor_features <- purrr::imap(features, \(x, idx) .get_precursor_features(x, idx, object))
  precursor_abundance <- purrr::imap(precursor_features, \(x, idx) {
    .get_precursor_abundance(x, idx, object, usage = FALSE)
  })
  precursor_usage <- purrr::imap(precursor_features, \(x, idx) .get_precursor_abundance(x, idx, object, usage = TRUE))

  # retrieve normalization factors
  norm_factors <- purrr::imap(features, \(x, idx) .get_norm_factors(x, idx, object))

  # plot
  plot_list <- purrr::map_depth(features, 2, \(x) vector("list", length(x)) |> purrr::set_names(x))
  for (assay_level in seq_along(plot_list)) {
    for (assay in seq_along(plot_list[[assay_level]])) {
      for (i in seq_along(plot_list[[assay_level]][[assay]])) {
        # hPTM/variant features
        feat <- features[[assay_level]][[assay]][[i]]
        summ <- summarized_usage[[assay_level]][[assay]] |> dplyr::filter(.data[["feature"]] == feat)
        est <- estimated_usage[[assay_level]][[assay]] |> dplyr::filter(.data[["feature"]] == feat)

        # precursors
        prec_ab <- precursor_abundance[[assay_level]][[assay]][[i]]
        prec_us <- precursor_usage[[assay_level]][[assay]][[i]]

        # normalization factors
        usage_level <- norm_factors[[assay_level]][[assay]][["usage_level"]]
        fact <- norm_factors[[assay_level]][[assay]][["scaling_factors"]][i, , drop = FALSE] |>
          tibble::as_tibble() |>
          tidyr::pivot_longer(tidyr::everything(), names_to = "sample", values_to = "value") |>
          dplyr::mutate(
            feature = "normalization_factor",
            level = dplyr::case_when(
              usage_level == "histone" ~ "Histone Abundance",
              usage_level == "histone_family" ~ "Histone Family Abundance",
              usage_level == "histone_group" ~ "Variant Group Abundance",
            )
          ) |>
          dplyr::relocate(feature)

        # break range for ggbreak
        break_min <- dplyr::bind_rows(summ, est, prec_us) |> tidyr::drop_na() |> dplyr::pull(value) |> max()
        break_max <- dplyr::bind_rows(prec_ab, fact) |> tidyr::drop_na() |> dplyr::pull(value) |> min()

        plot_data <- dplyr::bind_rows(summ, est, prec_ab, prec_us, fact) |>
          .format_for_plotting(object)

        # Generate the plot
        plot_list[[assay_level]][[assay]][[i]] <- .plot_usage(
          plot_data,
          title = feat,
          show_legend = show_legend,
          breaks = c(break_min, break_max),
          groups = colData(object)[
            colData(object) |> with(order(colData(object)$group, colData(object)$quantCols))
          ]$group |>
            table() |>
            as.integer() +
            0.5
        )
      }
    }
  }

  return(plot_list)
}

# ------------------------------------------------------------------------------
# Internal Helper Functions
# ------------------------------------------------------------------------------

.get_features <- function(assay_names, assay_level, object, model, significant_only) {
  purrr::map(assay_names, function(x) {
    rd <- rowData(object[[x]])
    # rowData columns with data.frames containing adjPval
    contrast_df <- paste0(names(model), "/") |>
      paste(collapse = "|") |>
      stringr::str_starts(names(rd), pattern = _)
    contrast_df <- purrr::map(names(rd)[contrast_df], \(x) tibble::as_tibble(rd[[x]], rownames = "name")) |>
      dplyr::bind_rows()
    if (significant_only) {
      contrast_df <- contrast_df |>
        dplyr::filter(.data[["adjPval"]] <= 0.05)
    }
    contrast_df <- contrast_df |>
      dplyr::arrange(.data[["adjPval"]]) |>
      dplyr::pull(.data[["name"]]) |>
      unique()
  }) |>
    # store assay name for later retrieval
    purrr::set_names(assay_names)
}

.get_summarized_usage <- function(features, assay_level, object) {
  purrr::imap(features, function(x, idx) {
    SummarizedExperiment::assay(object[[idx]][x]) |>
      tibble::as_tibble(rownames = "feature") |>
      tidyr::pivot_longer(!feature, names_to = "sample", values_to = "value") |>
      dplyr::mutate(level = if (assay_level == "ptm") "PTM Usage" else "Variant Group Usage")
  })
}

.get_estimated_usage <- function(features, assay_level, object, model) {
  purrr::imap(features, function(x, idx) {
    purrr::imap(model, function(model_formula, model_name) {
      design_mat <- model.matrix(as.formula(model_formula), data = droplevels(colData(object)))
      purrr::map(rowData(object[[idx]][x])[[model_name]], function(model) {
        coef <- msqrob2::getCoef(model)
        # fix ridge regression as this prepends "ridge" to the factor name
        names(coef) <- stringr::str_remove(names(coef), "^ridge")
        fixed_coef <- intersect(names(coef), colnames(design_mat))
        # only calculate fixed effects as random effects are of no interest
        design_mat[, fixed_coef, drop = FALSE] %*%
          coef[fixed_coef] |>
          as.numeric() |>
          purrr::set_names(rownames(colData(object)))
      }) |>
        tibble::as_tibble() |>
        dplyr::mutate(sample = rownames(colData(object))) |>
        tidyr::pivot_longer(!sample, names_to = "feature", values_to = "value") |>
        dplyr::relocate(feature) |>
        dplyr::mutate(
          level = if (assay_level == "ptm") {
            paste0("PTM Estimated (", model_name, ")")
          } else {
            paste0("Variant Group Estimated (", model_name, ")")
          }
        )
    }) |>
      dplyr::bind_rows()
  })
}

.link_data <- function(object, assay, assay_level = NULL, usage = TRUE) {
  links <- QFeatures::assayLinks(object, assay)
  link_fcol <- links[[1]]@fcol
  # PTM assay has a "deconvoluted" precursor assay in between
  if (is.null(assay_level) || assay_level == "variant") {
    precursor_assay <- if (isTRUE(usage)) names(links)[[2]] else names(links)[[3]]
  } else if (assay_level == "ptm") {
    precursor_assay <- if (isTRUE(usage)) names(links)[[3]] else names(links)[[4]]
  }

  list(precursor_assay = precursor_assay, fcol = link_fcol)
}

.get_precursor_features <- function(assay_features, assay_level, object) {
  purrr::imap(assay_features, function(x, idx) {
    link_data <- .link_data(object, idx)
    rd <- rowData(object[[link_data$precursor_assay]])
    purrr::map(x, \(y) rd[rd[[link_data$fcol]] == y, ]$precursor) |>
      purrr::set_names(x)
  })
}

.get_precursor_abundance <- function(features, assay_level, object, usage) {
  purrr::imap(features, function(x, idx) {
    link_data <- .link_data(object, idx, assay_level = assay_level, usage = usage)
    precursor_assay <- SummarizedExperiment::assay(object[[link_data$precursor_assay]])
    purrr::map(x, function(y) {
      precursor_assay[y, , drop = FALSE] |>
        tibble::as_tibble(rownames = "feature") |>
        tidyr::pivot_longer(!feature, names_to = "sample", values_to = "value") |>
        dplyr::mutate(level = if (isTRUE(usage)) "Precursor Usage" else "Precursor Abundance")
    })
  })
}

.get_norm_factors <- function(assay_features, assay_level, object) {
  purrr::imap(assay_features, function(x, idx) {
    # retrieve metadata
    precursor_assay <- .link_data(object, idx, assay_level)$precursor_assay
    metadata <- metadata(object[[precursor_assay]])

    # match factors to features
    factor_names <- rowData(object[[idx]][x])[[metadata$usage_level]]
    if (metadata$usage_level == "histone") {
      factor_names <- as.integer(factor_names)
    }
    scaling_factors <- metadata$scaling_factors[factor_names, ]

    list(
      scaling_factors = scaling_factors,
      usage_level = metadata$usage_level
    )
  })
}

.format_for_plotting <- function(df, object) {
  df |>
    dplyr::mutate(
      # grouping for individual lines
      group_id = paste(level, feature),
      # increase point size and linewidth for important assays
      lw = dplyr::case_when(
        level %in% c("Histone Abundance", "Histone Family Abundance", "Variant Group Abundance") ~ 1,
        stringr::str_detect(level, "Estimated") ~ 1,
        .default = .5
      ),
      point = dplyr::case_when(
        stringr::str_detect(level, "Estimated") ~ NA_real_,
        .default = value
      ),
      # define shape grouping for unique precursor points
      shape_group = dplyr::case_when(
        stringr::str_detect(level, "Precursor") ~ stringr::str_wrap(feature, 30, whitespace_only = FALSE),
        .default = NA_character_
      ),
      # reorder levels for legend, estimated is kept as last so it always overlaps
      level = suppressWarnings(forcats::fct_relevel(
        level,
        "Precursor Abundance",
        "Precursor Usage",
        "Histone Abundance",
        "Histone Family Abundance",
        "Variant Group Abundance",
        "Variant Group Usage",
        "PTM Usage"
      )),
      # run order should go by group then by sample name
      sample = factor(
        sample,
        rownames(colData(object))[colData(object) |> with(order(colData(object)$group, colData(object)$quantCols))]
      )
    )
}

.plot_usage <- function(df, title, show_legend, breaks = NULL, groups = NULL) {
  color_map <- c(
    "Histone Abundance" = "#5e81ac",
    "Histone Family Abundance" = "#5e81ac",
    "Variant Group Abundance" = "#5e81ac",
    "Precursor Abundance" = "#3b4252",
    "Precursor Usage" = "#4c566a",
    "Variant Group Usage" = "#bf616a",
    "PTM Usage" = "#bf616a"
  )
  alpha_map <- c(
    "Histone Abundance" = 1,
    "Histone Family Abundance" = 1,
    "Variant Group Abundance" = 1,
    "Variant Group Usage" = 1,
    "Precursor Abundance" = .5,
    "Precursor Usage" = .3,
    "PTM Usage" = 1
  )
  estimated_names <- df |>
    dplyr::filter(stringr::str_detect(level, "Estimated")) |>
    dplyr::pull(level) |>
    unique()
  color_map <- c(
    color_map,
    c("#b48ead", "#a3be8c", "#d08770", "#ebcb8b", "#8fbcbb", "#88c0d0")[seq_along(estimated_names)] |>
      purrr::set_names(as.character(estimated_names))
  )
  alpha_map <- c(alpha_map, rep.int(1, length(estimated_names)) |> purrr::set_names(as.character(estimated_names)))

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = sample,
      y = value,
      group = group_id,
      color = level,
      alpha = level
    )
  ) +

    ggplot2::geom_point(ggplot2::aes(y = point, shape = shape_group), size = 1, show.legend = TRUE) +
    ggplot2::geom_line(ggplot2::aes(linewidth = lw)) +
    {
      if (!is.null(groups)) ggplot2::geom_vline(xintercept = groups, linetype = "dashed")
    } +
    {
      if (!all(is.null(breaks)) && breaks[[1]] < breaks[[2]]) {
        ggbreak::scale_y_break(c(breaks[[1]], breaks[[2]]), scales = 1 / 1.618)
      }
    } +

    ggplot2::scale_colour_manual(values = color_map, na.value = "grey50") +
    ggplot2::scale_alpha_manual(values = alpha_map, na.value = 1) +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_shape_manual(values = c(16:25, 0:15, 97:122, 65:90), na.translate = FALSE) +

    ggplot2::labs(
      title = stringr::str_wrap(title, 60, whitespace_only = FALSE),
      x = "Run",
      y = NULL,
      colour = "Assay",
      alpha = "Assay",
      shape = "Precursor"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
      legend.position = if (show_legend) "right" else "none"
    )

  return(p)
}
