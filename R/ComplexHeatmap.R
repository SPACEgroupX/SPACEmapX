# To Do List:
# - change coloring variable names to match options

#' Plot Complex Heatmap for inferCNV Output
#'
#' This function generates a heatmap using the ComplexHeatmap package to visualize inferCNV results,
#' supporting reference and observation groups, one additional annotation, and user-defined color schemes.
#'
#' Barcodes needs to be formatted in the following way: pt03_AP1_CAATGGATACGCTCGA.1 patient_section_10xBarcodeID.1 with Patient ID being 4 characters long.
#'
#' @param infercnv_obj An inferCNV object created by the `infercnv::create_infercnv_object` and processed by `run()` or `run_infercnv()`.
#' @param output_dir Directory where the output heatmap PDF will be saved.
#' @param heatmap_thresholds_file Path to the inferCNV-generated `heatmap_thresholds.txt` file.
#' @param dendrogram_file A Newick file for dendrogram, can be either a single tree (`phylo`) or multiple trees (`multiPhylo`).
#' @param obs_class Character string indicating how observations are grouped. Example `"Clone"`, `"Histology"`.
#' @param annotation_name Character string indicating what additional annotation is desired. Example `"Histology"`.
#' @param annotation_df A data frame or path to a TSV file with additional annotation information (annotation should be in the first column). Required if annotation_type is used.
#' @param output_format File format for the saved heatmap: `"png"` (default) or `"pdf"`.
#' @param output_name File name for the saved heatmap. Default is `"complex_heatmap.png"`.
#' @param class_cols Optional. Named list or vector of colors for classes groups within object. Format: `c("A" = "#hex", ...)`.
#' @param anno_cols Optional. Named list or vector of colors for additional annotation groups. Format: `c("A" = "#hex", ...)`.
#' @param ref_pie_chart Optional. Adds a section distribution pie chart of the reference
#'
#' @return A image or pdf file is saved to `output_dir`, visualizing CNV patterns across references and observations with dendrograms and annotations.
#'
#' @details
#' - For `dendrogram_file`, a single tree (`phylo`) will apply a global dendrogram, while a multiple tree (`multiPhylo`) input will order spots by original grouping.
#' - The function will auto-detect and adjust for reference presence.
#' - Colors for sample, histology, and clone groups are automatic but can be overridden.
#' - Chromosome and CNV coloring is derived from inferCNV settings and thresholds.
#'
#' @examples
#' plot_complex_heatmap(
#'   infercnv_obj = infercnv_final,
#'   output_dir = "results/",
#'   heatmap_thresholds_file = "results/infercnv.heatmap_thresholds.txt",
#'   dendrogram_file = "results/infercnv.observations_dendrogram.txt",
#'   obs_class = "clone",
#'   annotation_name = "histology",
#'   annotation_df = "metadata/histologies.tsv",
#'   class_cols = c("A" = "#FF0000", "B" = "#0000FF"),
#'   anno_cols = c("A" = "#E41A1C", "B" = "#377EB8"),
#'   ref_pie_chart = TRUE
#' )

plot_complex_heatmap <- function(infercnv_obj,
                                 output_dir,
                                 heatmap_thresholds_file,
                                 dendrogram_file,
                                 obs_class = "class",
                                 annotation_name = NULL,
                                 annotation_df = NULL,
                                 output_format = "png",
                                 output_name = "complex_heatmap.png",
                                 class_cols = NULL,
                                 anno_cols = NULL,
                                 ref_pie_chart = FALSE) {

  library(ComplexHeatmap)
  library(tidyverse)
  library(RColorBrewer)
  library(circlize)
  library(grid)
  library(gridExtra)
  library(ape)
  library(ggplot2)
  library(ggrepel)
  library(dendextend)

  if (!is.null(annotation_name) && is.null(annotation_df)) {
    stop("`annotation_df` must be supplied if `annotation_name` is used")
  }

  if (!output_format %in% c("png", "pdf")) {
    stop("`output_format` must be either 'png' or 'pdf'.")
  }

  # Fetch expressions data
  expr <- infercnv_obj@expr.data

  # Gene ordering
  gn <- rownames(expr)
  geneOrder <- infercnv_obj@gene_order
  geneOrder <- geneOrder %>% filter(!str_detect(chr, "chrM"))
  sub_geneOrder <- geneOrder[intersect(gn, rownames(geneOrder)), ]
  expr <- expr[intersect(gn, rownames(geneOrder)), ]

  # Colors
  ## Sample
  unique_samples <- unique(str_remove(str_remove(colnames(expr), "_[^_]*$"), "^[^_]*_"))
  if (length(unique_samples) <= 9) {
    sample_colors <- setNames(RColorBrewer::brewer.pal(length(unique_samples), "Set3"), unique_samples)
  } else {
    sample_colors <- setNames(
      colorRampPalette(RColorBrewer::brewer.pal(9, "Set3"))(length(unique_samples)), unique_samples)
  }
  sample_colors <- sample_colors[!is.na(names(sample_colors))]

  ## Chromosome
  chr_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique(geneOrder$chr)))

  ## CNV
  cnv_colors <- c("#00008B", "#24249B", "#4848AB", "#6D6DBC", "#9191CC", "#B6B6DD", "#DADAEE", "#FFFFFF",
                  "#EEDADA", "#DDB6B6", "#CC9191", "#BC6D6D", "#AB4848", "#9B2424", "#8B0000")

  ## Histology
  hist_colors_default <- c("#0077B6", "#F4C0BA", "#F7952B", "#E54C2E", "#DD0303", "#7C2323", "#FFEA00", "#29BF12", "#212529", "#6C757D", "#AFE8F4")

  ## Clone
  # Extra colors "#CCCCCC", "#736F72", "#000000"
  #clone_colors_full <- c("#CAF0F8", "#80FFDB", "#00B4D8", "#01497C", "#155D27", "#29BF12", "#ABFF4F", "#FCF300", "#FF8000", "#E54C2E", "#D40000", "#780000", "#9B3F11", "#A24CCD", "#5A189A")
  clone_colors_default <- c("#CAF0F8", "#80FFDB", "#00B4D8", "#01497C", "#155D27", "#29BF12", "#ABFF4F", "#FCF300", "#FF8000", "#E54C2E")
  clone_x_colors_default <- c("#D40000", "#780000", "#9B3F11", "#A24CCD", "#5A189A")

  cnv_breakpoints <- read_csv(heatmap_thresholds_file, col_names = FALSE, show_col_types = FALSE)[[1]]
  cnv_col_fun <- function(x) {
    idx <- findInterval(x, cnv_breakpoints)
    idx[idx < 1] <- 1
    idx[idx > length(cnv_colors)] <- length(cnv_colors)
    cnv_colors[idx]
  }

  # Annotations - Observations
  obs_loc <- infercnv_obj@observation_grouped_cell_indices
  names(obs_loc) <- make.names(names(obs_loc))

  suppressWarnings({
    list2env(obs_loc, envir = environment())
  })

  display_name_map_obs <- names(infercnv_obj@observation_grouped_cell_indices)
  names(display_name_map_obs) <- names(obs_loc)

  obs_anno_list <- list()
  for (annotation in names(obs_loc)) {
    obs_anno_list[[annotation]] <- data.frame(
      CB = colnames(expr)[get(annotation)],
      class = display_name_map_obs[annotation])
  }
  obs_anno_df <- do.call(rbind, obs_anno_list)

  obs_expr <- expr[, colnames(expr) %in% obs_anno_df$CB]
  obs_anno_df <- obs_anno_df[match(colnames(obs_expr), obs_anno_df$CB), ]

  rownames(obs_anno_df) <- obs_anno_df$CB
  obs_anno_df$CB <- NULL

  obs_anno_df$Section <- str_remove(str_remove(rownames(obs_anno_df), "_[^_]*$"), "^[^_]*_")

  # Calculate the heatmaps plot heights
  obs_plot_height = ncol(obs_expr) / ncol(expr) * 12
  ref_plot_height = (ncol(expr) - ncol(obs_expr)) / ncol(expr) * 12

  # Dendrogram
  tree <- read.tree(dendrogram_file)
  if (inherits(tree, "phylo")) {

    dend <- as.dendrogram(as.hclust(tree))
    dend <- rev(dend)

    # Re-order based on dendrogram
    dend_labels <- labels(dend)
    missing_labels <- setdiff(dend_labels, colnames(obs_expr))
    if (length(missing_labels) > 0) {
      stop("Some labels in the dendrogram are not found in obs_expr: ", paste(missing_labels, collapse = ", "))
    }
    obs_expr <- obs_expr[, labels(dend)]
    obs_anno_df <- obs_anno_df[labels(dend), ]

    row_split <- NULL
  } else if (inherits(tree, "multiPhylo")) {

    clone_dend_list <- lapply(tree, function(single_tree) {
      as.dendrogram(as.hclust(single_tree))
    })

    clone_order_list <- lapply(clone_dend_list, labels)
    clone_spot_order <- unlist(clone_order_list)
    names(clone_order_list) <- names(infercnv_obj@observation_grouped_cell_indices)

    clone_names <- names(clone_order_list)
    clone_lengths <- sapply(clone_order_list, length)
    row_split <- factor(rep(clone_names, times = clone_lengths), levels = clone_names)

    # Re-order based on dendrogram
    missing_labels <- setdiff(clone_spot_order, colnames(obs_expr))
    if (length(missing_labels) > 0) {
      stop("Some labels in the dendrogram are not found in obs_expr: ", paste(missing_labels, collapse = ", "))
    }
    obs_expr <- obs_expr[, clone_spot_order]
    obs_anno_df <- obs_anno_df[clone_spot_order, ]

    dend <- function(mat) hclust(dist(mat), method = "ward.D2")
  }

  # Set ht_options
  ht_opt$ROW_ANNO_PADDING <- unit(3, "mm")
  ht_opt$DENDROGRAM_PADDING <- unit(3, "mm")
  ht_opt$HEATMAP_LEGEND_PADDING <- unit(5, "mm")
  ht_opt$ANNOTATION_LEGEND_PADDING <- unit(10, "mm")
  ht_opt$message <- FALSE

  # Obs Heatmap object
  top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = chr_colors, col = NA),
                                                 labels = NULL),
                                show_annotation_name = FALSE)

  bottom_anno <- HeatmapAnnotation(
    foo = anno_block(gp = gpar(fill = "NA",
                               col = "NA"),
                     labels = unique(geneOrder$chr),
                     labels_rot = 90,
                     labels_gp = gpar(cex = 1))
  )

  if (is.null(annotation_name)) {
    obs_anno_df <- obs_anno_df %>% rename({{obs_class}} := class)

    if (!is.null(class_cols)) {
      hist_colors <- class_cols
    } else {
      hist_colors <- setNames(hist_colors_default,
                              sort(unique(obs_anno_df[[obs_class]])))
      hist_colors <- hist_colors[!is.na((names(hist_colors)))]
    }

    anno_type <- c("Section", obs_class)
    anno_col <- list(sample_colors, hist_colors)
    names(anno_col) <- anno_type
    anno_leg <- list(
      list(
        direction = "horizontal",
        nrow = 1,
        grid_height = unit(8, "mm"),
        grid_width  = unit(8, "mm"),
        labels_gp   = gpar(fontsize = 16),
        title_gp    = gpar(fontsize = 16,
                           fontface = "bold")),
      list(
        direction = "horizontal",
        nrow = 1,
        grid_height = unit(8, "mm"),
        grid_width  = unit(8, "mm"),
        labels_gp   = gpar(fontsize = 16),
        title_gp    = gpar(fontsize = 16,
                           fontface = "bold"))
    )
    names(anno_leg) <- anno_type

    obs_left_anno <- rowAnnotation(
      df = obs_anno_df[, anno_type],
      col = anno_col,
      width = unit(19, "mm"),
      simple_anno_size_adjust = TRUE,
      gap = unit(3, "mm"),
      annotation_legend_param = anno_leg
    )
  } else {
    obs_anno_df <- obs_anno_df %>% rename({{obs_class}} := class)

    if (!is.null(class_cols)) {
      clone_colors <- class_cols
    } else {
      clone_levels <- sort(unique(obs_anno_df[[obs_class]]))
      non_x_clones_levels <- clone_levels[!str_detect(clone_levels, "^X")]
      non_x_clone_color_indices <- round(seq(1, length(clone_colors_default), length.out = length(non_x_clones_levels)))
      non_x_clone_colors <- setNames(clone_colors_default[non_x_clone_color_indices], non_x_clones_levels)
      x_clones_levels <- clone_levels[str_detect(clone_levels, "^X")]
      x_clone_color_indices <- round(seq(1, length(clone_x_colors_default), length.out = length(x_clones_levels)))
      x_clone_colors <- setNames(clone_x_colors_default[x_clone_color_indices], x_clones_levels)
      clone_colors <- c(non_x_clone_colors, x_clone_colors)[clone_levels]
    }

    if (is.character(annotation_df) && file.exists(annotation_df)) {
      hist_anno <- read.table(annotation_df, header = FALSE, row.names = 1, sep = "\t", stringsAsFactors = FALSE, col.names = c("Spot_ID_temp", annotation_name))
    } else if (is.data.frame(annotation_df)) {
      hist_anno <- annotation_df
      if (!annotation_name %in% colnames(hist_anno)) {
        stop("annotation_df data frame must contain a column named the same as 'annotation_name'")
      }
    } else {
      stop("`annotation_df` must be either a file path to a .tsv or a data frame.")
    }
    if (!all(rownames(obs_anno_df) %in% rownames(hist_anno))) {
      stop("annotation_df data frame rownames must contain all inferCNV object observation rownames")
    }

    hist_anno <- hist_anno[rownames(obs_anno_df), , drop = FALSE]

    obs_anno_df[[annotation_name]] <- hist_anno[rownames(obs_anno_df), annotation_name]

    if (!is.null(anno_cols)) {
      hist_colors <- anno_cols
    } else {
      hist_colors <- setNames(hist_colors_default,
                              sort(unique(obs_anno_df[[annotation_name]])))
      hist_colors <- hist_colors[!is.na((names(hist_colors)))]
    }

    anno_type <- c("Section", annotation_name, obs_class)
    anno_col <- list(sample_colors, hist_colors, clone_colors)
    names(anno_col) <- anno_type
    anno_leg <- list(
      list(
        direction = "horizontal",
        nrow = 1,
        grid_height = unit(8, "mm"),
        grid_width  = unit(8, "mm"),
        labels_gp   = gpar(fontsize = 16),
        title_gp    = gpar(fontsize = 16,
                           fontface = "bold")),
      list(
        direction = "horizontal",
        nrow = 1,
        grid_height = unit(8, "mm"),
        grid_width  = unit(8, "mm"),
        labels_gp   = gpar(fontsize = 16),
        title_gp    = gpar(fontsize = 16,
                           fontface = "bold")),
      list(
        direction = "horizontal",
        nrow = 1,
        grid_height = unit(8, "mm"),
        grid_width  = unit(8, "mm"),
        labels_gp   = gpar(fontsize = 16),
        title_gp    = gpar(fontsize = 16,
                           fontface = "bold"))
    )
    names(anno_leg) <- anno_type

    obs_left_anno <- rowAnnotation(
      df = obs_anno_df[, anno_type],
      col = anno_col,
      width = unit(30, "mm"),
      simple_anno_size_adjust = TRUE,
      gap = unit(3, "mm"),
      annotation_legend_param = anno_leg
    )
  }

  obs_ht <- Heatmap(
    t(obs_expr),
    col = cnv_col_fun,
    cluster_rows = dend,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_row_dend = TRUE,
    row_dend_side = "left",
    column_split = factor(sub_geneOrder$chr, levels = unique(sub_geneOrder$chr)),
    column_gap = unit(0, "mm"),
    row_split = row_split,
    cluster_row_slices = FALSE,
    row_gap = unit(0, "mm"),
    border = TRUE,
    show_heatmap_legend = FALSE,
    bottom_annotation = bottom_anno,
    top_annotation = top_anno,
    left_annotation = obs_left_anno,
    row_title = "Observations (Spots)",
    row_title_side = "right",
    row_title_rot = 90,
    row_title_gp = gpar(fontsize = 18, fontface = "bold", lineheight = 3),
    row_dend_width = unit(obs_plot_height, "cm"),
    height = unit(8, "inch"),
    use_raster = TRUE,
    raster_device = "agg_png",
    raster_quality = 6,
    raster_resize_mat = FALSE,
    raster_by_magick = FALSE
  )

  # Function for heatmap legend
  plot_infercnv_legend_static <- function(expr_data,
                                          cnv_breakpoints,
                                          cnv_colors,
                                          denscol = "blue",
                                          densadj = 0.25,
                                          key.xlab = "Modified Expression") {

    round_limits <- function(x) {
      values <- sort(unique(c(seq(0, 10, 0.10), seq(0, 10, 0.15))))
      nearest <- values[which.min(abs(values - x))]
      return(nearest)
    }

    x.range <- c(round_limits(cnv_breakpoints[1]), round_limits(cnv_breakpoints[length(cnv_breakpoints)]))

    breaks <- cnv_breakpoints
    colors <- cnv_colors

    expr_data_clamped <- pmax(expr_data, min(breaks))
    expr_data_clamped <- pmin(expr_data_clamped, max(breaks))

    h <- hist(expr_data_clamped, breaks = breaks, plot = FALSE)
    hist_data <- data.frame(
      x = rep(h$breaks, each = 2),
      y = c(0, rep(h$counts, each = 2), 0)
    )
    hist_data$y <- hist_data$y / max(hist_data$y, na.rm = TRUE) * 0.98

    color_data <- data.frame(
      xstart   = breaks[-length(breaks)],
      xend     = breaks[-1],
      midpoint = 0.5 * (breaks[-length(breaks)] + breaks[-1])
    )

    p <- ggplot() +
      geom_rect(data = color_data,
                aes(xmin = xstart, xmax = xend, ymin = 0, ymax = 1, fill = midpoint),
                color = NA) +
      scale_fill_gradientn(colors = colors,
                           limits = x.range,
                           breaks = c(x.range[1], mean(x.range), x.range[2]),
                           labels = round(c(x.range[1], mean(x.range), x.range[2]), 2),
                           guide = "none") +
      coord_cartesian(expand = FALSE) +
      xlab(key.xlab) + ylab("") +
      geom_step(data = hist_data, aes(x = x, y = y), direction = "hv", color = denscol) +
      ggtitle("Distribution of Expression") +
      theme_minimal(base_size = 12) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none")

    return(p)
  }

  # Heatmap legend
  legend_plot <- plot_infercnv_legend_static(expr, cnv_breakpoints, cnv_colors)
  legend_plot_grob <- ggplotGrob(legend_plot)

  # Annotations - Reference
  ## if there is reference
  if (is_empty(infercnv_obj@reference_grouped_cell_indices) == FALSE) {
    ref_loc <- infercnv_obj@reference_grouped_cell_indices

    names(ref_loc) <- make.names(names(ref_loc))

    suppressWarnings({
      list2env(ref_loc, envir = environment())
    })

    ref_anno_list <- list()
    for (annotation in names(ref_loc)) {
      ref_anno_list[[annotation]] <- data.frame(
        CB = colnames(expr)[get(annotation)],
        class = annotation
      )
    }
    ref_anno_df <- do.call(rbind, ref_anno_list)

    ref_expr <- expr[, colnames(expr) %in% ref_anno_df$CB]
    ref_anno_df <- ref_anno_df[match(colnames(ref_expr), ref_anno_df$CB), ]

    rownames(ref_anno_df) <- ref_anno_df$CB
    ref_anno_df$CB <- NULL

    ref_anno_df$Section <- str_remove(str_remove(rownames(ref_anno_df), "_[^_]*$"), "^[^_]*_")

    if (is.null(annotation_name)) {
      ref_anno_df[[obs_class]] <- "Benign"

      ref_left_anno <- rowAnnotation(
        df = ref_anno_df[, anno_type],
        col = anno_col,
        width = unit(19, "mm"),
        simple_anno_size_adjust = TRUE,
        gap = unit(3, "mm"),
        annotation_legend_param = anno_leg
      )
    } else {
      ref_anno_df[[obs_class]] <- "0"
      ref_anno_df[[annotation_name]] <- "Benign"

      ref_left_anno <- rowAnnotation(
        df = ref_anno_df[, anno_type],
        col = anno_col,
        width = unit(30, "mm"),
        simple_anno_size_adjust = TRUE,
        gap = unit(3, "mm"),
        annotation_legend_param = anno_leg
      )
    }

    # Ref heatmap object
    ref_ht <- Heatmap(
      t(ref_expr),
      col = cnv_col_fun,
      cluster_rows = function(mat) hclust(dist(mat), method = "ward.D2"),
      show_row_dend = FALSE,
      cluster_columns = FALSE,
      show_column_names = FALSE,
      show_row_names = FALSE,
      left_annotation = ref_left_anno,
      column_split = factor(sub_geneOrder$chr, levels = unique(sub_geneOrder$chr)),
      column_gap = unit(0, "mm"),
      column_title = NULL,
      border = TRUE,
      row_title = "References (Spots)",
      row_title_side = "right",
      row_title_rot = 90,
      row_title_gp = gpar(fontsize = 18, fontface = "bold", lineheight = 3),
      show_heatmap_legend = FALSE,
      height = unit(ref_plot_height, "inch"),
      use_raster = TRUE,
      raster_device = "agg_png",
      raster_quality = 6,
      raster_resize_mat = FALSE,
      raster_by_magick = FALSE
    )

    # Ref pie plot
    ref_section_counts <- ref_anno_df %>%
      count(Section) %>%
      rename(group = Section, value = n)

    ref_section_counts <- ref_section_counts %>%
      arrange(desc(group)) %>%
      mutate(
        prop = value / sum(value) * 100,
        ypos_center = cumsum(prop) - 0.5 * prop
      )

    ref_pie_plot <- ggplot(ref_section_counts, aes(x = "", y = prop, fill = group)) +
      geom_bar(stat = "identity", width = 1, color = "white") +
      coord_polar("y", start = 0) +
      theme_void() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14.4),
        plot.margin = unit(c(0.5,0,0,0), "cm")
      ) +
      geom_label_repel(aes(y = ypos_center, label = paste0(round(prop, 1), "%")),
                       nudge_x = 1,
                       segment.colour = "black",
                       fill = "white",
                       color = "black",
                       size = 4,
                       box.padding = unit(0.5, "lines"),
                       point.padding = unit(0.5, "lines"),
                       show.legend = FALSE
      ) +
      scale_fill_manual(values = sample_colors, name = "Sections") +
      labs(title = "Distribution of reference")

    ref_pie_grob <- ggplotGrob(ref_pie_plot)

  } else {
    ref_ht <- Heatmap(
      matrix(NA, nrow = 1, ncol = nrow(obs_expr)),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_heatmap_legend = FALSE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      column_split = factor(sub_geneOrder$chr, levels = unique(sub_geneOrder$chr)),
      column_gap = unit(0, "mm"),
      column_title = NULL,
      border = FALSE,
      rect_gp = gpar(fill = NA, col = NA),
      height = unit(ref_plot_height, "inch")
    )
  }

  ht_list <- ref_ht %v% obs_ht

  if (output_format == "png") {
    png(file.path(output_dir, output_name), width = 20, height = 15, units = "in", res = 300)
    draw(ht_list,
         padding = unit(c(1, 1, 2, 1), "cm"),
         column_title = "Genomic Region",
         column_title_side = "bottom",
         column_title_gp = gpar(fontsize = 18, fontface = "bold", lineheight = 3),
         ht_gap = unit(2, "cm"),
         annotation_legend_side = "bottom"
    )
    if (isTRUE(ref_pie_chart)) {
      pie_vp <- viewport(x = 0.03, y = 0.87, width = 0.15, height = 0.15, just = c("left", "top"))
      pushViewport(pie_vp)
      grid.draw(ref_pie_grob)
      popViewport()
    }
    legend_vp <- viewport(x = 0.05, y = 0.97, width = 0.10, height = 0.10, just = c("left", "top"))
    pushViewport(legend_vp)
    grid.draw(legend_plot_grob)
    popViewport()
    dev.off()
  } else if (output_format == "pdf") {
    pdf(file.path(output_dir, output_name), width = 20, height = 15)
    draw(ht_list,
         padding = unit(c(1, 1, 2, 1), "cm"),
         column_title = "Genomic Region",
         column_title_side = "bottom",
         column_title_gp = gpar(fontsize = 18, fontface = "bold", lineheight = 3),
         ht_gap = unit(2, "cm"),
         annotation_legend_side = "bottom"
    )
    if (isTRUE(ref_pie_chart)) {
      pie_vp <- viewport(x = 0.03, y = 0.87, width = 0.15, height = 0.15, just = c("left", "top"))
      pushViewport(pie_vp)
      grid.draw(ref_pie_grob)
      popViewport()
    }
    legend_vp <- viewport(x = 0.05, y = 0.97, width = 0.10, height = 0.10, just = c("left", "top"))
    pushViewport(legend_vp)
    grid.draw(legend_plot_grob)
    popViewport()
    dev.off()
  }
}
