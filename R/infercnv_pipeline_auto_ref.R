#!/usr/bin/env Rscript
# infercnv_bestk_auto.R
#
# * Per-patient, per-mode (filtered / nofilter)
# * First try to read existing k_by_cnvscore_k_results.rds
# * Otherwise continue k-means from max(k)+1 until a valid k is found
# * Plot normalized size vs. CNV variance line+points, save PNG + CSV
# * Use the cleanest cluster at best-k as reference barcodes
# * Run inferCNV on all epithelial spots with that reference
# --------------------------------------------------------------------

suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(infercnv)
  library(parallel)
  library(ggplot2)
  library(viridis)
  library(readr)
})

# ───── Argument parsing ─────────────────────────────────────────────
parser <- ArgumentParser()
parser$add_argument("--patients",  required=TRUE,
                    help="Comma-separated list of patient IDs, e.g., pt03,pt13")
parser$add_argument("--section_outdir",
                    default="/srv/home/jintong.shi/Prostate_Oxford/data/visium_seurat_object/250708",
                    help="Directory containing per-section Seurat objects")
parser$add_argument("--output_base",
                    default="/srv/home/jintong.shi/Prostate_Oxford/processed_data/visium/infercnv",
                    help="Base directory for inferCNV outputs")
parser$add_argument("--gene_pos_file",
                    default="/srv/home/jintong.shi/RefData_public/gencode_v38_gene_pos.txt",
                    help="Gene position file for inferCNV")
parser$add_argument("--min_feats",   type="integer", default=500,
                    help="Minimum features per spot for filtering")
parser$add_argument("--min_counts",  type="integer", default=1000,
                    help="Minimum counts per spot for filtering")
parser$add_argument("--cutoff",      type="double",  default=0.1,
                    help="inferCNV cutoff")
parser$add_argument("--ncores",      type="integer", default=8,
                    help="Number of threads for inferCNV")
parser$add_argument("--prop_min",    type="double",  default=0.10,
                    help="Minimum cluster size as fraction of benign total")
parser$add_argument("--k_window",    type="integer", default=2,
                    help="Window size for valley detection (k_window=2 means 2 before and after)")
args <- parser$parse_args()

patients       <- unlist(strsplit(args$patients,",")) |> trimws()
section_outdir <- args$section_outdir
output_base    <- args$output_base
gene_pos_file  <- args$gene_pos_file
min_feats      <- args$min_feats
min_counts     <- args$min_counts
cutoff         <- args$cutoff
ncores         <- args$ncores
prop_min       <- args$prop_min
k_window       <- args$k_window     # Number of ks before/after for valley test

# ───── Read gene position file once ─────────────────────────────────
geneFile <- read.table(gene_pos_file, header=FALSE, sep="\t", stringsAsFactors = FALSE) 

# ───── Helper: run inferCNV ──────────────────────────────────────────
run_infercnv <- function(expr_mat, annot_df, ref_names, out_dir){
  obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = expr_mat,
    gene_order_file   = gene_pos_file,
    annotations_file  = annot_df,
    delim             = "\t",
    ref_group_names   = ref_names,
    chr_exclude       = "chrM"
  )
  infercnv::run(
    obj, cutoff=cutoff, out_dir=out_dir,
    denoise=TRUE, HMM=FALSE,
    cluster_by_groups = FALSE,
    write_phylo = TRUE,
    analysis_mode = "cells",
    num_threads = ncores
  )
}

run_infercnv_noplot <- function(expr_mat, annot_df, ref_names, out_dir){
  # Same as above but skip the default plotting step
  obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = expr_mat,
    gene_order_file   = gene_pos_file,
    annotations_file  = annot_df,
    delim             = "\t",
    ref_group_names   = ref_names,
    chr_exclude       = "chrM"
  )
  infercnv::run(
    obj, cutoff=cutoff, out_dir=out_dir,
    denoise=TRUE, HMM=FALSE,
    cluster_by_groups = FALSE,
    write_phylo = TRUE,
    analysis_mode = "cells",
    num_threads = ncores,
    no_plot = TRUE
  )
  # After running, attempt custom plotting if threshold/dendrogram files are missing
  inferCNV_object <- readRDS(file.path(out_dir,'run.final.infercnv_obj'))
  thresh_file <- file.path(out_dir,"infercnv.heatmap_thresholds.txt")
  den_file    <- file.path(out_dir,"infercnv.observations_dendrogram.txt")
  if (file.exists(thresh_file) && file.exists(den_file)) {
    message("→ thresholds & dendrogram exist, skipping plot_cnv()")
  } else {
    tryCatch({
      plot_cnv(
        inferCNV_object,
        out_dir=out_dir,
        plot_chr_scale=FALSE,
        contig_cex = 2,
        title = 'infercnv',
        obs_title="Observations (Cells)",
        ref_title="References (Cells)",
        cluster_by_groups=FALSE,
        x.center=1,
        x.range="auto",
        hclust_method='ward.D',
        color_safe_pal=FALSE,
        output_filename='infercnv',
        output_format="png",
        dynamic_resize=1,
        write_phylo = TRUE
      )
    }, error = function(e) {
      message("⚠ plot_cnv() error, proceeding: ", e$message)
    })
  }
  # Then plot the ComplexHeatmap version
  start_time <- Sys.time()
  plot_complex_heatmap(
    infercnv_obj = inferCNV_object,
    output_dir = out_dir,
    heatmap_thresholds_file = "infercnv.heatmap_thresholds.txt",
    dendrogram_file = "infercnv.observations_dendrogram.txt",
    obs_class = "histology",
    ref_pie_chart = TRUE
  )
  end_time <- Sys.time()
  message("Complex heatmap took ", round(difftime(end_time, start_time, units="mins"), 2), " minutes")
}

# ───── Main loop ────────────────────────────────────────────────────
for (pt in patients) {
  for (mode in c("filtered")) {# "nofilter" also possible
    message("\n========== Processing ", pt, " / ", mode, " ==========")
    
    # 0) set up paths and read Seurat object
    sel_dir   <- file.path(output_base, pt, mode, "benign_no-ref", "k_selection")
    dir.create(sel_dir, recursive=TRUE, showWarnings=FALSE)
    kres_file <- file.path(sel_dir, "k_by_cnvscore_k_results.rds")
    
    seu <- readRDS(file.path(section_outdir, paste0(pt,"_allseu.rds.gz")))
    epi_keep <- subset(seu, subset = histology_V2 %in%
                         c("Benign","Cancer","GG2","GG3","GG4","GG5",
                           "Metastatic","IDC","HGPIN","PIN","GG1"))
    if (mode=="filtered") {
      epi_keep <- subset(epi_keep,
                         subset = nFeature_Spatial>=min_feats &
                                  nCount_Spatial>=min_counts)
    }
    
    # 1) extract benign count matrix
    benign <- subset(epi_keep, subset = histology_V2 == "Benign")
    mat_b  <- GetAssayData(benign, assay="Spatial", slot="counts")
    n_b    <- ncol(mat_b)
    if (n_b < 100) {
      message("   * Fewer than 100 benign spots, skipping")
      next
    }
    
    # 2) run benign-no-ref inferCNV once if needed
    dir_bnr <- file.path(output_base, pt, mode, "benign_no-ref")
    if (!file.exists(file.path(dir_bnr,"run.final.infercnv_obj"))) {
      message("   * Running benign-no-ref inferCNV...")
      ann_b <- data.frame(sample=benign$section_name, row.names=colnames(mat_b))
      run_infercnv(mat_b, ann_b, NULL, dir_bnr)
    }
    inf_b    <- readRDS(file.path(dir_bnr,"run.final.infercnv_obj"))
    expr.data <- inf_b@expr.data
    cnv_mat   <- t(expr.data)                 # cells × genes
    cnv_var   <- colMeans((expr.data - 1)^2)  # per-cell CNV score
    
    # 3) k-means loop (reuse old results if available)
    kres   <- list()
    stat   <- tibble()  # will store k / size / cnv
    found  <- FALSE
    k_start <- 2
    
    # 3.1 load existing kres if present
    if (file.exists(kres_file)) {
      message("   * Loading existing k-results...")
      kres_loaded <- readRDS(kres_file)
      for (kn in names(kres_loaded)) {
        kres[[kn]] <- kres_loaded[[kn]]
        best_cl <- kres_loaded[[kn]]$best_cl
        sums    <- kres_loaded[[kn]]$sums
        stat <- bind_rows(stat,
                          tibble(k   = kres_loaded[[kn]]$k,
                                 size= sums$size[sums$cluster==best_cl],
                                 cnv = sums$avg [sums$cluster==best_cl]))
      }
      # check for valley in existing data
      if (nrow(stat) >= (1 + 2 * k_window)) {
        for (idx in (k_window+1):(nrow(stat)-k_window)) {
          cur   <- stat$cnv[idx]
          prevs <- stat$cnv[(idx-k_window):(idx-1)]
          nexts <- stat$cnv[(idx+1):(idx+k_window)]
          if (all(cur < prevs) && all(cur < nexts) &&
              stat$size[idx] >= prop_min * n_b) {
            k_best <- stat$k[idx]
            found  <- TRUE
            message("   * Found best-k in existing data: k", k_best)
            break
          }
        }
      }
      k_start <- ifelse(found, NA, max(stat$k) + 1)
    }
    
    # 3.2 if not found, extend k-means until valley found
    k <- k_start
    while (!found) {
      message("      Computing k = ", k)
      km <- kmeans(cnv_mat, centers=k)
      df <- tibble(
        barcode = rownames(cnv_mat),
        cluster = factor(km$cluster, levels=1:k),
        cnv     = cnv_var
      )
      sums    <- df %>% group_by(cluster) %>% summarize(
        avg  = mean(cnv),
        size = n(),
        .groups = "drop"
      )
      best_cl <- sums %>% slice_min(avg) %>% pull(cluster) %>% as.character()
      kres[[paste0("k",k)]] <- list(k=k, df=df, sums=sums, best_cl=best_cl)
      stat <- bind_rows(stat,
                        tibble(k=k,
                               size = sums$size[sums$cluster==best_cl],
                               cnv  = sums$avg [sums$cluster==best_cl]))
      # valley detection
      if (nrow(stat) >= (1+2*k_window)) {
        idx <- nrow(stat) - k_window
        if (idx > k_window) {
          cur   <- stat$cnv[idx]
          prevs <- stat$cnv[(idx-k_window):(idx-1)]
          nexts <- stat$cnv[(idx+1):(idx+k_window)]
          if (all(cur < prevs) && all(cur < nexts) &&
              stat$size[idx] >= prop_min * n_b) {
            k_best <- stat$k[idx]
            found  <- TRUE
            message("   * Found best-k by extension: k", k_best)
          }
        }
      }
      k <- k + 1
      if (k > 40) {
        message("   ** Reached k=40 without finding valley, skipping")
        break
      }
    }
    
    # save updated kres
    saveRDS(kres, kres_file)
    if (!found) next
    
    # 4) plot normalized size vs CNV variance up to best_k + window
    stat_norm <- stat %>%
      mutate(
        cnv_norm  = cnv   / max(cnv),
        size_prop = size  / n_b
      )
    stat_plot <- stat_norm %>% filter(k <= k_best + k_window)
    write_csv(stat_plot,
              file.path(sel_dir, paste0("cnvVar_up_to_bestKplus",k_window,".csv")))
    p <- ggplot(stat_plot, aes(x=k, y=cnv_norm, group=1)) +
      geom_line(color="royalblue") + geom_point(size=3) +
      geom_text(aes(label=k), vjust=-0.5) +
      xlab("k") + ylab("Normalized CNV variance") +
      ggtitle(sprintf("%s / %s", pt, mode)) +
      theme_minimal()
    ggsave(file.path(sel_dir, paste0("cnvVar_by_k_window",k_window,".png")),
           p, width=6, height=4)
    message("   * Size vs CNV-var plot saved")
    
    # 5) extract best-cluster barcodes as reference
    best_r     <- kres[[paste0("k",k_best)]]
    bcand      <- best_r$df$barcode[ best_r$df$cluster==best_r$best_cl ]
    ref_out_dir <- file.path(output_base, pt, mode,
                             paste0("all_epi_with_ref_k",k_best))
    dir.create(ref_out_dir, recursive=TRUE, showWarnings=FALSE)
    write_lines(bcand, file.path(ref_out_dir,"ref_cells.txt"))
    
    # 6) run all-epi inferCNV with that reference
    if (file.exists(file.path(ref_out_dir,'run.final.infercnv_obj'))) {
      message("   * inferCNV-with-ref already exists (k=",k_best,")")
    } else {
      expr_all <- GetAssayData(epi_keep, assay="Spatial", slot="counts")
      meta_all <- epi_keep@meta.data %>%
        mutate(group = ifelse(rownames(.) %in% bcand,
                              "pure_benign", histology_V2))
      ann_all  <- data.frame(group=meta_all$group,
                             row.names=rownames(meta_all))
      message("   * Running all-epi inferCNV-with-ref (k=",k_best,")")
      if (pt %in% c('pt02','pt10','pt01')) {
        run_infercnv_noplot(expr_all, ann_all, "pure_benign", ref_out_dir)
      } else {
        run_infercnv(expr_all, ann_all, "pure_benign", ref_out_dir)
      }
    }
    message("   * Finished: ", ref_out_dir)
  } # end mode loop
}   # end patient loop

message("\n===== All patients/modes processed =====")
