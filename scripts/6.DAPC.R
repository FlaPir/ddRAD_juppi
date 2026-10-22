########################################
# 1.DAPC                               #
# Author: Flavio Piras                 #
# Reference: ....                      #
#                                      #
########################################


#RStudio 2025.05.1+513 "Mariposa Orchid" Release (ab7c1bc795c7dcff8f26215b832a3649a19fc16c, 2025-06-01) for windows

library(vcfR)
library(adegenet)


# 1. SELECT K MODE
# -----------------------------------------------------
# "metadata" -> Supplied a priori from sample_locality.txt
# "manual_2" -> Inferred forcing K = 2 via find.clusters()
# "auto_bic" -> Inferred choosing optimal K via BIC
MODO_K <- "auto_bic"

# Working and output directories
dir_input  <- "path/to/your/input_folder"
dir_output <- "path/to/your/output_folder"

if(!dir.exists(dir_output)) dir.create(dir_output, recursive=TRUE)
setwd(dir_input)

# 2. IMPORT AND DATA PREPARATION
# -----------------------------------------------------
vcf <- read.vcfR("Dataset_B.vcf")
gl_data <- vcfR2genlight(vcf)

info_pop <- read.csv("sample_locality.txt", sep="\t", stringsAsFactors=FALSE) #two column file tab separated; header Sample\tLocality; e.g. BO11_23\tBologna
colnames(info_pop) <- c("Sample","Locality")

# Trim whitespace
info_pop[["Sample"]]   <- trimws(info_pop[["Sample"]])
info_pop[["Locality"]] <- trimws(info_pop[["Locality"]])

# Sample name matching (Extracts "BO11_23" from "Bologna_BO11_23")
vcf_raw_names <- trimws(indNames(gl_data))
vcf_clean_names <- sapply(strsplit(vcf_raw_names, "_"), function(x) {
  paste(tail(x, 2), collapse = "_")
})

# Match clean IDs
matched_idx <- match(vcf_clean_names, info_pop[["Sample"]])
info_pop_matched <- info_pop[matched_idx, ]

# Max PCs allowed for clustering based on sample count
max_pca_find <- min(200, nInd(gl_data) - 1)

# Population / Group Assignment
if (MODO_K == "manual_2") {
  cat("\n[INFO] Mode: K = 2 (Inferred via find.clusters)\n")
  set.seed(123)
  grp <- find.clusters(gl_data, n.clust = 2, n.pca = max_pca_find, choose.n.clust = FALSE)
  pop(gl_data) <- grp[["grp"]]
  
} else if (MODO_K == "auto_bic") {
  cat("\n[INFO] Mode: Optimal K (Inferred via BIC with find.clusters)\n")
  set.seed(123)
  grp <- find.clusters(gl_data, max.n.clust = 10, n.pca = max_pca_find, choose.n.clust = FALSE, criterion = "min")
  pop(gl_data) <- grp[["grp"]]
  
} else {
  cat("\n[INFO] Mode: Metadata Populations (Supplied a priori)\n")
  pop_levels <- unique(info_pop_matched[["Locality"]])
  pop_levels <- na.omit(pop_levels)
  if ("Lido_delle_Nazioni_J" %in% pop_levels) {
    pop_levels <- c(setdiff(pop_levels, "Lido_delle_Nazioni_J"), "Lido_delle_Nazioni_J")
  }
  pop(gl_data) <- factor(info_pop_matched[["Locality"]], levels=pop_levels)
}

# Safety check for population assignment
cat("\n--- SAMPLE COUNT PER POPULATION ---\n")
print(table(pop(gl_data), useNA = "always"))

if (any(is.na(pop(gl_data)))) {
  cat("Warning: Removing", sum(is.na(pop(gl_data))), "samples with NA population...\n")
  gl_data <- gl_data[!is.na(pop(gl_data))]
  pop(gl_data) <- droplevels(pop(gl_data))
}

pop_levels <- levels(pop(gl_data))
K <- length(pop_levels)

# 3. COLOR PALETTE
# -----------------------------------------------------
if (MODO_K == "metadata") {
  reds <- colorRampPalette(c("#F6B7B7","#C44536","#7B0000"))(K)
  group_cols <- setNames(reds, pop_levels)
  if ("Lido_delle_Nazioni_J" %in% pop_levels) {
    group_cols["Lido_delle_Nazioni_J"] <- "#1F4E79"
  }
} else {
  cols <- colorRampPalette(c("#1F4E79", "#C44536", "#2E7D32", "#E65100", "#6A1B9A"))(K)
  group_cols <- setNames(cols, pop_levels)
}

# 4. CROSS-VALIDATION & FIT DAPC
# -----------------------------------------------------
tiff(file.path(dir_output, "DAPC_xval.tif"), width=7, height=6, units="in", res=600, compression="lzw")
set.seed(123)

max_pca_allowed <- min(300, nInd(gl_data) - 1)

xval <- xvalDapc(tab(gl_data, NA.method="mean"), pop(gl_data),
                 n.pca.max=max_pca_allowed, training.set=0.9, result="groupMean",
                 center=TRUE, scale=FALSE, n.rep=50)
dev.off()

n_pc_opt <- as.numeric(xval[["Number of PCs Achieving Highest Mean Success"]])
if(is.na(n_pc_opt) || n_pc_opt < (K+2)) {
  n_pc <- K + 2
} else {
  n_pc <- n_pc_opt
}

n_da <- max(1, K - 1)
dapc_final <- dapc(gl_data, n.pca=n_pc, n.da=n_da)

# PRINT REVISION REPORT
# -----------------------------------------------------
cat("\n=====================================================\n")
cat("          REVISION PARAMETER REPORT                  \n")
cat("=====================================================\n")
cat(" -> Group Assignment Type               :", ifelse(MODO_K=="metadata", "Supplied (a priori)", "Inferred (a posteriori)"), "\n")
cat(" -> Number of Populations/Clusters (K)  :", K, "\n")
cat(" -> PCs Suggested by Cross-Validation   :", ifelse(is.na(n_pc_opt), "N/A", n_pc_opt), "\n")
cat(" -> PCs ACTUALLY RETAINED IN DAPC       :", dapc_final[["n.pca"]], "\n")
cat(" -> Discriminant Functions (DA) Used    :", dapc_final[["n.da"]], "\n")
cat("=====================================================\n\n")

# =====================================================
# FILE 1: DAPC SCATTER 2D (Only if K >= 3)
# =====================================================
if (dapc_final[["n.da"]] >= 2) {
  tiff(file.path(dir_output, "DAPC_scat.tif"), width=10, height=7, units="in", res=600, compression="lzw")
  par(mar=c(5, 5, 4, 13))
  
  coords_2d <- dapc_final[["ind.coord"]][, 1:2]
  centroids_2d <- dapc_final[["grp.coord"]][, 1:2]
  pop_factor <- pop(gl_data)
  pop_levs <- levels(pop_factor)
  
  all_x <- c(coords_2d[,1], centroids_2d[,1])
  all_y <- c(coords_2d[,2], centroids_2d[,2])
  
  extend_lims <- function(vec, pad=0.2) {
    r <- range(vec)
    diff_r <- diff(r)
    return(c(r[1] - diff_r * pad, r[2] + diff_r * pad))
  }
  
  prop_eig <- dapc_final[["eig"]] / sum(dapc_final[["eig"]]) * 100
  
  plot(NA, xlim=extend_lims(all_x), ylim=extend_lims(all_y),
       xlab=paste0("LD1 (", round(prop_eig[1], 2), "% discriminant variation)"),
       ylab=paste0("LD2 (", round(prop_eig[2], 2), "% discriminant variation)"),
       main=paste0("DAPC Analysis (K = ", K, ")"), cex.lab=1.1, axes=FALSE)
  axis(1); axis(2); box(lwd=1.2)
  
  # Spider legs
  for(i in pop_levs){
    ind_idx <- which(pop_factor == i)
    segments(coords_2d[ind_idx, 1], coords_2d[ind_idx, 2],
             centroids_2d[i, 1], centroids_2d[i, 2],
             col = adjustcolor(group_cols[i], 0.3), lwd = 1.1)
  }
  
  # Ellipses
  for(i in pop_levs){
    sub_coords <- coords_2d[pop_factor == i, , drop=FALSE]
    if(nrow(sub_coords) > 3){
      mu <- colMeans(sub_coords); sigma <- cov(sub_coords)
      if(det(sigma) > 1e-10){
        Q <- try(chol(sigma), silent=TRUE)
        if(!inherits(Q, "try-error")){
          angles <- seq(0, 2*pi, length.out=100)
          unit.circle <- cbind(cos(angles), sin(angles))
          ell_points <- t(t(unit.circle %*% Q) + mu)
          lines(ell_points, col=adjustcolor(group_cols[i], 0.6), lwd=1.5)
        }
      }
    }
  }
  
  points(coords_2d[,1], coords_2d[,2], col=adjustcolor(group_cols[pop_factor], 0.8), pch=19, cex=1.3)
  points(centroids_2d[,1], centroids_2d[,2], pch=23, bg=group_cols[rownames(centroids_2d)], col="black", cex=2.5, lwd=1.5)
  
  # External legend
  legend("topleft", inset=c(1.02, 0), xpd=TRUE, legend=pop_levs, col=group_cols,
         pch=19, pt.cex=1.4, bty="n", cex=0.9, title="Populations / Clusters")
  
  dev.off()
} else {
  cat("ℹ️ K = 2: 2D Scatterplot omitted (only LD1 axis exists).\n")
}

# =====================================================
# FILE 2: DENSITY LD1
# =====================================================
tiff(file.path(dir_output, "DAPC_dens.tif"), width=8, height=5, units="in", res=600, compression="lzw")
par(mar=c(5,5,3,2))

pop_factor <- pop(gl_data)
ld1 <- dapc_final[["ind.coord"]][,1]
dens_list <- lapply(levels(pop_factor), function(p) density(ld1[pop_factor==p], bw="nrd0", adjust=1.2))

ymax <- max(sapply(dens_list, function(d) max(d[["y"]])))
xmin <- min(ld1) - (diff(range(ld1))*0.25)
xmax <- max(ld1) + (diff(range(ld1))*0.25)

plot(NA, xlim=c(xmin, xmax), ylim=c(0, ymax*1.1),
     xlab="Discriminant Function 1 (LD1)", ylab="Density",
     main=paste0("Population Separation along LD1 (K = ", K, ")"), frame=FALSE)
abline(h=0, col="grey80")

for(i in seq_along(dens_list)){
  polygon(dens_list[[i]][["x"]], dens_list[[i]][["y"]],
          col=adjustcolor(group_cols[levels(pop_factor)[i]], 0.4),
          border=group_cols[levels(pop_factor)[i]], lwd=2)
}
dev.off()

# =====================================================
# FILE 3: STANDALONE LEGEND
# =====================================================
tiff(file.path(dir_output, "DAPC_leg.tif"), width=4, height=6, units="in", res=600, compression="lzw")
plot.new()
legend("center", legend=levels(pop(gl_data)), col=group_cols,
       pch=19, pt.cex=1.5, bty="n", cex=1.1, title="Populations / Clusters")
dev.off()

# =====================================================
# FILE 4: SUPPLEMENTARY - UNSUPERVISED PCA
# =====================================================
cat("\nCalculating unsupervised PCA for Supplementary Material...")
pca_unsup <- glPca(gl_data, nf = 3)
pca_var <- (pca_unsup[["eig"]] / sum(pca_unsup[["eig"]])) * 100

tiff(file.path(dir_output, "PCA_unsupervised_Supp.tif"), width=10, height=7, units="in", res=600, compression="lzw")
par(mar=c(5, 5, 4, 13))

plot(pca_unsup[["scores"]][, 1], pca_unsup[["scores"]][, 2],
     col = adjustcolor(group_cols[pop(gl_data)], 0.8),
     pch = 19, cex = 1.4,
     xlab = paste0("PC1 (", round(pca_var[1], 2), "% total genetic variance)"),
     ylab = paste0("PC2 (", round(pca_var[2], 2), "% total genetic variance)"),
     main = "Unsupervised PCA (No A Priori Group Maximization)",
     axes = FALSE)
axis(1); axis(2); box(lwd=1.2)

# External legend
legend("topleft", inset=c(1.02, 0), xpd=TRUE, legend=levels(pop(gl_data)),
       col=group_cols, pch=19, pt.cex=1.4, bty="n", cex=0.9, title="Populations / Clusters")

dev.off()

# =====================================================
# FILE 5: BIC EVALUATION (When MODO_K is "auto_bic" or "manual_2")
# =====================================================
if (MODO_K %in% c("auto_bic", "manual_2")) {
  
  cat("\n=====================================================\n")
  cat("           BIC STATISTICS FOR K SELECTION            \n")
  cat("=====================================================\n")
  print(grp[["Kstat"]])
  cat("=====================================================\n")
  
  tiff(file.path(dir_output, "BIC_curve_Supp.tif"), width=6, height=5, units="in", res=600, compression="lzw")
  plot(grp[["Kstat"]], type="b", col="#1F4E79", pch=19, lwd=2,
       xlab="Number of Clusters (K)", ylab="BIC",
       main="BIC Value by Number of Clusters (K)")
  dev.off()
}



dapc_final[["n.pca"]]

cat("\n\n✅ EXECUTION COMPLETED SUCCESSFULLY!")
cat("\nOutput files saved in folder:", dir_output, "\n")
