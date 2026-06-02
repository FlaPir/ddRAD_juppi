#RStudio 2025.05.1+513 "Mariposa Orchid" Release (ab7c1bc795c7dcff8f26215b832a3649a19fc16c, 2025-06-01) for windows

library(vcfR)
library(adegenet)
library(RColorBrewer)

setwd("path/to/VCF")
vcf <- read.vcfR("dataset.vcf")
gl_data <- vcfR2genlight(vcf)

info_pop <- read.csv("path/to/popmap",  #sample - population file
                    sep="\t", stringsAsFactors=FALSE)
colnames(info_pop) <- c("Sample","Locality")
info_pop$Sample <- trimws(info_pop$Sample)
info_pop$Population <- trimws(info_pop$Locality)
info_pop <- info_pop[match(indNames(gl_data), info_pop$Sample), ]

pop_levels <- unique(info_pop$Locality)

pop(gl_data) <- factor(info_pop$Population, levels=pop_levels)
K <- length(pop_levels)

if (K <= 9) {
  my_cols <- brewer.pal(K, "Set1") 
} else {
  my_cols <- colorRampPalette(brewer.pal(9, "Set1"))(K)
}
group_cols <- setNames(my_cols, pop_levels)
if("Lido_delle_Nazioni_J" %in% names(group_cols)) {
  group_cols["Lido_delle_Nazioni_J"] <- "#1F4E79" 
}

tiff("DAPC_CrossValidation.tif", 
      width=7, 
      height=6, 
      units="in", 
      res=600, 
      compression="lzw")
      
set.seed(123)
dev.off()

n_pc <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`)
if(is.na(n_pc) || n_pc < (K+2)) { n_pc <- K+2 }
dapc_final <- dapc(gl_data, n.pca=n_pc, n.da=K-1)

#scatterplot

tiff("DAPC_scat.tif", width=8, height=8, units="in", res=600, compression="lzw")
par(mar=c(5, 5, 4, 2)) 

coords_2d <- dapc_final$ind.coord[, 1:2]
centroids_2d <- dapc_final$grp.coord[, 1:2]
pop_factor <- pop(gl_data)
pop_levs <- levels(pop_factor)

all_x <- c(coords_2d[,1], centroids_2d[,1])
all_y <- c(coords_2d[,2], centroids_2d[,2])

extend_lims <- function(vec, pad=0.3) {
  r <- range(vec)
  diff_r <- diff(r)
  return(c(r[1] - diff_r * pad, r[2] + diff_r * pad))
}

plot(NA, xlim=extend_lims(all_x), ylim=extend_lims(all_y), 
     xlab=paste0("LD1 (", round(dapc_final$eig[1]/sum(dapc_final$eig)*100, 2), "%)"),
     ylab=paste0("LD2 (", round(dapc_final$eig[2]/sum(dapc_final$eig)*100, 2), "%)"),
     main="DAPC Analysis", cex.lab=1.2, axes=FALSE)
axis(1); axis(2); box(lwd=1.2)

for(i in pop_levs){
  ind_idx <- which(pop_factor == i)
  segments(coords_2d[ind_idx, 1], coords_2d[ind_idx, 2],
           centroids_2d[i, 1], centroids_2d[i, 2],
           col = adjustcolor(group_cols[i], 0.3), lwd = 1.1)
}

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

points(coords_2d[,1], 
        coords_2d[,2], 
        col=adjustcolor(group_cols[pop_factor], 
        0.8), 
        pch=19, 
        cex=1.3)
points(centroids_2d[,1], 
        centroids_2d[,2], 
        pch=23, 
        bg=group_cols[rownames(centroids_2d)], 
        col="black", 
        cex=2.5, 
        lwd=1.5)
dev.off()

#density plot

tiff("DAPC_dens.tif", width=8, height=5, units="in", res=600, compression="lzw")
par(mar=c(5,5,3,2))

ld1 <- dapc_final$ind.coord[,1]
dens_list <- lapply(levels(pop_factor), function(p) density(ld1[pop_factor==p], bw="nrd0", adjust=1.2))

ymax <- max(sapply(dens_list, function(d) max(d$y)))
xmin <- min(ld1) - (diff(range(ld1))*0.25) # Padding orizzontale
xmax <- max(ld1) + (diff(range(ld1))*0.25)

plot(NA, xlim=c(xmin, xmax), ylim=c(0, ymax*1.1),
     xlab="Discriminant Function 1 (LD1)", ylab="Density",
     main="Population Separation along LD1", frame=FALSE)
abline(h=0, col="grey80")

for(i in seq_along(dens_list)){
  polygon(dens_list[[i]]$x, dens_list[[i]]$y,
          col=adjustcolor(group_cols[levels(pop_factor)[i]], 0.4),
          border=group_cols[levels(pop_factor)[i]], lwd=2)
}
dev.off()

#legend

tiff("DAPC_PROVA_leg.tif", 
      width=4, 
      height=6, 
      units="in", 
      res=600, 
      compression="lzw")
plot.new()
legend("center", 
      legend=levels(pop(gl_data)), 
      col=group_cols, 
      pch=19, 
      pt.cex=1.5, 
      bty="n", 
      cex=1.1, 
      title="Populations")
dev.off()

