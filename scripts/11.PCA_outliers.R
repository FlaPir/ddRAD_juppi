#RStudio 2025.05.1+513 "Mariposa Orchid" Release (ab7c1bc795c7dcff8f26215b832a3649a19fc16c, 2025-06-01) for windows

library(adegenet)
library(vcfR)
library(ggplot2)
library(dplyr)

setwd("path/to/datasetA.vcf")

vcf <- read.vcfR("datasetA.vcf")
obj <- vcfR2genlight(vcf)

vcf_fix <- as.data.frame(getFIX(vcf), stringsAsFactors = FALSE)
snp_map <- vcf_fix %>%
  mutate(
    SNP_ID = if_else(is.na(ID) | ID == ".", paste(CHROM, POS, sep="_"), ID),
    Position = as.numeric(POS)
  ) %>%
  dplyr::select(SNP_ID, Chromosome = CHROM, Position)

pca_results <- glPca(obj, nf = 2)
pc1_loadings <- pca_results$loadings[,1]
z_scores <- as.numeric(scale(pc1_loadings))
pca_all_snps <- data.frame(
  Index = 1:length(pc1_loadings),
  SNP_ID = snp_map$SNP_ID,
  Chromosome = snp_map$Chromosome,
  Position = snp_map$Position,
  Loading = pc1_loadings,
  Z_score = z_scores
)
threshold <- 3
pca_outliers_only <- pca_all_snps %>% 
  filter(abs(Z_score) > threshold) %>% 
  arrange(desc(abs(Z_score)))

#visualisation

p_pca <- ggplot(pca_all_snps, aes(x = Index, y = Z_score)) +
  geom_point(alpha = 0.4, color = "grey") +
  geom_point(data = pca_outliers_only, color = "orange", size = 1.5) +
  geom_hline(yintercept = c(threshold), linetype = "dashed", color = "red") +
  labs(
    title = "PCA Loading Plot - Outliers",
    subtitle = paste("Z-score >", threshold),
    x = "SNP",
    y = "Loadings Z-score"
  ) +
  theme_minimal()

ggsave("PCA_Top_Outliers.tif", 
        plot = p_pca, 
        width = xx, 
        height = xx, 
        units = "in", 
        dpi = 600, 
        compression = "lzw")

#csv and bed file creation

bed_outliers <- pca_outliers_only %>%
  mutate(
    Start = Position - 1,
    End = Position
  ) %>%
  dplyr::select(Chromosome, Start, End, SNP_ID, Z_score)

write.csv(pca_outliers_only, "TOP_pca_outliers.csv", row.names = FALSE)
write.table(bed_outliers, "TOP_pca_outliers.bed", 
            sep="\t", 
            quote=FALSE, 
            row.names=FALSE, 
            col.names=FALSE)
