library(readxl)
library(writexl)
library(reshape)
library(ggplot2)

# This part of the code (from ### to ###) processes the excel table from Goncalves et al 2022 and calculates correlations.
# It was copied from  https://github.com/NatashaKochanova/Chromosome-instability-correlation-analysis/
###
setwd("/Users/natalia.kochanova/Desktop/Chromosome structure lab/Bioinformatics/R/202207_Correlation analysis/Supplementary/")
matrix <- readxl::read_excel("1-s2.0-S1535610822002744-mmc3.xlsx", sheet = 1)
matrix <- as.data.frame(matrix)
matrix2 <- matrix[, -1]
rownames(matrix2) <- matrix[, 1]
matrix3 <- matrix2[-1, ]
colnames(matrix3) <- matrix2[1, ]
matrix3[is.na(matrix3)] <- 0
matrix3 <- lapply(matrix3, as.numeric)
matrix3 <- as.data.frame(matrix3)

m3_cor <- cor(matrix3, method = "spearman")
m3_cor[is.na(m3_cor)] <- 0
###
m3_cor <- as.matrix(m3_cor)

# Subset correlations for CENPB/C and KI67/CPC

m3_cor_CENP_B_C <- m3_cor[, grep("CENPB|CENPC", colnames(m3_cor))]
m3_cor_CENP_B_C <- as.data.frame(m3_cor_CENP_B_C)
m3_cor_CENP_B_C$IDs <- rownames(m3_cor)
write_xlsx(m3_cor_CENP_B_C, "corr_matrix_CENP_B_C.xlsx")

m3_cor_CPC_KI67 <- m3_cor[, grep("AURKB|BIRC5|BOREA|INCE|KI67", colnames(m3_cor))]
m3_cor_CPC_KI67 <- as.data.frame(m3_cor_CPC_KI67)
m3_cor_CPC_KI67$IDs <- rownames(m3_cor)
write_xlsx(m3_cor_CPC_KI67, "corr_matrix_CPC_KI67.xlsx")

# Subset and plot correlation ranges for for CENPs and NDC80

m3_cor_CENPsNDC80cm <- m3_cor[, grep("CENP|NDC80|NUF2|SPC24|SPC25", colnames(m3_cor))]
m3_cor_CENPsNDC80cm <- as.data.frame(m3_cor_CENPsNDC80cm)

## Exclude CENP-F
m3_cor_CENPsNDC80cm <- m3_cor_CENPsNDC80cm[, -1]
m3_cor_CENPsNDC80cm <- as.matrix(m3_cor_CENPsNDC80cm)
row.names <- c(
  "P07199.CENPB_HUMAN", "Q03188.CENPC_HUMAN", "Q02224.CENPE_HUMAN",
  "Q9H3R5.CENPH_HUMAN", "Q9BS16.CENPK_HUMAN", "Q9NSP4.CENPM_HUMAN",
  "Q6IPU0.CENPP_HUMAN", "Q7L2Z9.CENPQ_HUMAN", "Q71F23.CENPU_HUMAN",
  "Q8N2Z9.CENPS_HUMAN", "A8MT69.CENPX_HUMAN", "Q7Z7K6.CENPV_HUMAN",
  "O14777.NDC80_HUMAN", "Q9BZD4.NUF2_HUMAN", "Q8NBT2.SPC24_HUMAN",
  "Q9HBM1.SPC25_HUMAN"
)
m3_cor_CENPsNDC80cm <- m3_cor_CENPsNDC80cm[, row.names]
m3_cor_CENPsNDC80cm <- as.data.frame(m3_cor_CENPsNDC80cm)
m3_cor_CENPsNDC80cm$IDs <- rownames(m3_cor)
write_xlsx(m3_cor_CENPsNDC80cm, "CENPs_NDC80_corrs_Goncalves.xlsx")

m3_cor_CENPsNDC80cm <- as.matrix(m3_cor_CENPsNDC80cm)
colnames(m3_cor_CENPsNDC80cm) <- gsub(".*\\.(.*)\\_.*", "\\1", colnames(m3_cor_CENPsNDC80cm))

m3_cor_CENPsNDC80cm_m <- melt(m3_cor_CENPsNDC80cm)
m3_cor_CENPsNDC80cm_m$value <- as.numeric(m3_cor_CENPsNDC80cm_m$value)

p_cd <- ggplot(m3_cor_CENPsNDC80cm_m, mapping = aes(x = `Var2`, y = `value`)) +
  geom_boxplot() +
  ylim(-1, 1) +
  xlab("Protein") +
  ylab("Correlation range") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf("CorrRangePlot_CENPs_Goncalves et all_new.pdf", width = 7, height = 4)
print(p_cd)
dev.off()

# Subset and plot protein levels ranges for CENPs and NDC80

m3_CENPsNDC80cm <- matrix3[, grep("CENP|NDC80|NUF2|SPC24|SPC25", colnames(matrix3))]
m3_CENPsNDC80cm <- as.data.frame(m3_CENPsNDC80cm)
m3_CENPsNDC80cm <- m3_CENPsNDC80cm[, -1]
m3_CENPsNDC80cm <- as.matrix(m3_CENPsNDC80cm)
row.names <- c(
  "P07199.CENPB_HUMAN", "Q03188.CENPC_HUMAN", "Q02224.CENPE_HUMAN",
  "Q9H3R5.CENPH_HUMAN", "Q9BS16.CENPK_HUMAN", "Q9NSP4.CENPM_HUMAN",
  "Q6IPU0.CENPP_HUMAN", "Q7L2Z9.CENPQ_HUMAN", "Q71F23.CENPU_HUMAN",
  "Q8N2Z9.CENPS_HUMAN", "A8MT69.CENPX_HUMAN", "Q7Z7K6.CENPV_HUMAN",
  "O14777.NDC80_HUMAN", "Q9BZD4.NUF2_HUMAN", "Q8NBT2.SPC24_HUMAN",
  "Q9HBM1.SPC25_HUMAN"
)
m3_CENPsNDC80cm <- m3_CENPsNDC80cm[, row.names]
colnames(m3_CENPsNDC80cm) <- gsub(".*\\.(.*)\\_.*", "\\1", colnames(m3_CENPsNDC80cm))

m3_CENPsNDC80cm_m <- melt(m3_CENPsNDC80cm)
m3_CENPsNDC80cm_m$value <- as.numeric(m3_CENPsNDC80cm_m$value)

p_cd <- ggplot(m3_CENPsNDC80cm_m, mapping = aes(x = `X2`, y = `value`)) +
  geom_boxplot() +
  xlab("Protein") +
  ylab("Protein levels range") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf("ProtLevelsRangePlot_CENPs_Goncalves et all_new.pdf", width = 7, height = 4)
print(p_cd)
dev.off()

# Subset and plot a correlation heatmap for CENPs and NDC80

m3_cor_CENPsNDC80 <- m3_cor[grep("CENP|NDC80|NUF2|SPC24|SPC25", rownames(m3_cor)), grep("CENP|NDC80|NUF2|SPC24|SPC25", colnames(m3_cor))]
m3_cor_CENPsNDC80 <- m3_cor_CENPsNDC80[-c(1), -c(1)]
row.names <- c(
  "P07199.CENPB_HUMAN", "Q03188.CENPC_HUMAN", "Q02224.CENPE_HUMAN",
  "Q9H3R5.CENPH_HUMAN", "Q9BS16.CENPK_HUMAN", "Q9NSP4.CENPM_HUMAN",
  "Q6IPU0.CENPP_HUMAN", "Q7L2Z9.CENPQ_HUMAN", "Q71F23.CENPU_HUMAN",
  "Q8N2Z9.CENPS_HUMAN", "A8MT69.CENPX_HUMAN", "Q7Z7K6.CENPV_HUMAN",
  "O14777.NDC80_HUMAN", "Q9BZD4.NUF2_HUMAN", "Q8NBT2.SPC24_HUMAN",
  "Q9HBM1.SPC25_HUMAN"
)
row.names <- rev(row.names)
m3_cor_CENPsNDC80 <- m3_cor_CENPsNDC80[row.names, row.names]
colnames(m3_cor_CENPsNDC80) <- gsub(".*\\.(.*)\\_.*", "\\1", colnames(m3_cor_CENPsNDC80))
rownames(m3_cor_CENPsNDC80) <- gsub(".*\\.(.*)\\_.*", "\\1", rownames(m3_cor_CENPsNDC80))

row.names2 <- gsub(".*\\.(.*)\\_.*", "\\1", row.names)

m3_cor_CENPsNDC80_melt <- melt.matrix(m3_cor_CENPsNDC80)
colnames(m3_cor_CENPsNDC80_melt) <- c("X1", "X2", "Correlation")
m3_cor_CENPsNDC80_melt$X1 <- factor(m3_cor_CENPsNDC80_melt$X1, levels = row.names2)
m3_cor_CENPsNDC80_melt$X2 <- factor(m3_cor_CENPsNDC80_melt$X2, levels = rev(row.names2))
p <- ggplot(m3_cor_CENPsNDC80_melt, aes(x = X2, y = X1, fill = Correlation)) +
  geom_tile() +
  coord_fixed() +
  scale_fill_gradientn(limits = c(-1, 1), colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("") +
  ylab("")
p
pdf("CorrHeatmap_Goncalves et al_CENPs_new.pdf", width = 9, height = 7)
print(p)
dev.off()

# Process Knol et al dataset and calculate correlations

setwd("/Users/natalia.kochanova/Desktop/Chromosome structure lab/Bioinformatics/R/202506_Correlation analysis_centromere/Other cancer datasets/")
library(readxl)
data_knol <- read_excel("mmc3_KnolEtAl2025.xlsx")
data_knol$ID <- paste(data_knol$`Protein Group`, data_knol$`UniProt Entry name`, sep = ".")
data_knol2 <- data_knol[, -c(1:10)]
rownames(data_knol2) <- data_knol$ID
data_knol3 <- data_knol2[, -1173]
data_knol3 <- as.matrix(data_knol3)
rownames(data_knol3) <- data_knol$ID
data_knol3 <- t(data_knol3)
class(data_knol3) <- "numeric"
data_knol3[is.na(data_knol3)] <- 0

dk3_cor <- cor(data_knol3, method = "spearman")
dk3_cor[is.na(dk3_cor)] <- 0
dk3_cor <- as.matrix(dk3_cor)


# Subset correlations for CENPA/B/C

dk3_cor_CENP_A_B_C <- dk3_cor[, grep("CENPA|CENPB|CENPC", colnames(dk3_cor))]
dk3_cor_CENP_A_B_C <- as.data.frame(dk3_cor_CENP_A_B_C)
dk3_cor_CENP_A_B_C$IDs <- rownames(dk3_cor)
write_xlsx(dk3_cor_CENP_A_B_C, "corr_matrix_knoll_CENP_A_B_C.xlsx")

# Subset correlations for Ki-67 and CPC

dk3_cor_CPC_KI67 <- dk3_cor[, grep("AURKB|BIRC5|BOREA|INCE|KI67", colnames(dk3_cor))]
dk3_cor_CPC_KI67 <- as.data.frame(dk3_cor_CPC_KI67)
dk3_cor_CPC_KI67$IDs <- rownames(dk3_cor)
write_xlsx(dk3_cor_CPC_KI67, "corr_matrix_CPC_KI67_knoll.xlsx")

# Subset and plot correlation ranges for for CENPs and NDC80

dk3_cor_CENPsNDC80cm <- dk3_cor[, grep("CENP|NDC80|NUF2|SPC24|SPC25", colnames(dk3_cor))]
dk3_cor_CENPsNDC80cm <- as.data.frame(dk3_cor_CENPsNDC80cm)
## Exclude CENP-F
dk3_cor_CENPsNDC80cm <- dk3_cor_CENPsNDC80cm[, -c(2, 13)]
dk3_cor_CENPsNDC80cm <- as.matrix(dk3_cor_CENPsNDC80cm)
row.names <- c(
  "P49450.CENPA_HUMAN", "P07199.CENPB_HUMAN", "Q03188.CENPC_HUMAN", "Q02224.CENPE_HUMAN",
  "Q9H3R5.CENPH_HUMAN", "Q92674.CENPI_HUMAN", "Q9BS16.CENPK_HUMAN", "Q9NSP4.CENPM_HUMAN",
  "Q8N0S6.CENPL_HUMAN", "Q96H22.CENPN_HUMAN",
  "Q9BU64.CENPO_HUMAN", "Q6IPU0.CENPP_HUMAN", "Q7L2Z9.CENPQ_HUMAN", "Q71F23.CENPU_HUMAN",
  "Q96BT3.CENPT_HUMAN", "Q8N2Z9.CENPS_HUMAN", "A8MT69.CENPX_HUMAN", "Q7Z7K6.CENPV_HUMAN",
  "O14777.NDC80_HUMAN", "Q9BZD4.NUF2_HUMAN", "Q8NBT2.SPC24_HUMAN", "Q9HBM1.SPC25_HUMAN"
)
dk3_cor_CENPsNDC80cm <- dk3_cor_CENPsNDC80cm[, row.names]
dk3_cor_CENPsNDC80cm <- as.data.frame(dk3_cor_CENPsNDC80cm)
dk3_cor_CENPsNDC80cm$IDs <- rownames(dk3_cor)
write_xlsx(dk3_cor_CENPsNDC80cm, "CENPs_NDC80_corrs_Knol.xlsx")

dk3_cor_CENPsNDC80cm <- as.matrix(dk3_cor_CENPsNDC80cm)
colnames(dk3_cor_CENPsNDC80cm) <- gsub(".*\\.(.*)\\_.*", "\\1", colnames(dk3_cor_CENPsNDC80cm))

dk3_cor_CENPsNDC80cm_m <- melt(dk3_cor_CENPsNDC80cm)
dk3_cor_CENPsNDC80cm_m$value <- as.numeric(dk3_cor_CENPsNDC80cm_m$value)

p_cd <- ggplot(dk3_cor_CENPsNDC80cm_m, mapping = aes(x = `X2`, y = `value`)) +
  geom_boxplot() +
  ylim(-1, 1) +
  xlab("Protein") +
  ylab("Correlation range") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf("CorrRangePlot_CENPs_Knol et all_new.pdf", width = 7, height = 4)
print(p_cd)
dev.off()

# Subset and plot protein levels ranges for CENPs and NDC80
dk3_CENPsNDC80cm <- data_knol3[, grep("CENP|NDC80|NUF2|SPC24|SPC25", colnames(data_knol3))]
dk3_CENPsNDC80cm <- as.data.frame(dk3_CENPsNDC80cm)
dk3_CENPsNDC80cm <- dk3_CENPsNDC80cm[, -c(2, 13)]
dk3_CENPsNDC80cm <- as.matrix(dk3_CENPsNDC80cm)

# Calculate how many zeros/NAs in the CENP columns
r <- colSums(dk3_CENPsNDC80cm == 0)
r <- unlist(r)
r2 <- r[-c(4, 7, 8)]
mean(r2)

## Impute data_knol3 with NAs not substituted ba zeros
dk3_CENPsNDC80cm <- data_knol3[, grep("CENP|NDC80|NUF2|SPC24|SPC25", colnames(data_knol3))]
dk3_CENPsNDC80cm <- dk3_CENPsNDC80cm[, -c(2, 13)]
dk3_CENPsNDC80cm <- dk3_CENPsNDC80cm[, -c(4, 7, 8)]
r <- sum(is.na(dk3_CENPsNDC80cm)) / ncol(dk3_CENPsNDC80cm)


row.names <- c(
  "P49450.CENPA_HUMAN", "P07199.CENPB_HUMAN", "Q03188.CENPC_HUMAN", "Q02224.CENPE_HUMAN",
  "Q9H3R5.CENPH_HUMAN", "Q92674.CENPI_HUMAN", "Q9BS16.CENPK_HUMAN", "Q9NSP4.CENPM_HUMAN",
  "Q8N0S6.CENPL_HUMAN", "Q96H22.CENPN_HUMAN",
  "Q9BU64.CENPO_HUMAN", "Q6IPU0.CENPP_HUMAN", "Q7L2Z9.CENPQ_HUMAN", "Q71F23.CENPU_HUMAN",
  "Q96BT3.CENPT_HUMAN", "Q8N2Z9.CENPS_HUMAN", "A8MT69.CENPX_HUMAN", "Q7Z7K6.CENPV_HUMAN",
  "O14777.NDC80_HUMAN", "Q9BZD4.NUF2_HUMAN", "Q8NBT2.SPC24_HUMAN", "Q9HBM1.SPC25_HUMAN"
)
dk3_CENPsNDC80cm <- dk3_CENPsNDC80cm[, row.names]
colnames(dk3_CENPsNDC80cm) <- gsub(".*\\.(.*)\\_.*", "\\1", colnames(dk3_CENPsNDC80cm))

dk3_CENPsNDC80cm_m <- melt(dk3_CENPsNDC80cm)
dk3_CENPsNDC80cm_m$value <- as.numeric(dk3_CENPsNDC80cm_m$value)

p_cd <- ggplot(dk3_CENPsNDC80cm_m, mapping = aes(x = `X2`, y = `value`)) +
  geom_boxplot() +
  xlab("Protein") +
  ylab("Protein levels range") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf("ProtLevelsRangePlot_CENPs_Knol et all_new.pdf", width = 7, height = 4)
print(p_cd)
dev.off()


# Subset and plot a correlation heatmap for CENPs and NDC80

dk3_cor_CENPsNDC80cm <- dk3_cor[, grep("CENP|NDC80|NUF2|SPC24|SPC25", colnames(dk3_cor))]
dk3_cor_CENPsNDC80cm <- as.data.frame(dk3_cor_CENPsNDC80cm)
dk3_cor_CENPsNDC80cm$IDs <- rownames(dk3_cor)
write_xlsx(dk3_cor_CENPsNDC80cm, "corr_matrix_knoll_CENPsNDC80.xlsx")

dk3_cor_CENPsNDC80 <- dk3_cor[grep("CENP|NDC80|NUF2|SPC24|SPC25", rownames(dk3_cor)), grep("CENP|NDC80|NUF2|SPC24|SPC25", colnames(dk3_cor))]
dk3_cor_CENPsNDC80 <- dk3_cor_CENPsNDC80[-c(2, 13), -c(2, 13)]
row.names <- c(
  "P49450.CENPA_HUMAN", "P07199.CENPB_HUMAN", "Q03188.CENPC_HUMAN", "Q02224.CENPE_HUMAN",
  "Q9H3R5.CENPH_HUMAN", "Q92674.CENPI_HUMAN", "Q9BS16.CENPK_HUMAN", "Q9NSP4.CENPM_HUMAN",
  "Q8N0S6.CENPL_HUMAN", "Q96H22.CENPN_HUMAN",
  "Q9BU64.CENPO_HUMAN", "Q6IPU0.CENPP_HUMAN", "Q7L2Z9.CENPQ_HUMAN", "Q71F23.CENPU_HUMAN",
  "Q96BT3.CENPT_HUMAN", "Q8N2Z9.CENPS_HUMAN", "A8MT69.CENPX_HUMAN", "Q7Z7K6.CENPV_HUMAN",
  "O14777.NDC80_HUMAN", "Q9BZD4.NUF2_HUMAN", "Q8NBT2.SPC24_HUMAN", "Q9HBM1.SPC25_HUMAN"
)
row.names <- rev(row.names)
dk3_cor_CENPsNDC80 <- dk3_cor_CENPsNDC80[row.names, row.names]
colnames(dk3_cor_CENPsNDC80) <- gsub(".*\\.(.*)\\_.*", "\\1", colnames(dk3_cor_CENPsNDC80))
rownames(dk3_cor_CENPsNDC80) <- gsub(".*\\.(.*)\\_.*", "\\1", rownames(dk3_cor_CENPsNDC80))
row.names2 <- gsub(".*\\.(.*)\\_.*", "\\1", row.names)

dk3_cor_CENPsNDC80_melt <- melt.matrix(dk3_cor_CENPsNDC80)
colnames(dk3_cor_CENPsNDC80_melt) <- c("X1", "X2", "Correlation")

dk3_cor_CENPsNDC80_melt$X1 <- factor(dk3_cor_CENPsNDC80_melt$X1, levels = row.names2)
dk3_cor_CENPsNDC80_melt$X2 <- factor(dk3_cor_CENPsNDC80_melt$X2, levels = rev(row.names2))

p <- ggplot(dk3_cor_CENPsNDC80_melt, aes(x = X2, y = X1, fill = Correlation)) +
  geom_tile() +
  coord_fixed() +
  scale_fill_gradientn(limits = c(-1, 1), colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("") +
  ylab("")
p
pdf("CorrHeatmap_CENPs_Knol et all_new.pdf", width = 10, height = 8)
print(p)
dev.off()
