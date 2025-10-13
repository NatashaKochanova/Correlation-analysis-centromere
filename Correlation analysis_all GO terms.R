setwd("/Users/natalia.kochanova/Desktop/Chromosome structure lab/Bioinformatics/R/202509_Correlation analysis_GO terms/")

library(readxl)
setwd("/datastore/home/nkochano/Desktop/")
setwd("Correlations_GO_SL_int/")

### Copied from Correlation-analysis-chromosome-instability https://github.com/NatashaKochanova/Chromosome-instability-correlation-analysis (#### to ####)
####

# Reading and preparing the proteomic data from Goncalves et al:
matrix <- readxl::read_excel("1-s2.0-S1535610822002744-mmc3.xlsx", sheet = 1)
matrix <- as.data.frame(matrix)
matrix2 <- matrix[, -1]
rownames(matrix2) <- matrix[, 1]
matrix3 <- matrix2[-1, ]
colnames(matrix3) <- matrix2[1, ]
matrix3[is.na(matrix3)] <- 0
matrix3 <- lapply(matrix3, as.numeric)
matrix3 <- as.data.frame(matrix3)
# Making a correlation matrix:
m3_cor <- cor(matrix3, method = "spearman")
m3_cor[is.na(m3_cor)] <- 0

# The same for proteomic data from Knol et al:
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

####
# Map GO terms for Goncalves et al (ids mapped previously in https://github.com/NatashaKochanova/Chromosome-instability-correlation-analysis)
GO_mapping <- read_xlsx("idmapping_active_true_2025_09_02.xlsx")
GOs <- lapply(GO_mapping$`Gene Ontology IDs`, function(x) strsplit(x, "GO:"))
GOs <- unique(unlist(GOs))
GOs <- gsub(";", "", GOs)
GOs <- gsub(" ", "", GOs)
GOs <- GOs[-c(1, 16)]

correlations_GO <- vector("list", length(GOs))
names(correlations_GO) <- GOs

for (i in 1:length(GOs)) {
  names(correlations_GO)[i] <- GOs[i]
  ids <- GO_mapping[grep(GOs[i], GO_mapping$`Gene Ontology IDs`), "From"]
  m3cor_ids <- m3_cor[unlist(lapply(ids$From, function(x) grep(x, rownames(m3_cor), fixed = TRUE))), unlist(lapply(ids$From, function(x) grep(x, colnames(m3_cor), fixed = TRUE))), drop = FALSE]
  if (length(ids) > 1) {
    m3cor_ids <- m3cor_ids[!duplicated(rownames(m3cor_ids)), !duplicated(colnames(m3cor_ids)), drop = FALSE]
  }
  m3cor_ids[upper.tri(m3cor_ids)] <- 1
  m3cor_ids <- as.vector(m3cor_ids)
  m3cor_ids <- m3cor_ids[m3cor_ids != 1]
  correlations_GO[[i]] <- median(m3cor_ids)
  print(i / length(GOs))
  i <- i + 1
}

correlations_GO <- unlist(correlations_GO)
attr(correlations_GO, "names") <- NULL
write.table(correlations_GO, "medians_GO.txt")
correlations_GO <- read.table("medians_GO.txt")

m3_cor_all <- m3_cor
m3_cor_all[upper.tri(m3_cor_all)] <- 1
m3_cor_all <- as.vector(m3_cor_all)
m3_cor_all <- m3_cor_all[m3_cor_all != 1]
med <- median(m3_cor_all)
length(unique(m3_cor_all))
correlations_GO <- as.numeric(correlations_GO)
correlations_GO <- as.data.frame(correlations_GO)
med_corrs <- median(correlations_GO$correlations_GO, na.rm = T)

ggr_plot <- ggplot(data = correlations_GO) +
  geom_density(aes(x = correlations_GO), alpha = 0.5, fill = c("skyblue")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Medians distribution") +
  geom_vline(xintercept = med_corrs, colour = "black") +
  geom_vline(xintercept = med, colour = "darkblue")

pdf("GO_medians_Goncalves_new.pdf", width = 15, height = 6)
print(ggr_plot)
dev.off()

# Map GO terms Knol et al
library(stringr)
IDs_Knol <- lapply(data_knol$`Protein Group`, function(x) str_split(x, ";"))
IDs.Knol.df <- data.frame(unlist(IDs_Knol), unlist(IDs_Knol))
writexl::write_xlsx(IDs.Knol.df, "IDs_Knol.xlsx")

GO_mapping_Knol <- read_excel("idmapping_2025_09_21_Knol_all.xlsx")
GOs <- lapply(GO_mapping_Knol$`Gene Ontology IDs`, function(x) strsplit(x, "GO:"))
GOs <- unique(unlist(GOs))
GOs <- gsub(";", "", GOs)
GOs <- gsub(" ", "", GOs)
GOs <- GOs[-c(1, 112)]

correlations_GO <- vector("list", length(GOs))
names(correlations_GO) <- GOs

for (i in 1:length(GOs)) {
  names(correlations_GO)[i] <- GOs[i]
  ids <- GO_mapping_Knol[grep(GOs[i], GO_mapping_Knol$`Gene Ontology IDs`), "From"]
  knol_ids <- dk3_cor[unlist(lapply(ids$From, function(x) grep(x, rownames(dk3_cor), fixed = TRUE))), unlist(lapply(ids$From, function(x) grep(x, colnames(dk3_cor), fixed = TRUE))), drop = FALSE]
  if (length(ids) > 1) {
    knol_ids <- knol_ids[!duplicated(rownames(knol_ids)), !duplicated(colnames(knol_ids)), drop = FALSE]
  }
  knol_ids[upper.tri(knol_ids)] <- 1
  knol_ids <- as.vector(knol_ids)
  knol_ids <- knol_ids[knol_ids != 1]
  correlations_GO[[i]] <- median(knol_ids)
  print(i / length(GOs))
  i <- i + 1
}

correlations_GO <- unlist(correlations_GO)
attr(correlations_GO, "names") <- NULL
write.table(correlations_GO, "medians_GO_Knol.txt")

m3_cor_all <- dk3_cor
m3_cor_all[upper.tri(m3_cor_all)] <- 1
m3_cor_all <- as.vector(m3_cor_all)
m3_cor_all <- m3_cor_all[m3_cor_all != 1]
med <- median(m3_cor_all)
length(unique(m3_cor_all))
correlations_GO <- as.numeric(correlations_GO)
correlations_GO <- as.data.frame(correlations_GO)
med_corrs <- median(correlations_GO$correlations_GO, na.rm = T)

ggr_plot <- ggplot(data = correlations_GO) +
  geom_density(aes(x = correlations_GO), alpha = 0.5, fill = c("skyblue")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Medians distribution") +
  geom_vline(xintercept = med_corrs, colour = "black") +
  geom_vline(xintercept = med, colour = "darkblue")

pdf("GO_medians_Knol.pdf", width = 15, height = 6)
print(ggr_plot)
dev.off()


#### Filtering synthetic lethal interactions
### Copied from Correlation-analysis-chromosome-instability - https://github.com/NatashaKochanova/Chromosome-instability-correlation-analysis (#### to ####)
####

# Synthetic lethality (SL) database
sl <- read.csv("Human_SL.csv")


# Filter CRISPR and low throughput screens interactions
sl_selected <- sl[sl$r.source == "CRISPR/CRISPRi" | sl$r.source == "Low Throughput" | sl$r.source == "High Throughput|Low Throughput", ]
sl_selected_1 <- sl_selected$n1.name
sl_selected_2 <- sl_selected$n2.name
sl_selected_ids <- c(sl_selected_1, sl_selected_2)
sl_selected_ids <- unique(sl_selected_ids)
length(sl_selected_ids)
c <- data.frame(sl_selected_ids, sl_selected_ids)
write_xlsx(c, "ids_sl_cr&lt.xlsx")

ids <- readxl::read_excel("idmapping_2024_02_06_cr_lt.xlsx")
ids <- ids[ids$Reviewed == "reviewed", ]

ids_sl_cr_lt <- c()
for (i in 1:nrow(ids)) {
  if (length(grep(ids[i, "Entry"], names(matrix3), value = TRUE)) >= 1) {
    ids_sl_cr_lt <- append(ids_sl_cr_lt, grep(ids[i, "Entry"], names(matrix3), value = TRUE))
  }
}
ids_sl_cr_lt <- unique(ids_sl_cr_lt)
length(ids_sl_cr_lt) # there are 677 of them

####

# Select SL pairs where both proteins are in the dataset:
sl_cr_lt_sel <- data.frame()
for (i in 1:nrow(sl_selected)) {
  if (length(grep(ids[ids$From == sl_selected$n1.name[i], "Entry"], names(matrix3), value = TRUE)) >= 1) {
    if (length(grep(ids[ids$From == sl_selected$n2.name[i], "Entry"], names(matrix3), value = TRUE)) >= 1) {
      sl_cr_lt_sel <- rbind(sl_cr_lt_sel, sl_selected[i, ])
    }
  }
}

nrow(sl_cr_lt_sel)



# Filter all SL correlations
correlations_sl <- data.frame()
for (i in 1:nrow(sl_cr_lt_sel)) {
  prot_1 <- grep(ids[ids$From == sl_cr_lt_sel$n1.name[i], "Entry"], names(matrix3), value = TRUE)
  prot_2 <- grep(ids[ids$From == sl_cr_lt_sel$n2.name[i], "Entry"], names(matrix3), value = TRUE)
  cor <- m3_cor[prot_1, prot_2]
  v <- c(prot_1, prot_2, cor)
  correlations_sl <- rbind(correlations_sl, v)
}

library(reshape2)
m3_cor_lt <- m3_cor
m3_cor_lt[upper.tri(m3_cor_lt)] <- NA
correlations_all <- melt(m3_cor_lt)
names(correlations_all) <- c("Protein_1", "Protein_2", "Correlation")
correlations_all <- correlations_all[!is.na(correlations_all$Correlation), ]
correlations_all <- correlations_all[!(correlations_all$Correlation == 1), ]
correlations_all$Distribution <- "All"
names(correlations_sl) <- c("Protein_1", "Protein_2", "Correlation")
correlations_sl$Distribution <- "Synthetic lethal"
correlations_all$Protein_1 <- as.character(correlations_all$Protein_1)
correlations_all$Protein_2 <- as.character(correlations_all$Protein_2)
correlations_sl$Correlation <- as.numeric(correlations_sl$Correlation)
library(data.table)
correlations_all_all <- rbindlist(list(correlations_all, correlations_sl))
library(ggplot2)
library(ggridges)
med <- median(correlations_sl$Correlation)
med2 <- median(correlations_all$Correlation)
GO_distr_all_plot_gd_all <- ggplot(correlations_all_all, aes(y = as.factor(Distribution), x = as.numeric(Correlation))) +
  geom_density_ridges(fill = "skyblue", alpha = 0.5) +
  ylab("Synthetic lethality") +
  xlab("Correlations distribution") +
  geom_vline(xintercept = med) +
  geom_vline(xintercept = med2)
GO_distr_all_plot_gd_all

pdf("GO_SL_Goncalves.pdf", width = 9, height = 7)
print(GO_distr_all_plot_gd_all)
dev.off()

## For Knol et al
####
ids <- readxl::read_excel("idmapping_2024_02_06_cr_lt.xlsx")
ids <- ids[ids$Reviewed == "reviewed", ]

ids_sl_cr_lt <- c()
for (i in 1:nrow(ids)) {
  if (length(grep(ids[i, "Entry"], colnames(data_knol3), value = TRUE)) >= 1) {
    ids_sl_cr_lt <- append(ids_sl_cr_lt, grep(ids[i, "Entry"], colnames(data_knol3), value = TRUE))
  }
}
ids_sl_cr_lt <- unique(ids_sl_cr_lt)
length(ids_sl_cr_lt) # there are 1125 of them
####

# Select SL pairs where both proteins are in the dataset:
sl_cr_lt_sel <- data.frame()
for (i in 1:nrow(sl_selected)) {
  if (length(grep(ids[ids$From == sl_selected$n1.name[i], "Entry"], colnames(data_knol3), value = TRUE)) >= 1) {
    if (length(grep(ids[ids$From == sl_selected$n2.name[i], "Entry"], colnames(data_knol3), value = TRUE)) >= 1) {
      sl_cr_lt_sel <- rbind(sl_cr_lt_sel, sl_selected[i, ])
    }
  }
}

nrow(sl_cr_lt_sel)



# Filter all SL correlations
library(reshape2)
correlations_sl <- data.frame()

for (i in 1:nrow(sl_cr_lt_sel)) {
  prot_1 <- grep(ids[ids$From == sl_cr_lt_sel$n1.name[i], "Entry"], colnames(data_knol3), value = TRUE)
  prot_2 <- grep(ids[ids$From == sl_cr_lt_sel$n2.name[i], "Entry"], colnames(data_knol3), value = TRUE)
  cor <- dk3_cor[prot_1, prot_2, drop = FALSE]
  if (length(cor) > 1) {
    v <- reshape2::melt(cor)
    names(v) <- c("Protein_1", "Protein_2", "Correlation")
  } else if (length(cor) == 1) {
    v <- c(prot_1, prot_2, cor)
    names(v) <- c("Protein_1", "Protein_2", "Correlation")
  }

  correlations_sl <- rbind(correlations_sl, v)
  names(correlations_sl) <- c("Protein_1", "Protein_2", "Correlation")
}
correlations_sl <- unique(correlations_sl)
correlations_sl <- correlations_sl[correlations_sl$Correlation != 1, ]

library(reshape2)
m3_cor_lt <- dk3_cor
m3_cor_lt[upper.tri(m3_cor_lt)] <- NA
correlations_all <- reshape2::melt(m3_cor_lt)
names(correlations_all) <- c("Protein_1", "Protein_2", "Correlation")
correlations_all <- correlations_all[!is.na(correlations_all$Correlation), ]
correlations_all <- correlations_all[!(correlations_all$Correlation == 1), ]
correlations_all$Distribution <- "All"
names(correlations_sl) <- c("Protein_1", "Protein_2", "Correlation")
correlations_sl$Distribution <- "Synthetic lethal"
correlations_all$Protein_1 <- as.character(correlations_all$Protein_1)
correlations_all$Protein_2 <- as.character(correlations_all$Protein_2)
correlations_sl$Correlation <- as.numeric(correlations_sl$Correlation)
library(data.table)
correlations_all_all <- rbindlist(list(correlations_all, correlations_sl))
library(ggplot2)
library(ggridges)
correlations_all_all$Correlation <- as.numeric(correlations_all_all$Correlation)
all.c <- as.vector(correlations_all_all[correlations_all_all$Distribution == "Synthetic lethal", "Correlation"])
all.c2 <- as.vector(correlations_all_all[correlations_all_all$Distribution == "All", "Correlation"])
med <- median(unlist(all.c))
med2 <- median(unlist(all.c2))
GO_distr_all_plot_gd_all <- ggplot(correlations_all_all, aes(y = as.factor(Distribution), x = as.numeric(Correlation))) +
  geom_density_ridges(fill = "skyblue", alpha = 0.5) +
  ylab("Synthetic lethality") +
  xlab("Correlations distribution") +
  geom_vline(xintercept = med) +
  geom_vline(xintercept = med2)
GO_distr_all_plot_gd_all
pdf("GO_SL_Knol_new.pdf", width = 9, height = 7)
print(GO_distr_all_plot_gd_all)
dev.off()
