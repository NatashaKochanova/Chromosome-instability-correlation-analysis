library(readxl)
library(writexl)
library(ggplot2)
library(ggridges)
library(data.table)
library(dplyr)
library(tidyr)
library(reshape2)

# Reading and preparing the proteomic data:
#
setwd("/Users/natalia.kochanova/Desktop/Chromosome structure lab/Bioinformatics/R/202312_Correlation analysis_chromosome instability/Final_code/")
#
setwd("/Users/natalia.kochanova/Desktop/Chromosome structure lab/Bioinformatics/R/202312_Correlation analysis_chromosome instability/")

matrix <- readxl::read_excel("1-s2.0-S1535610822002744-mmc3.xlsx", sheet = 1)
matrix <- as.data.frame(matrix)
matrix2 <- matrix[, -1]
rownames(matrix2) <- matrix[, 1]
matrix3 <- matrix2[-1, ]
colnames(matrix3) <- matrix2[1, ]
matrix3[is.na(matrix3)] <- 0
matrix3 <- lapply(matrix3, as.numeric)
matrix3 <- as.data.frame(matrix3)

# Extracting the Uniprot IDs from the matrix and assigning GO terms to them:

matrix_uniprot_IDs <- colnames(matrix3)
matrix_uniprot_IDs <- sub("\\..*", "", matrix_uniprot_IDs)
mUIDs <- data.frame(matrix_uniprot_IDs, matrix_uniprot_IDs)
write_xlsx(mUIDs, "matrix_uniprot_IDs.xlsx")

# Assigning GO terms through the Uniprot Website and reading the table into R:

IDs_GO <- read_excel("idmapping_2024_02_03_GO_matrix_IDs.xlsx")

# Making a correlation matrix:

m3_cor <- cor(matrix3, method = "spearman")
m3_cor[is.na(m3_cor)] <- 0


# Checking which chromosome instability genes are in the human cancer dataset:
# (the list of chromosome instability genes was converted to Uniprot IDs through the Uniprot website)

uniprot_ids <- readxl::read_excel("idmapping_2023_12_22.xlsx")
ids <- c()
for (i in 1:nrow(uniprot_ids)) {
  if (length(grep(uniprot_ids[i, "Entry"], names(matrix3), value = TRUE)) >= 1) {
    ids <- append(ids, grep(uniprot_ids[i, "Entry"], names(matrix3), value = TRUE))
  }
}
ids <- unique(ids)
length(ids) # there are 21 of them


## Subsetting GO IDs:

centromere_ids <- IDs_GO[grepl("centromeric|centromere", IDs_GO$`Gene Ontology (GO)`), ]
kinetochore_ids <- IDs_GO[grepl("kinetochore", IDs_GO$`Gene Ontology (GO)`), ]
cellcycle_ids <- IDs_GO[grepl("cell cycle", IDs_GO$`Gene Ontology (GO)`), ]
celldivision_ids <- IDs_GO[grepl("cell division", IDs_GO$`Gene Ontology (GO)`), ]


## Calculating distributions of correlations for chromosomal instability proteins against GO term proteins and not:

distributions_all <- list()

for (i in 1:length(ids)) {
  prot <- grep(ids[i], names(matrix3), value = TRUE)

  cor_prot_1 <- m3_cor[prot, ]
  cor_prot_1 <- as.data.frame(cor_prot_1)
  cor_prot_1$ID <- rownames(cor_prot_1)
  cor_prot_1 <- cor_prot_1[order(-cor_prot_1[, 1]), ]

  centromere_distr_df <- cor_prot_1
  centromere_distr_df$Distribution <- "Centromere"
  centromere_distr_df$GO <- "No"
  for (j in 1:nrow(centromere_distr_df)) {
    if (sub("\\..*", "", centromere_distr_df[j, "ID"]) %in% centromere_ids$From) {
      centromere_distr_df[j, "GO"] <- "Yes"
    }
  }
  centromere_distr_df$Protein <- prot
  distributions_all <- rbind(distributions_all, centromere_distr_df)

  kinetochore_distr_df <- cor_prot_1
  kinetochore_distr_df$Distribution <- "Kinetochore"
  kinetochore_distr_df$GO <- "No"
  for (j in 1:nrow(kinetochore_distr_df)) {
    if (sub("\\..*", "", kinetochore_distr_df[j, "ID"]) %in% kinetochore_ids$From) {
      kinetochore_distr_df[j, "GO"] <- "Yes"
    }
  }
  kinetochore_distr_df$Protein <- prot
  distributions_all <- rbind(distributions_all, kinetochore_distr_df)

  cell_cycle_distr_df <- cor_prot_1
  cell_cycle_distr_df$Distribution <- "Cell Cycle"
  cell_cycle_distr_df$GO <- "No"
  for (j in 1:nrow(cell_cycle_distr_df)) {
    if (sub("\\..*", "", cell_cycle_distr_df[j, "ID"]) %in% cellcycle_ids$From) {
      cell_cycle_distr_df[j, "GO"] <- "Yes"
    }
  }
  cell_cycle_distr_df$Protein <- prot
  distributions_all <- rbind(distributions_all, cell_cycle_distr_df)

  cell_division_distr_df <- cor_prot_1
  cell_division_distr_df$Distribution <- "Cell Division"
  cell_division_distr_df$GO <- "No"
  for (j in 1:nrow(cell_division_distr_df)) {
    if (sub("\\..*", "", cell_division_distr_df[j, "ID"]) %in% celldivision_ids$From) {
      cell_division_distr_df[j, "GO"] <- "Yes"
    }
  }
  cell_division_distr_df$Protein <- prot
  distributions_all <- rbind(distributions_all, cell_division_distr_df)

  i <- i + 1
}


## Plotting yes-no GO distributions for individual proteins (ranked by the mean difference in the yes-no means across GO terms):

distributions_all$Distribution_Protein <- paste(distributions_all$Distribution, distributions_all$Protein)

HMv <- aggregate(cor_prot_1 ~ Distribution_Protein + GO, data = distributions_all, FUN = mean)
HMv_1 <- HMv %>%
  group_by(Distribution_Protein) %>%
  do(m = (.[.$GO == "Yes", "cor_prot_1"]) - (.[.$GO == "No", "cor_prot_1"])) %>%
  summarise(Distribution_Protein, Mean = m)
HMv_2 <- extract(HMv_1, Distribution_Protein, into = c("Distribution", "Protein"), "(.*)\\s+([^ ]+)$")
colnames(HMv_2) <- c("Distribution", "Protein", "value")
HMv_2 <- data.frame(HMv_2$Distribution, HMv_2$Protein, HMv_2$value$cor_prot_1)
colnames(HMv_2) <- c("Distribution", "Protein", "value")
HMv_3 <- aggregate(value ~ Protein, data = HMv_2, FUN = sum)
HMv_3 <- HMv_3 %>% arrange(desc(-value))


distributions_all$Protein <- factor(distributions_all$Protein,
  levels = HMv_3$Protein
)

GO_distr_all_plot_gd_prot_yesno <- ggplot(distributions_all, aes(y = as.factor(Protein), x = as.numeric(cor_prot_1), fill = GO)) +
  geom_density_ridges(alpha = 0.5) +
  facet_grid(. ~ Distribution) +
  ylab("Protein") +
  xlab("Сorrelations distribution") +
  scale_fill_manual(values = c("#14C7BA", "#E64A00")) +
  geom_vline(xintercept = 0)


## Comprising a table with GO correlations for all proteins in the human cancer dataset:

## Preparing IDs_GO df:

IDs_GO$ID <- "a"
for (i in 1:nrow(IDs_GO)) {
  IDs_GO[i, "ID"] <- paste(IDs_GO[i, "From"], ".", IDs_GO[i, "Entry Name"], sep = "")
}

IDs_GO_2 <- IDs_GO[, c("Gene Ontology (GO)", "ID")]

## Preparing separate dfs for GO terms:

centromere_ids <- IDs_GO_2[grepl("centromeric|centromere", IDs_GO_2$`Gene Ontology (GO)`), ]
kinetochore_ids <- IDs_GO_2[grepl("kinetochore", IDs_GO_2$`Gene Ontology (GO)`), ]
cellcycle_ids <- IDs_GO_2[grepl("cell cycle", IDs_GO_2$`Gene Ontology (GO)`), ]
celldivision_ids <- IDs_GO_2[grepl("cell division", IDs_GO_2$`Gene Ontology (GO)`), ]

## Making a dataframe with all correlations:
m3_cor_lt <- m3_cor
m3_cor_lt[upper.tri(m3_cor_lt)] <- NA
correlations_all <- melt(m3_cor_lt)
names(correlations_all) <- c("Protein_1", "Protein_2", "Correlation")
correlations_all <- correlations_all[!is.na(correlations_all$Correlation), ]
correlations_all$Distribution <- "All"
correlations_all$Protein_1 <- as.character(correlations_all$Protein_1)
correlations_all$Protein_2 <- as.character(correlations_all$Protein_2)


## Making them datatable format for fast merging:

corr_all_DT <- setDT(correlations_all)
IDs_GO_centromere <- setDT(centromere_ids)
IDs_GO_kinetochore <- setDT(kinetochore_ids)
IDs_GO_cellcycle <- setDT(cellcycle_ids)
IDs_GO_celldivision <- setDT(celldivision_ids)

## Merging the GO and correlations tables fast:

merge_centr_1 <- merge(corr_all_DT, IDs_GO_centromere, by.x = "Protein_1", by.y = "ID", all.x = FALSE)
merge_centr_2 <- merge(corr_all_DT, IDs_GO_centromere, by.x = "Protein_2", by.y = "ID", all.x = FALSE)
centromere_GO <- rbindlist(list(merge_centr_1, merge_centr_2), use.names = TRUE)
centromere_GO <- unique(centromere_GO)
centromere_GO$Distribution <- "Centromere"
centromere_GO$Correlations <- "All"

merge_kin_1 <- merge(corr_all_DT, IDs_GO_kinetochore, by.x = "Protein_1", by.y = "ID", all.x = FALSE)
merge_kin_2 <- merge(corr_all_DT, IDs_GO_kinetochore, by.x = "Protein_2", by.y = "ID", all.x = FALSE)
kinetochore_GO <- rbindlist(list(merge_kin_1, merge_kin_2), use.names = TRUE)
kinetochore_GO <- unique(kinetochore_GO)
kinetochore_GO$Distribution <- "Kinetochore"
kinetochore_GO$Correlations <- "All"

merge_cellcycle_1 <- merge(corr_all_DT, IDs_GO_cellcycle, by.x = "Protein_1", by.y = "ID", all.x = FALSE)
merge_cellcycle_2 <- merge(corr_all_DT, IDs_GO_cellcycle, by.x = "Protein_2", by.y = "ID", all.x = FALSE)
cellcycle_GO <- rbindlist(list(merge_cellcycle_1, merge_cellcycle_2), use.names = TRUE)
cellcycle_GO <- unique(cellcycle_GO)
cellcycle_GO$Distribution <- "Cell Cycle"
cellcycle_GO$Correlations <- "All"

merge_celldivision_1 <- merge(corr_all_DT, IDs_GO_celldivision, by.x = "Protein_1", by.y = "ID", all.x = FALSE)
merge_celldivision_2 <- merge(corr_all_DT, IDs_GO_celldivision, by.x = "Protein_2", by.y = "ID", all.x = FALSE)
celldivision_GO <- rbindlist(list(merge_celldivision_1, merge_celldivision_2), use.names = TRUE)
celldivision_GO <- unique(celldivision_GO)
celldivision_GO$Distribution <- "Cell Division"
celldivision_GO$Correlations <- "All"

corr_GO_all <- rbindlist(list(centromere_GO, kinetochore_GO, cellcycle_GO, celldivision_GO), use.names = TRUE)

# Comprising a datatable, containing correlations against GO proteins both for chromosome instability proteins and for all proteins in the human cancer dataset:

distributions_all$Correlations <- "Chromosome instability genes"
distributions_all_yes <- distributions_all[distributions_all$GO == "Yes", ]
distributions_all_yes <- distributions_all_yes[, -6]
names(distributions_all_yes) <- c("Correlation", "Protein_2", "Distribution", "GO", "Protein_1", "Correlations")
names(corr_GO_all) <- c("Protein_1", "Protein_2", "Correlation", "Distribution", "GO", "Correlations")
corr_GO_all_all <- rbindlist(list(distributions_all_yes, corr_GO_all), use.names = TRUE)

## Plotting correlations distributions for chromosome instability proteins vs all proteins in the human cancer proteome:
## Calculating statistics:

GO_distr_GO_all_vs_inst <- ggplot(corr_GO_all_all, aes(y = as.factor(Distribution), x = as.numeric(Correlation), fill = as.factor(Correlations))) +
  geom_density_ridges(alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("GO term") +
  xlab("Correlations distribution") +
  scale_fill_manual(values = c("#14C7BA", "#E64A00")) +
  geom_vline(xintercept = 0)

ks.test(as.numeric(unlist(corr_GO_all_all[corr_GO_all_all$Distribution == "Centromere" & corr_GO_all_all$Correlations == "All", "Correlation"])), as.numeric(unlist(corr_GO_all_all[corr_GO_all_all$Distribution == "Centromere" & corr_GO_all_all$Correlations == "Chromosome instability genes", "Correlation"])))
ks.test(as.numeric(unlist(corr_GO_all_all[corr_GO_all_all$Distribution == "Kinetochore" & corr_GO_all_all$Correlations == "All", "Correlation"])), as.numeric(unlist(corr_GO_all_all[corr_GO_all_all$Distribution == "Kinetochore" & corr_GO_all_all$Correlations == "Chromosome instability genes", "Correlation"])))
ks.test(as.numeric(unlist(corr_GO_all_all[corr_GO_all_all$Distribution == "Cell Cycle" & corr_GO_all_all$Correlations == "All", "Correlation"])), as.numeric(unlist(corr_GO_all_all[corr_GO_all_all$Distribution == "Cell Cycle" & corr_GO_all_all$Correlations == "Chromosome instability genes", "Correlation"])))
ks.test(as.numeric(unlist(corr_GO_all_all[corr_GO_all_all$Distribution == "Cell Division" & corr_GO_all_all$Correlations == "All", "Correlation"])), as.numeric(unlist(corr_GO_all_all[corr_GO_all_all$Distribution == "Cell Division" & corr_GO_all_all$Correlations == "Chromosome instability genes", "Correlation"])))

## Calculating and plotting the fractions of correlations > 0.2:

corr_GO_all_all$Distribution_correlations <- paste(corr_GO_all_all$Distribution, corr_GO_all_all$Correlations)
corr_fraction_GOall <- corr_GO_all_all %>%
  group_by(Distribution_correlations) %>%
  do(f = nrow(.[.$Correlation > 0.2, ]) / nrow(.)) %>%
  summarise(Distribution_correlations, Fraction = f)
corr_fraction_GOall2 <- colsplit(corr_fraction_GOall$Distribution_correlations, " ", c("Distribution", "Correlations"))
corr_fraction_GOall <- cbind(corr_fraction_GOall, corr_fraction_GOall2)
corr_fraction_GOall[c(1:2), "Distribution"] <- "Cell Cycle"
corr_fraction_GOall[c(3:4), "Distribution"] <- "Cell Division"
corr_fraction_GOall[c(1, 3), "Correlations"] <- "All"
corr_fraction_GOall[c(2, 4), "Correlations"] <- "Chromosome instability genes"

p_all_fractions <- ggplot(corr_fraction_GOall, aes(y = Fraction, x = Distribution, fill = Correlations)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = c("#14C7BA", "#E64A00"))


### Distributions of correlations against proteins, to which centromeric, kinetochore, etc proteins are synthetic lethal:

setwd("/Users/natalia.kochanova/Desktop/Chromosome structure lab/Bioinformatics/R/202401_Corr_analysis_synthetic lethality/new_corr/")

### Reading the synthetic lethality (SL) database:

sl <- read.csv("Human_SL.csv")

### Filtering CRISPR and low throughput screens interactions:

sl_selected <- sl[sl$r.source == "CRISPR/CRISPRi" | sl$r.source == "Low Throughput" | sl$r.source == "High Throughput|Low Throughput", ]
sl_selected_1 <- sl_selected$n1.name
sl_selected_2 <- sl_selected$n2.name
sl_selected_ids <- c(sl_selected_1, sl_selected_2)
sl_selected_ids <- unique(sl_selected_ids)
length(sl_selected_ids)
c <- data.frame(sl_selected_ids, sl_selected_ids)
write_xlsx(c, "ids_sl_cr&lt.xlsx")

ids_sl <- readxl::read_excel("idmapping_2024_02_06_cr_lt.xlsx")
ids_sl <- ids_sl[ids_sl$Reviewed == "reviewed", ]

ids_sl_cr_lt <- c()
for (i in 1:nrow(ids_sl)) {
  if (length(grep(ids_sl[i, "Entry"], names(matrix3), value = TRUE)) >= 1) {
    ids_sl_cr_lt <- append(ids_sl_cr_lt, grep(ids_sl[i, "Entry"], names(matrix3), value = TRUE))
  }
}
ids_sl_cr_lt <- unique(ids_sl_cr_lt)
length(ids_sl_cr_lt) # there are 677 of them

### Putting uniprot IDs in the SL table:

sl_selected$n1.name.ID <- 1
sl_selected$n2.name.ID <- 2
for (i in 1:nrow(sl_selected)) {
  prot <- sl_selected[i, "n1.name"]
  sl_selected[i, "n1.name.ID"] <- as.character(ids_sl[ids_sl$From == prot, "Entry"])
  prot2 <- sl_selected[i, "n2.name"]
  sl_selected[i, "n2.name.ID"] <- as.character(ids_sl[ids_sl$From == prot2, "Entry"])
}

centromere_ids <- IDs_GO[grepl("centromeric|centromere", IDs_GO$`Gene Ontology (GO)`), ]
kinetochore_ids <- IDs_GO[grepl("kinetochore", IDs_GO$`Gene Ontology (GO)`), ]
cellcycle_ids <- IDs_GO[grepl("cell cycle", IDs_GO$`Gene Ontology (GO)`), ]
celldivision_ids <- IDs_GO[grepl("cell division", IDs_GO$`Gene Ontology (GO)`), ]


### Filtering correlations of chromosome instability proteins against proteins, synthetic lethal with GO terms:

distributions_all$Synthetic_lethal <- "No"

centromere_sl <- distributions_all[distributions_all$Distribution == "Centromere", ]
for (i in 1:nrow(centromere_sl)) {
  prot <- sub("\\..*", "", centromere_sl[i, "ID"])
  if (prot %in% sl_selected$n1.name.ID) {
    n1 <- sl_selected[sl_selected$n1.name.ID == prot, ]
    for (j in 1:nrow(n1)) {
      if (n1[j, "n2.name.ID"] %in% centromere_ids$From) {
        centromere_sl[i, "Synthetic_lethal"] <- "Yes"
      }
    }
    if (prot %in% sl_selected$n2.name.ID) {
      n2 <- sl_selected[sl_selected$n2.name.ID == prot, ]
      for (j in 1:nrow(n2)) {
        if (n2[j, "n1.name.ID"] %in% centromere_ids$From) {
          centromere_sl[i, "Synthetic_lethal"] <- "Yes"
        }
      }
    }
  }
}
nrow(centromere_sl[centromere_sl$Synthetic_lethal == "Yes", ])

kinetochore_sl <- distributions_all[distributions_all$Distribution == "Kinetochore", ]
for (i in 1:nrow(kinetochore_sl)) {
  prot <- sub("\\..*", "", kinetochore_sl[i, "ID"])
  if (prot %in% sl_selected$n1.name.ID) {
    n1 <- sl_selected[sl_selected$n1.name.ID == prot, ]
    for (j in 1:nrow(n1)) {
      if (n1[j, "n2.name.ID"] %in% kinetochore_ids$From) {
        kinetochore_sl[i, "Synthetic_lethal"] <- "Yes"
      }
    }
    if (prot %in% sl_selected$n2.name.ID) {
      n2 <- sl_selected[sl_selected$n2.name.ID == prot, ]
      for (j in 1:nrow(n2)) {
        if (n2[j, "n1.name.ID"] %in% kinetochore_ids$From) {
          kinetochore_sl[i, "Synthetic_lethal"] <- "Yes"
        }
      }
    }
  }
}
nrow(kinetochore_sl[kinetochore_sl$Synthetic_lethal == "Yes", ])

cell_cycle_sl <- distributions_all[distributions_all$Distribution == "Cell Cycle", ]
for (i in 1:nrow(cell_cycle_sl)) {
  prot <- sub("\\..*", "", cell_cycle_sl[i, "ID"])
  if (prot %in% sl_selected$n1.name.ID) {
    n1 <- sl_selected[sl_selected$n1.name.ID == prot, ]
    for (j in 1:nrow(n1)) {
      if (n1[j, "n2.name.ID"] %in% cellcycle_ids$From) {
        cell_cycle_sl[i, "Synthetic_lethal"] <- "Yes"
      }
    }
    if (prot %in% sl_selected$n2.name.ID) {
      n2 <- sl_selected[sl_selected$n2.name.ID == prot, ]
      for (j in 1:nrow(n2)) {
        if (n2[j, "n1.name.ID"] %in% cellcycle_ids$From) {
          cell_cycle_sl[i, "Synthetic_lethal"] <- "Yes"
        }
      }
    }
  }
}
nrow(cell_cycle_sl[cell_cycle_sl$Synthetic_lethal == "Yes", ])

cell_division_sl <- distributions_all[distributions_all$Distribution == "Cell Division", ]
for (i in 1:nrow(cell_division_sl)) {
  prot <- sub("\\..*", "", cell_division_sl[i, "ID"])
  if (prot %in% sl_selected$n1.name.ID) {
    n1 <- sl_selected[sl_selected$n1.name.ID == prot, ]
    for (j in 1:nrow(n1)) {
      if (n1[j, "n2.name.ID"] %in% celldivision_ids$From) {
        cell_division_sl[i, "Synthetic_lethal"] <- "Yes"
      }
    }
    if (prot %in% sl_selected$n2.name.ID) {
      n2 <- sl_selected[sl_selected$n2.name.ID == prot, ]
      for (j in 1:nrow(n2)) {
        if (n2[j, "n1.name.ID"] %in% celldivision_ids$From) {
          cell_division_sl[i, "Synthetic_lethal"] <- "Yes"
        }
      }
    }
  }
}
nrow(cell_division_sl[cell_division_sl$Synthetic_lethal == "Yes", ])

distributions_all_sl <- centromere_sl
distributions_all_sl <- rbind(distributions_all_sl, kinetochore_sl)
distributions_all_sl <- rbind(distributions_all_sl, cell_cycle_sl)
distributions_all_sl <- rbind(distributions_all_sl, cell_division_sl)


### The same, but for all proteins in the cancer proteome dataset:

### Making a table containing IDs from human cancer proteome data and synthetic lethality with GO terms:

rn_IDs <- rownames(m3_cor)
rn_IDs <- cbind(rn_IDs, rn_IDs, rn_IDs, rn_IDs, rn_IDs)
colnames(rn_IDs) <- c("ID", "SL_Centr", "SL_Kin", "SL_CC", "SL_CD")
rn_IDs <- as.data.frame(rn_IDs)
rn_IDs$SL_Centr <- "No"
rn_IDs$SL_Kin <- "No"
rn_IDs$SL_CC <- "No"
rn_IDs$SL_CD <- "No"


for (i in 1:nrow(rn_IDs)) {
  prot <- sub("\\..*", "", rn_IDs[i, "ID"])
  if (prot %in% sl_selected$n1.name.ID) {
    n1 <- sl_selected[sl_selected$n1.name.ID == prot, ]
    for (j in 1:nrow(n1)) {
      if (n1[j, "n2.name.ID"] %in% centromere_ids$From) {
        rn_IDs[i, "SL_Centr"] <- "Yes"
      }
    }
    if (prot %in% sl_selected$n2.name.ID) {
      n2 <- sl_selected[sl_selected$n2.name.ID == prot, ]
      for (j in 1:nrow(n2)) {
        if (n2[j, "n1.name.ID"] %in% centromere_ids$From) {
          rn_IDs[i, "SL_Centr"] <- "Yes"
        }
      }
    }
  }
}
nrow(rn_IDs[rn_IDs$SL_Centr == "Yes", ])

for (i in 1:nrow(rn_IDs)) {
  prot <- sub("\\..*", "", rn_IDs[i, "ID"])
  if (prot %in% sl_selected$n1.name.ID) {
    n1 <- sl_selected[sl_selected$n1.name.ID == prot, ]
    for (j in 1:nrow(n1)) {
      if (n1[j, "n2.name.ID"] %in% kinetochore_ids$From) {
        rn_IDs[i, "SL_Kin"] <- "Yes"
      }
    }
    if (prot %in% sl_selected$n2.name.ID) {
      n2 <- sl_selected[sl_selected$n2.name.ID == prot, ]
      for (j in 1:nrow(n2)) {
        if (n2[j, "n1.name.ID"] %in% kinetochore_ids$From) {
          rn_IDs[i, "SL_Kin"] <- "Yes"
        }
      }
    }
  }
}
nrow(rn_IDs[rn_IDs$SL_Kin == "Yes", ])

for (i in 1:nrow(rn_IDs)) {
  prot <- sub("\\..*", "", rn_IDs[i, "ID"])
  if (prot %in% sl_selected$n1.name.ID) {
    n1 <- sl_selected[sl_selected$n1.name.ID == prot, ]
    for (j in 1:nrow(n1)) {
      if (n1[j, "n2.name.ID"] %in% celldivision_ids$From) {
        rn_IDs[i, "SL_CD"] <- "Yes"
      }
    }
    if (prot %in% sl_selected$n2.name.ID) {
      n2 <- sl_selected[sl_selected$n2.name.ID == prot, ]
      for (j in 1:nrow(n2)) {
        if (n2[j, "n1.name.ID"] %in% celldivision_ids$From) {
          rn_IDs[i, "SL_CD"] <- "Yes"
        }
      }
    }
  }
}
nrow(rn_IDs[rn_IDs$SL_CD == "Yes", ])

for (i in 1:nrow(rn_IDs)) {
  prot <- sub("\\..*", "", rn_IDs[i, "ID"])
  if (prot %in% sl_selected$n1.name.ID) {
    n1 <- sl_selected[sl_selected$n1.name.ID == prot, ]
    for (j in 1:nrow(n1)) {
      if (n1[j, "n2.name.ID"] %in% cellcycle_ids$From) {
        rn_IDs[i, "SL_CC"] <- "Yes"
      }
    }
    if (prot %in% sl_selected$n2.name.ID) {
      n2 <- sl_selected[sl_selected$n2.name.ID == prot, ]
      for (j in 1:nrow(n2)) {
        if (n2[j, "n1.name.ID"] %in% cellcycle_ids$From) {
          rn_IDs[i, "SL_CC"] <- "Yes"
        }
      }
    }
  }
}
nrow(rn_IDs[rn_IDs$SL_CC == "Yes", ])

rn_IDs <- setDT(rn_IDs)

### Merging all correlations with synthetic-lethal-with-GO table:

corr_all_DT_SL <- merge(corr_all_DT, rn_IDs, by.x = "Protein_1", by.y = "ID", all.x = TRUE)
corr_all_DT_SL <- merge(corr_all_DT_SL, rn_IDs, by.x = "Protein_2", by.y = "ID", all.x = TRUE)

corr_all_DT_SL[, SL_Centr := paste0(SL_Centr.x, SL_Centr.y)]
corr_all_DT_SL[, SL_Kin := paste0(SL_Kin.x, SL_Kin.y)]
corr_all_DT_SL[, SL_CC := paste0(SL_CC.x, SL_CC.y)]
corr_all_DT_SL[, SL_CD := paste0(SL_CD.x, SL_CD.y)]

corr_all_DT_SL_centr <- corr_all_DT_SL[grep("Yes", corr_all_DT_SL$SL_Centr), ]
corr_all_DT_SL_kin <- corr_all_DT_SL[grep("Yes", corr_all_DT_SL$SL_Kin), ]
corr_all_DT_SL_CC <- corr_all_DT_SL[grep("Yes", corr_all_DT_SL$SL_CC), ]
corr_all_DT_SL_CD <- corr_all_DT_SL[grep("Yes", corr_all_DT_SL$SL_CD), ]


corr_all_DT_SL_centr$Distribution <- "Centromere"
corr_all_DT_SL_kin$Distribution <- "Kinetochore"
corr_all_DT_SL_CC$Distribution <- "Cell Cycle"
corr_all_DT_SL_CD$Distribution <- "Cell Division"

corr_all_DT_SL_yes <- rbindlist(list(corr_all_DT_SL_centr, corr_all_DT_SL_kin, corr_all_DT_SL_CC, corr_all_DT_SL_CD))
corr_all_DT_SL_yes <- corr_all_DT_SL_yes[, c(1:4)]
corr_all_DT_SL_yes$Correlations <- "All"

### Making the datatable, containing correlations of both chromosome instability proteins and all proteins from the cancer human proteome against proteins, synthetic lethal with GOs:

distributions_all_sl_chinst <- distributions_all_sl
distributions_all_sl_chinst_yes <- distributions_all_sl_chinst[distributions_all_sl_chinst$Synthetic_lethal == "Yes", ]
distributions_all_sl_chinst_yes <- distributions_all_sl_chinst_yes[, c(1, 2, 3, 5, 6)]
names(distributions_all_sl_chinst_yes) <- c("Correlation", "Protein_2", "Distribution", "Protein_1", "Correlations")
distributions_all_sl_chinst_yes$Correlations <- "Chromosome instability genes"
distributions_all_sl_chinst_yes <- setDT(distributions_all_sl_chinst_yes)

corr_SL_all <- rbindlist(list(corr_all_DT_SL_yes, distributions_all_sl_chinst_yes), use.names = T)

rm(corr_all_DT)
rm(corr_all_DT_SL)
rm(correlations_all)

### Making the plots and calculating statistics:

GO_distr_plot_gd_all_vs_chrinst <- ggplot(corr_SL_all, aes(y = as.factor(Distribution), x = as.numeric(Correlation), fill = Correlations)) +
  geom_density_ridges(alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("GO term") +
  xlab("Correlations distribution") +
  scale_fill_manual(values = c("#14C7BA", "#E64A00")) +
  geom_vline(xintercept = 0)
ks.test(as.numeric(unlist(corr_SL_all[corr_SL_all$Distribution == "Centromere" & corr_SL_all$Correlations == "All", "Correlation"])), as.numeric(unlist(corr_SL_all[corr_SL_all$Distribution == "Centromere" & corr_SL_all$Correlations == "Chromosome instability genes", "Correlation"])))
ks.test(as.numeric(unlist(corr_SL_all[corr_SL_all$Distribution == "Kinetochore" & corr_SL_all$Correlations == "All", "Correlation"])), as.numeric(unlist(corr_SL_all[corr_SL_all$Distribution == "Kinetochore" & corr_SL_all$Correlations == "Chromosome instability genes", "Correlation"])))
ks.test(as.numeric(unlist(corr_SL_all[corr_SL_all$Distribution == "Cell Cycle" & corr_SL_all$Correlations == "All", "Correlation"])), as.numeric(unlist(corr_SL_all[corr_SL_all$Distribution == "Cell Cycle" & corr_SL_all$Correlations == "Chromosome instability genes", "Correlation"])))
ks.test(as.numeric(unlist(corr_SL_all[corr_SL_all$Distribution == "Cell Division" & corr_SL_all$Correlations == "All", "Correlation"])), as.numeric(unlist(corr_SL_all[corr_SL_all$Distribution == "Cell Division" & corr_SL_all$Correlations == "Chromosome instability genes", "Correlation"])))

corr_SL_all$Distribution_correlations <- paste(corr_SL_all$Distribution, corr_SL_all$Correlations)
corr_fraction_GOSLall <- corr_SL_all %>%
  group_by(Distribution_correlations) %>%
  do(f = nrow(.[.$Correlation > 0.2, ]) / nrow(.)) %>%
  summarise(Distribution_correlations, Fraction = f)
corr_fraction_GOSLall2 <- colsplit(corr_fraction_GOSLall$Distribution_correlations, " ", c("Distribution", "Correlations"))
corr_fraction_GOSLall <- cbind(corr_fraction_GOSLall, corr_fraction_GOSLall2)
corr_fraction_GOSLall[c(1:2), "Distribution"] <- "Cell Cycle"
corr_fraction_GOSLall[c(3:4), "Distribution"] <- "Cell Division"
corr_fraction_GOSLall[c(1, 3), "Correlations"] <- "All"
corr_fraction_GOSLall[c(2, 4), "Correlations"] <- "Chromosome instability genes"

p_all_fractions_SL <- ggplot(corr_fraction_GOSLall, aes(y = Fraction, x = Distribution, fill = Correlations)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = c("#14C7BA", "#E64A00"))

### For individual chromosome instability proteins - all their correlations (synthetic lethal yes or no). Ranking proteins by mean difference in the means across GO terms:

distributions_all_sl$Distribution_Protein <- paste(distributions_all_sl$Distribution, distributions_all_sl$Protein)

HMv <- aggregate(cor_prot_1 ~ Distribution_Protein + Synthetic_lethal, data = distributions_all_sl, FUN = mean)
HMv_1 <- HMv %>%
  group_by(Distribution_Protein) %>%
  do(m = (.[.$Synthetic_lethal == "Yes", "cor_prot_1"]) - (.[.$Synthetic_lethal == "No", "cor_prot_1"])) %>%
  summarise(Distribution_Protein, Mean = m)
HMv_2 <- extract(HMv_1, Distribution_Protein, into = c("Distribution", "Protein"), "(.*)\\s+([^ ]+)$")
colnames(HMv_2) <- c("Distribution", "Protein", "value")
HMv_2 <- data.frame(HMv_2$Distribution, HMv_2$Protein, HMv_2$value$cor_prot_1)
colnames(HMv_2) <- c("Distribution", "Protein", "value")
HMv_3 <- aggregate(value ~ Protein, data = HMv_2, FUN = sum)
HMv_3 <- HMv_3 %>% arrange(desc(-value))

distributions_all_sl$Protein <- factor(distributions_all_sl$Protein,
  levels = HMv_3$Protein
)
GO_distr_all_plot_gd_prot_yesno_sl <- ggplot(distributions_all_sl, aes(y = as.factor(Protein), x = as.numeric(cor_prot_1), fill = Synthetic_lethal)) +
  geom_density_ridges(alpha = 0.5) +
  facet_grid(. ~ Distribution) +
  ylab("Protein") +
  xlab("Сorrelations distribution") +
  scale_fill_manual(values = c("#14C7BA", "#E64A00")) +
  geom_vline(xintercept = 0)
