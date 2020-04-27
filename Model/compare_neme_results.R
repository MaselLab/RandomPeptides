# Comparing our significantly increasing/decreasing to Neme et al.

# Load packages.
library(tidyverse)

# Load peptide data.
peptide.data <- read.table(file = "Scripts/RandomPeptides/Data/supplemental_table_1.tsv", header = T, stringsAsFactors = F)
peptide.data

# Standard errors on weights.
peptide.data$std.err <- sqrt(1/peptide.data$Weight.nb)
# Significantly increasing.
length(peptide.data[peptide.data$Fitness.nb - 1.96*peptide.data$std.err > 1,]$PeptideID)
# Signficantly decreasing.
length(peptide.data[peptide.data$Fitness.nb + 1.96*peptide.data$std.err < 1,]$PeptideID)
# Different than what Neme et al. found.
length(
  peptide.data[
    (peptide.data$Fitness.nb + 1.96*peptide.data$std.err > 1) & (peptide.data$fitness == 0),
    ]$PeptideID
)

# Need to load in the peptides Neme et al. found beneficial.
rafik.results <- read_tsv(file = "Data/rafik_results.tsv")
rafik.results
pep.increasing <- rafik.results[!is.na(rafik.results$`Up any`),]$ID
length(pep.increasing)

277 - length(peptide.data[(peptide.data$Fitness.nb - 1.96*peptide.data$std.err > 1) &
                            (peptide.data$PeptideID %in% pep.increasing),]$PeptideID)

# I suspect that the peptides that are no longer "increasing" for us fall into one of two groups:
# 1) peptides that have counts for time points 2-4 > time point 1, but where 2-4 might be decreasing, or
# 2) peptides that increase slowly, where time points 1-2 aren't that different, but either 3 or 4 is
# elevated.
# Is this the case?
count.data <- read_tsv(file = "Scripts/RandomPeptides/Data/supplemental_table_1.tsv")
count.data

count.data$d1.total <- count.data$`d1-r1` + count.data$`d1-r2` + count.data$`d1-r3` + count.data$`d1-r4` + count.data$`d1-r5`
count.data$d2.total <- count.data$`d2-r1` + count.data$`d2-r2` + count.data$`d2-r3` + count.data$`d2-r4` + count.data$`d2-r5`
count.data$d3.total <- count.data$`d3-r1` + count.data$`d3-r2` + count.data$`d3-r3` + count.data$`d3-r4` + count.data$`d3-r5`
count.data$d4.total <- count.data$`d4-r1` + count.data$`d4-r2` + count.data$`d4-r3` + count.data$`d4-r4` + count.data$`d4-r5`

count.data[, c("PeptideID", "d1.total", "d2.total", "d3.total", "d4.total")]

# Getting the list of genes that are in Neme et al.'s list of "increasing" peptides, but not in our list.
different.genes <- peptide.data[(peptide.data$Fitness.nb - 1.96*peptide.data$std.err < 1) &
                                  (peptide.data$PeptideID %in% pep.increasing),]$PeptideID
different.genes
count.data[count.data$PeptideID %in% different.genes, c("PeptideID", "d1.total", "d2.total", "d3.total", "d4.total")][11:20,]

# Checking how many clusters are increasing/decreasing.
peptides.maxweight <- 
  peptide.data %>%
  group_by(Cluster) %>%
  filter(Weight.nb == max(Weight.nb))
peptides.maxweight

# Increasing.
length(peptides.maxweight[peptides.maxweight$Fitness.nb - 1.96*peptides.maxweight$std.err > 1,]$PeptideID) # 127
# Decreasing.
length(peptides.maxweight[peptides.maxweight$Fitness.nb + 1.96*peptides.maxweight$std.err < 1,]$PeptideID) # 480
