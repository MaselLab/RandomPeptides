# Supplementary table.

# Load packages.
library(tidyverse)

# Load full data.
peptide.data <- read.table(file = "Data/peptide_data_clusters_4-9-20.tsv", header = T, stringsAsFactors = F)

# Trimmed data table.
trimmed.data <- peptide.data[, c(
  "PeptideID", "AASeq", "Cluster", "Fitness.nb", "Weight.nb", "PredHel.x", "GC.avg",
  "ISD.iupred2", "ISD.delta", "Clustering.Six", "net.charge", "charge.pos", "charge.neg", "WaltzBinary", "CamSol.avg",
  "Ala", "Leu", "Cys", "Glu", "Gly", "Gln", "Phe", "Met", "Ile", "Val",
  "Pro", "Asp", "Asn", "Thr", "Tyr", "Trp", "Lys", "His", "Arg", "Ser"
)]
names(trimmed.data)[6]
names(trimmed.data)[6] <- "PredHel"

# Today's date.
#todays.date <- "3-27-2020"

# Getting sequencing counts.
count.data <- read_csv(file = "Data/Results_Exp7.csv")
count.data

count.data.trimmed <- count.data[count.data$PeptideID %in% peptide.data$PeptideID,]
count.data.trimmed

# Full data with counts.
full.data <- merge(trimmed.data, count.data.trimmed, by = "PeptideID")
full.data

#output.path <- paste("Scripts/RandomPeptides/Data/peptide_data_counts_", todays.date, ".tsv", sep = "")
output.path <- paste("Scripts/RandomPeptides/Data/supplemental_table_1.tsv", sep = "")
write_tsv(full.data, path = output.path)
