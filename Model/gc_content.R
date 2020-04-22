# Figuring out the GC content of the 1000 random peptides and the full million.

# Load packages.
library(tidyverse)
library(Biostrings)
library(wCorr)

# Importing random peptide data.
peptide.data <- read.table(file = "Data/peptide_data_clusters_2-14-20.tsv", header = T, stringsAsFactors = F)

# Variable for checking just the random region. If T, only the random region is checked for GC content. If F,
# the whole sequence is checked for GC content.
random.region <- T

# Importing nucleotide sequences and matching them up with the appropriate peptide ID.
nt.fasta <- readDNAStringSet("Data/ORF_NR99_DB.fa")
translate(nt.fasta[1])[[1]]
nt.seqs.df <- tibble("id" = names(nt.fasta),
                     "nts" = paste(nt.fasta),
                     "aas" = paste(translate(nt.fasta, if.fuzzy.codon = "solve")))
# The following command returned multiple matches.
nt.seqs.df[as.character(nt.seqs.df$aas) == as.character(peptide.data$AASeq[1]), "nts"]
# This means I can't actually figure out the exact nucleotide sequence of the peptides I have.
# A consensus will have to suffice.
length(unique(nt.seqs.df$aas))
length(nt.seqs.df$aas)
length(unique(nt.seqs.df$nts))
#merged.sequences.df <- merge(nt.seqs.df, peptide.data[, c("PeptideID", "AASeq")])

# Working out the substrings to get just random region.
nt.seqs.df$aas[1:20]
nt.fasta[1]
substr(nt.fasta[1], 1, 195)
substr(nt.fasta[1], 1 + 12, 195 - 33)
str_length(substr(nt.fasta[1], 13, 162))
translate(DNAString(substr(nt.fasta[1], 13, 162))) #This looks like the random region.
paste(substr(nt.fasta[1:10], 13, 162))

# Examining GC content of just the random region if random.region = T. If random.region = F, then
# this code is skipped and the full RNA sequence, including non-random regions, is checked for GC
# content. NOTE: The translated sequences in the tibble still include non-random regions.
if (random.region == T){
  nt.seqs.df <- tibble("id" = names(nt.fasta),
                       "nts" = paste(substr(nt.fasta, 13, 162)),
                       "aas" = paste(translate(nt.fasta, if.fuzzy.codon = "solve")))
}
if (random.region == T){
  nt.fasta <- DNAStringSet(nt.seqs.df$nts)
}

# Working out GC content of each nucleotide sequence.
alphabetFrequency(DNAString(nt.seqs.df$nts[1]))[c("G", "C", "A", "T")]

# Calculate mean GC content for a set of sequences.
mean.gc <- function(nts.set){
  require(Biostrings)
  nts.counts <- alphabetFrequency(nts.set)[,c("G", "C", "A", "T")]
  nts.total <- sum(nts.counts)
  gc.total <- NA
  # If nts.counts has only 1 row, it doesn't get stored as a matrix, and [,c()] notation
  # will return an error. So, to deal with the possibility that nts.set has only one sequence,
  # we check first if nrow(nts.counts) is null.
  if (is.null(nrow(nts.counts))){
    gc.total <- sum(nts.counts[c("G", "C")])
  } else {
    gc.total <- sum(nts.counts[,c("G", "C")])
  }
  return(gc.total / nts.total)
}

# Calculating mean GC for each peptide in our dataset (but not the full set of 2 million + sequences).
peptide.data$GC.avg <- NA
temp.seq.set <- NA
for (i in 1:length(peptide.data$PeptideID)) {
  temp.seq.set <- nt.seqs.df[as.character(nt.seqs.df$aas) == as.character(peptide.data$AASeq[i]), "nts"]
  peptide.data$GC.avg[i] <- mean.gc(DNAStringSet(temp.seq.set$nts))
}
head(peptide.data$GC.avg)
tail(peptide.data$GC.avg)

# Figuring out why some amino acid sequences weren't found in the supposedly "full" list from Neme et al.
# Importing more data. First, the full list of peptides.
full.aa.data <- read_tsv(file = "Data/PEP_NR_DB.tab", col_names = F)
names(full.aa.data) <- c("PeptideID", "AASeq")
full.aa.data

# Checking the full list of peptide sequences.
peptide.data[!peptide.data$AASeq %in% full.aa.data$AASeq, c("PeptideID", "AASeq", "Cluster")]
# Nothing is missing.

# Checking the full list of nucleotide sequences.
nt.seqs.df[!nt.seqs.df$aas %in% full.aa.data$AASeq,]
full.aa.data[!full.aa.data$AASeq %in% nt.seqs.df$aas,] # There are 44796 nucleotide seqs missing.
length(full.aa.data$PeptideID) - length(unique(nt.seqs.df$aas)) # Same number missing. OK.

# Checking Rafik's files with nucleotide sequences for the significantly increasing/decreasing. Maybe
# a few of the missing nucleotide sequences are in here.
aa.nt.set1 <- read_csv(file = "Data/peptide_aa_nt_seqs.csv")
aa.nt.set2 <- read_csv(file = "Data/peptide_aa_nt_seqs_pt2.csv")

aa.nt.set1
aa.nt.set2

aa.nt.set1[!aa.nt.set1$PeptideID %in% aa.nt.set2$PeptideID,]
aa.nt.set1[!aa.nt.set1$AASeq %in% aa.nt.set2$AASeq,]
aa.nt.set1[!aa.nt.set1$NTSeq %in% aa.nt.set2$NTSeq,]
length(peptide.data[!peptide.data$AASeq %in% aa.nt.set1$AASeq, "PeptideID"])
length(peptide.data[!peptide.data$AASeq %in% aa.nt.set2$AASeq, "PeptideID"])
length(peptide.data$PeptideID) - length(aa.nt.set2$PeptideID)
str_length(peptide.data$AASeq[1])
str_length(aa.nt.set1$AASeq[1])
str_length(aa.nt.set2$AASeq[2])

# Looks like set1 has AA seqs that are one AA too long. I think that's just a mistake, since everything
# else is the same. I'll stick with set2.
rm(aa.nt.set1)

# Pulling out the Peptide IDs that are not in nt.seqs.df.
id.missing <- peptide.data[!peptide.data$AASeq %in% nt.seqs.df$aas, "PeptideID"]
id.missing
id.missing %in% aa.nt.set2$PeptideID
# It's also not in the other set of peptide IDs.

# Excluding the one peptide with the missing NT sequence for the moment.
incomplete.peptide.data <- peptide.data[peptide.data$PeptideID != "PEPNR00000006745",]
length(incomplete.peptide.data$PeptideID)

# Calculating GC content.
for (i in 1:length(incomplete.peptide.data$PeptideID)) {
  temp.seq.set <- nt.seqs.df[as.character(nt.seqs.df$aas) == as.character(incomplete.peptide.data$AASeq[i]), "nts"]
  incomplete.peptide.data$GC.avg[i] <- mean.gc(DNAStringSet(temp.seq.set$nts))
}

# Selecting highest weight peptide per cluster. Which means it won't matter that the one peptide is missing,
# since it's not the max weight peptide for its cluster.
peptides.maxweight <- 
  incomplete.peptide.data %>%
  group_by(Cluster) %>%
  filter(Weight.nb == max(Weight.nb))

# Now checking average GC content.
mean(peptides.maxweight$GC.avg)

# Checking correlation with fitness.
weightedCorr(incomplete.peptide.data$Fitness.nb,
             incomplete.peptide.data$GC.avg,
             method = "Spearman", weights = incomplete.peptide.data$Weight.nb)

# Checking variance, min, and max.
stats::var(peptides.maxweight$GC.avg) # Seems low.
min(peptides.maxweight$GC.avg)
max(peptides.maxweight$GC.avg)
quantile(peptides.maxweight$GC.avg)

# Checking the full data set's GC content.
mean.gc(nt.fasta)

# What is the GC content of the unique peptide sequences?
unique.seqs.df <- nt.seqs.df %>% distinct(aas, .keep_all = T)
unique.seqs.df
mean.gc(DNAStringSet(unique.seqs.df$nts)) # Still about 59%

# Updating the supplemental table with GC content.
peptide.data$GC.avg <- NA
for (i in 1:length(peptide.data$PeptideID)) {
  if (peptide.data$PeptideID[i] %in% incomplete.peptide.data$PeptideID){
    peptide.data$GC.avg[i] <-
      incomplete.peptide.data[incomplete.peptide.data$PeptideID == peptide.data$PeptideID[i], "GC.avg"]
  }
}
peptide.data[is.na(peptide.data$GC.avg), c("PeptideID", "GC.avg")]

write_tsv(peptide.data, path = "Data/peptide_data_clusters_4-9-20.tsv")
