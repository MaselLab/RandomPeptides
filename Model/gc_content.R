# Figuring out the GC content of the 1000 random peptides.

# Load packages.
library(tidyverse)
library(Biostrings)
library(wCorr)

# Importing random peptide data.
peptide.data <- read.table(file = "Data/peptide_data_clusters_2-14-20.tsv", header = T, stringsAsFactors = F)

# Importing nucleotide sequences and matching them up with the appropriate peptide ID.
nt.fasta <- readDNAStringSet("Data/ORF_NR99_DB.fa")
translate(nt.fasta[1])[[1]]
nt.seqs.df <- tibble("id" = names(nt.fasta),
                     "nts" = paste(nt.fasta),
                     "aas" = paste(translate(nt.fasta, if.fuzzy.codon = "solve")))
# The following command returned multiple matches.
nt.seqs.df[as.character(nt.seqs.df$aas) == as.character(peptide.data$AASeq[1]), "nts"]
# This means I can't actually figure out the exact nucleotide sequence of the peptides I have. Uh oh.
length(unique(nt.seqs.df$aas))
length(nt.seqs.df$aas)
length(unique(nt.seqs.df$nts))
#merged.sequences.df <- merge(nt.seqs.df, peptide.data[, c("PeptideID", "AASeq")])

# Working out GC content of each nucleotide sequence.
alphabetFrequency(DNAString(nt.seqs.df$nts[1]))[c("G", "C", "A", "T")]

# GC content function. Takes an input DNAString and calculates the GC content.
# Does not account for unknown nucleotides!
# Only takes one sequence at time!
calculate.gc <- function(dna.seq){
  require(Biostrings)
  gc.totals <- alphabetFrequency(dna.seq)[c("G", "C")]
  gc.sum <- sum(gc.totals)
  nt.totals <- alphabetFrequency(dna.seq)[c("G", "C", "A", "T")]
  nt.sum <- sum(nt.totals)
  gc.content <- gc.sum / nt.sum
  return(gc.content)
}

calculate.gc(DNAString(nt.seqs.df$nts[1]))

# Calculate mean GC content for a set of sequences.
mean.gc <- function(seqs.nt) {
  require(Biostrings)
  seqs.df <- data.frame(
    "nts" = seqs.nt,
    "gc" = NA
  )
  for (i in 1:length(seqs.nt)) {
    seqs.df[i, "gc"] <- calculate.gc(DNAString(seqs.df[i, "nts"][1]))
  }
  gc.avg <- base::mean(seqs.df[, "gc"])
  #print(head(seqs.df))
  #print(tail(seqs.df))
  return(gc.avg)
}

# Test seq set.
test.seq.set <- nt.seqs.df[as.character(nt.seqs.df$aas) == as.character(peptide.data$AASeq[1]), "nts"]
test.seq.set
length(test.seq.set$nts)
test.seq.df <- data.frame("nts" = test.seq.set, "gc" = NA)
head(test.seq.df)
test.seq.df[1, "nt"]
calculate.gc(DNAString(test.seq.df[1, "nts"][1]))
test.seq.df[1, "gc"] <- calculate.gc(DNAString(test.seq.df[1, "nts"][1]))
test.seq.df[1,]
test.seq.df[is.na(test.seq.df$nts),]
mean.gc(test.seq.set$nts)

# Calculating mean GC for each peptide in our dataset (but not the full set of 2 million + sequences).
peptide.data$GC.avg <- NA
temp.seq.set <- NA
for (i in 1:length(peptide.data$PeptideID)) {
  temp.seq.set <- nt.seqs.df[as.character(nt.seqs.df$aas) == as.character(peptide.data$AASeq[i]), "nts"]
  peptide.data$GC.avg[i] <- mean.gc(temp.seq.set$nts)
}

# Calculating GC content for each nucleotide sequence.
# Nevermind, this takes too long.
# nt.seqs.df$gc <- NA
# for (i in 1:length(nt.seqs.df$id)) {
#   nt.seqs.df$gc[i] <- calculate.gc(DNAString(nt.seqs.df$nts[i]))
# }

# I need to figure this out.
# Importing more data. First, the full list of peptides.
full.aa.data <- read_tsv(file = "Data/PEP_NR_DB.tab", col_names = F)
names(full.aa.data) <- c("PeptideID", "AASeq")
full.aa.data

peptide.data[!peptide.data$AASeq %in% full.aa.data$AASeq, c("PeptideID", "AASeq", "Cluster")]
nt.seqs.df[!nt.seqs.df$aas %in% full.aa.data$AASeq,]
full.aa.data[!full.aa.data$AASeq %in% nt.seqs.df$aas,] # There are 44796 seqs missing from nt.seqs.df.
length(full.aa.data$PeptideID) - length(unique(nt.seqs.df$aas)) # Same number missing. OK.

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
id.missing %in% aa.nt.set2$PeptideID
# It's also not in the other set of peptide IDs.

# Excluding the one peptide with the missing NT sequence for the moment.
incomplete.peptide.data <- peptide.data[peptide.data$PeptideID != "PEPNR00000006745",]
length(incomplete.peptide.data$PeptideID)

# Calculating GC content.
for (i in 1:length(incomplete.peptide.data$PeptideID)) {
  temp.seq.set <- nt.seqs.df[as.character(nt.seqs.df$aas) == as.character(incomplete.peptide.data$AASeq[i]), "nts"]
  incomplete.peptide.data$GC.avg[i] <- mean.gc(temp.seq.set$nts)
}

# Selecting highest weight peptide per cluster.
peptides.maxweight <- 
  incomplete.peptide.data %>%
  group_by(Cluster) %>%
  filter(Weight.nb == max(Weight.nb))

# Now checking average GC content.
mean(peptides.maxweight$GC.avg)

weightedCorr(incomplete.peptide.data$Fitness.nb,
             incomplete.peptide.data$GC.avg,
             method = "Spearman", weights = incomplete.peptide.data$Weight.nb)

# Checking the full data set's GC content.
nt.fasta[1:10]
calculate.gc(nt.fasta[1])

mean.gc(nt.fasta[1:100])

# That takes too long. Remaking my function to be more efficient.
gc.mean.faster <- function(nts.set){
  require(Biostrings)
  nts.counts <- alphabetFrequency(nts.set)[,c("G", "C", "A", "T")]
  nts.total <- sum(nts.counts)
  gc.total <- sum(nts.counts[,c("G", "C")])
  return(gc.total / nts.total)
}

# Checking again.
start.time <- proc.time()
gc.mean.faster(nt.fasta[1:100000])
proc.time() - start.time

# Way better.
gc.mean.faster(nt.fasta) # 56% GC content for the whole dataset! That is strange...
