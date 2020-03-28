# Machine learning script.
# Luke Kosinski
# lkosinsk@gmail.com

# Load packages.
library(ggplot2)
library(caret)
library(tidyverse)
library(neuralnet)
library(nnet)
library(lme4)
#library(DECIPHER)

# Load peptide data.
peptide.data <- read.table(file = "Data/peptide_data_clusters_2-14-20.tsv", header = T, stringsAsFactors = F)
peptide.data <- as_tibble(peptide.data)
peptide.data
glimpse(peptide.data)
str(peptide.data)
summary(peptide.data)

# Taking the largest weight sequence per cluster.
peptide.reduced <- dplyr::select(peptide.data,
                              PeptideID,
                              Fitness.nb,
                              Weight.nb,
                              Cluster,
                              AASeq)
peptide.clusters.maxweight <- 
  peptide.reduced %>%
  group_by(Cluster) %>%
  filter(Weight.nb == max(Weight.nb))

peptide.clusters.maxweight

# Converting sequences into categorical data that a neural network can process.
# Joanna's suggestion: make 20 variables for each of the 50 random peptide positions,
# one for each possible amino acid in each position.
# More tycost-ecolical machine learning approach: use one-hot encoding to make 50 dummy variables.
# I'll try the one hot encoding first.
peptide.clusters.maxweight$Pos1 <- base::substr(peptide.clusters.maxweight$AASeq, 5, 5)
peptide.clusters.maxweight$Pos2 <- base::substr(peptide.clusters.maxweight$AASeq, 6, 6)
peptide.clusters.maxweight$Pos3 <- base::substr(peptide.clusters.maxweight$AASeq, 7, 7)
peptide.clusters.maxweight$Pos4 <- base::substr(peptide.clusters.maxweight$AASeq, 8, 8)
peptide.clusters.maxweight$Pos5 <- base::substr(peptide.clusters.maxweight$AASeq, 9, 9)
peptide.clusters.maxweight$Pos6 <- base::substr(peptide.clusters.maxweight$AASeq, 10, 10)
peptide.clusters.maxweight$Pos7 <- base::substr(peptide.clusters.maxweight$AASeq, 11, 11)
peptide.clusters.maxweight$Pos8 <- base::substr(peptide.clusters.maxweight$AASeq, 12, 12)
peptide.clusters.maxweight$Pos9 <- base::substr(peptide.clusters.maxweight$AASeq, 13, 13)
peptide.clusters.maxweight$Pos10 <- base::substr(peptide.clusters.maxweight$AASeq, 14, 14)
peptide.clusters.maxweight$Pos11 <- base::substr(peptide.clusters.maxweight$AASeq, 15, 15)
peptide.clusters.maxweight$Pos12 <- base::substr(peptide.clusters.maxweight$AASeq, 16, 16)
peptide.clusters.maxweight$Pos13 <- base::substr(peptide.clusters.maxweight$AASeq, 17, 17)
peptide.clusters.maxweight$Pos14 <- base::substr(peptide.clusters.maxweight$AASeq, 18, 18)
peptide.clusters.maxweight$Pos15 <- base::substr(peptide.clusters.maxweight$AASeq, 19, 19)
peptide.clusters.maxweight$Pos16 <- base::substr(peptide.clusters.maxweight$AASeq, 20, 20)
peptide.clusters.maxweight$Pos17 <- base::substr(peptide.clusters.maxweight$AASeq, 21, 21)
peptide.clusters.maxweight$Pos18 <- base::substr(peptide.clusters.maxweight$AASeq, 22, 22)
peptide.clusters.maxweight$Pos19 <- base::substr(peptide.clusters.maxweight$AASeq, 23, 23)
peptide.clusters.maxweight$Pos20 <- base::substr(peptide.clusters.maxweight$AASeq, 24, 24)
peptide.clusters.maxweight$Pos21 <- base::substr(peptide.clusters.maxweight$AASeq, 25, 25)
peptide.clusters.maxweight$Pos22 <- base::substr(peptide.clusters.maxweight$AASeq, 26, 26)
peptide.clusters.maxweight$Pos23 <- base::substr(peptide.clusters.maxweight$AASeq, 27, 27)
peptide.clusters.maxweight$Pos24 <- base::substr(peptide.clusters.maxweight$AASeq, 28, 28)
peptide.clusters.maxweight$Pos25 <- base::substr(peptide.clusters.maxweight$AASeq, 29, 29)
peptide.clusters.maxweight$Pos26 <- base::substr(peptide.clusters.maxweight$AASeq, 30, 30)
peptide.clusters.maxweight$Pos27 <- base::substr(peptide.clusters.maxweight$AASeq, 31, 31)
peptide.clusters.maxweight$Pos28 <- base::substr(peptide.clusters.maxweight$AASeq, 32, 32)
peptide.clusters.maxweight$Pos29 <- base::substr(peptide.clusters.maxweight$AASeq, 33, 33)
peptide.clusters.maxweight$Pos30 <- base::substr(peptide.clusters.maxweight$AASeq, 34, 34)
peptide.clusters.maxweight$Pos31 <- base::substr(peptide.clusters.maxweight$AASeq, 35, 35)
peptide.clusters.maxweight$Pos32 <- base::substr(peptide.clusters.maxweight$AASeq, 36, 36)
peptide.clusters.maxweight$Pos33 <- base::substr(peptide.clusters.maxweight$AASeq, 37, 37)
peptide.clusters.maxweight$Pos34 <- base::substr(peptide.clusters.maxweight$AASeq, 38, 38)
peptide.clusters.maxweight$Pos35 <- base::substr(peptide.clusters.maxweight$AASeq, 39, 39)
peptide.clusters.maxweight$Pos36 <- base::substr(peptide.clusters.maxweight$AASeq, 40, 40)
peptide.clusters.maxweight$Pos37 <- base::substr(peptide.clusters.maxweight$AASeq, 41, 41)
peptide.clusters.maxweight$Pos38 <- base::substr(peptide.clusters.maxweight$AASeq, 42, 42)
peptide.clusters.maxweight$Pos39 <- base::substr(peptide.clusters.maxweight$AASeq, 43, 43)
peptide.clusters.maxweight$Pos40 <- base::substr(peptide.clusters.maxweight$AASeq, 44, 44)
peptide.clusters.maxweight$Pos41 <- base::substr(peptide.clusters.maxweight$AASeq, 45, 45)
peptide.clusters.maxweight$Pos42 <- base::substr(peptide.clusters.maxweight$AASeq, 46, 46)
peptide.clusters.maxweight$Pos43 <- base::substr(peptide.clusters.maxweight$AASeq, 47, 47)
peptide.clusters.maxweight$Pos44 <- base::substr(peptide.clusters.maxweight$AASeq, 48, 48)
peptide.clusters.maxweight$Pos45 <- base::substr(peptide.clusters.maxweight$AASeq, 49, 49)
peptide.clusters.maxweight$Pos46 <- base::substr(peptide.clusters.maxweight$AASeq, 50, 50)
peptide.clusters.maxweight$Pos47 <- base::substr(peptide.clusters.maxweight$AASeq, 51, 51)
peptide.clusters.maxweight$Pos48 <- base::substr(peptide.clusters.maxweight$AASeq, 52, 52)
peptide.clusters.maxweight$Pos49 <- base::substr(peptide.clusters.maxweight$AASeq, 53, 53)
peptide.clusters.maxweight$Pos50 <- base::substr(peptide.clusters.maxweight$AASeq, 54, 54)
peptide.clusters.maxweight

peptide.catvars <- peptide.clusters.maxweight[,c(2:3, 6:55)]

seq.dmy <- dummyVars(" ~ .", data = peptide.catvars)
seq.dmy.df <- data.frame(predict(seq.dmy, newdata = peptide.catvars))
head(seq.dmy.df, 1)
seq.dmy.df$PeptideID <- peptide.clusters.maxweight$PeptideID

# Checking for linear dependencies.
seq.dmy.catonly <- seq.dmy.df[ ,c(-1, -2, -1003)]
seq.dmy.combos <- findLinearCombos(seq.dmy.catonly)
seq.dmy.combos$remove
seq.dmy.catonly.trimmed <- seq.dmy.catonly[, -seq.dmy.combos$remove]
seq.dmy.catonly.trimmed$PeptideID <- seq.dmy.df$PeptideID
seq.dmy.catonly.trimmed$Fitness.nb <- seq.dmy.df$Fitness.nb
seq.dmy.catonly.trimmed$Weight.nb <- seq.dmy.df$Weight.nb

# I can't use character data with the NeuralNet package. So, I'll go with AA counts.
# Also checking with the one-hot encoding for the dummy variables, which appears to be the
# same as Joanna's suggestion.
# peptide.counts <- select(peptide.data,
#                          PeptideID,
#                          Fitness.nb,
#                          Clustering.Six,
#                          WaltzBinary,
#                          Cluster,
#                          Weight.nb,
#                          Leu, Pro, Met, Trp, Ala,
#                          Val, Phe, Ile, Gly, Ser,
#                          Thr, Cys, Asn, Gln, Tyr,
#                          His, Asp, Glu, Lys, Arg)
# sum(is.na(peptide.counts))

# Making one more data frame of properties for each AA in each position.
# peptide.aa.props <- dplyr::select(peptide.data,
#                                   PeptideID,
#                                   Fitness.nb,
#                                   Weight.nb,
#                                   Cluster,
#                                   AASeq)
# peptide.aa.props$Pos1.RSA <- mean.metric.calculator(peptide.seqs$Pos1, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos2.RSA <- mean.metric.calculator(peptide.seqs$Pos2, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos3.RSA <- mean.metric.calculator(peptide.seqs$Pos3, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos4.RSA <- mean.metric.calculator(peptide.seqs$Pos4, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos5.RSA <- mean.metric.calculator(peptide.seqs$Pos5, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos6.RSA <- mean.metric.calculator(peptide.seqs$Pos6, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos7.RSA <- mean.metric.calculator(peptide.seqs$Pos7, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos8.RSA <- mean.metric.calculator(peptide.seqs$Pos8, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos9.RSA <- mean.metric.calculator(peptide.seqs$Pos9, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos10.RSA <- mean.metric.calculator(peptide.seqs$Pos10, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos11.RSA <- mean.metric.calculator(peptide.seqs$Pos11, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos12.RSA <- mean.metric.calculator(peptide.seqs$Pos12, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos13.RSA <- mean.metric.calculator(peptide.seqs$Pos13, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos14.RSA <- mean.metric.calculator(peptide.seqs$Pos14, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos15.RSA <- mean.metric.calculator(peptide.seqs$Pos15, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos16.RSA <- mean.metric.calculator(peptide.seqs$Pos16, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos17.RSA <- mean.metric.calculator(peptide.seqs$Pos17, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos18.RSA <- mean.metric.calculator(peptide.seqs$Pos18, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos19.RSA <- mean.metric.calculator(peptide.seqs$Pos19, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos20.RSA <- mean.metric.calculator(peptide.seqs$Pos20, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos21.RSA <- mean.metric.calculator(peptide.seqs$Pos21, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos22.RSA <- mean.metric.calculator(peptide.seqs$Pos22, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos23.RSA <- mean.metric.calculator(peptide.seqs$Pos23, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos24.RSA <- mean.metric.calculator(peptide.seqs$Pos24, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos25.RSA <- mean.metric.calculator(peptide.seqs$Pos25, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos26.RSA <- mean.metric.calculator(peptide.seqs$Pos26, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos27.RSA <- mean.metric.calculator(peptide.seqs$Pos27, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos28.RSA <- mean.metric.calculator(peptide.seqs$Pos28, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos29.RSA <- mean.metric.calculator(peptide.seqs$Pos29, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos30.RSA <- mean.metric.calculator(peptide.seqs$Pos30, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos31.RSA <- mean.metric.calculator(peptide.seqs$Pos31, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos32.RSA <- mean.metric.calculator(peptide.seqs$Pos32, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos33.RSA <- mean.metric.calculator(peptide.seqs$Pos33, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos34.RSA <- mean.metric.calculator(peptide.seqs$Pos34, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos35.RSA <- mean.metric.calculator(peptide.seqs$Pos35, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos36.RSA <- mean.metric.calculator(peptide.seqs$Pos36, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos37.RSA <- mean.metric.calculator(peptide.seqs$Pos37, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos38.RSA <- mean.metric.calculator(peptide.seqs$Pos38, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos39.RSA <- mean.metric.calculator(peptide.seqs$Pos39, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos40.RSA <- mean.metric.calculator(peptide.seqs$Pos40, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos41.RSA <- mean.metric.calculator(peptide.seqs$Pos41, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos42.RSA <- mean.metric.calculator(peptide.seqs$Pos42, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos43.RSA <- mean.metric.calculator(peptide.seqs$Pos43, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos44.RSA <- mean.metric.calculator(peptide.seqs$Pos44, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos45.RSA <- mean.metric.calculator(peptide.seqs$Pos45, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos46.RSA <- mean.metric.calculator(peptide.seqs$Pos46, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos47.RSA <- mean.metric.calculator(peptide.seqs$Pos47, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos48.RSA <- mean.metric.calculator(peptide.seqs$Pos48, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos49.RSA <- mean.metric.calculator(peptide.seqs$Pos49, metric = "RSA-hydrophilicity")
# peptide.aa.props$Pos50.RSA <- mean.metric.calculator(peptide.seqs$Pos50, metric = "RSA-hydrophilicity")
# 
# peptide.aa.props$Pos1.area <- mean.metric.calculator(peptide.seqs$Pos1, metric = "area")
# peptide.aa.props$Pos2.area <- mean.metric.calculator(peptide.seqs$Pos2, metric = "area")
# peptide.aa.props$Pos3.area <- mean.metric.calculator(peptide.seqs$Pos3, metric = "area")
# peptide.aa.props$Pos4.area <- mean.metric.calculator(peptide.seqs$Pos4, metric = "area")
# peptide.aa.props$Pos5.area <- mean.metric.calculator(peptide.seqs$Pos5, metric = "area")
# peptide.aa.props$Pos6.area <- mean.metric.calculator(peptide.seqs$Pos6, metric = "area")
# peptide.aa.props$Pos7.area <- mean.metric.calculator(peptide.seqs$Pos7, metric = "area")
# peptide.aa.props$Pos8.area <- mean.metric.calculator(peptide.seqs$Pos8, metric = "area")
# peptide.aa.props$Pos9.area <- mean.metric.calculator(peptide.seqs$Pos9, metric = "area")
# peptide.aa.props$Pos10.area <- mean.metric.calculator(peptide.seqs$Pos10, metric = "area")
# peptide.aa.props$Pos11.area <- mean.metric.calculator(peptide.seqs$Pos11, metric = "area")
# peptide.aa.props$Pos12.area <- mean.metric.calculator(peptide.seqs$Pos12, metric = "area")
# peptide.aa.props$Pos13.area <- mean.metric.calculator(peptide.seqs$Pos13, metric = "area")
# peptide.aa.props$Pos14.area <- mean.metric.calculator(peptide.seqs$Pos14, metric = "area")
# peptide.aa.props$Pos15.area <- mean.metric.calculator(peptide.seqs$Pos15, metric = "area")
# peptide.aa.props$Pos16.area <- mean.metric.calculator(peptide.seqs$Pos16, metric = "area")
# peptide.aa.props$Pos17.area <- mean.metric.calculator(peptide.seqs$Pos17, metric = "area")
# peptide.aa.props$Pos18.area <- mean.metric.calculator(peptide.seqs$Pos18, metric = "area")
# peptide.aa.props$Pos19.area <- mean.metric.calculator(peptide.seqs$Pos19, metric = "area")
# peptide.aa.props$Pos20.area <- mean.metric.calculator(peptide.seqs$Pos20, metric = "area")
# peptide.aa.props$Pos21.area <- mean.metric.calculator(peptide.seqs$Pos21, metric = "area")
# peptide.aa.props$Pos22.area <- mean.metric.calculator(peptide.seqs$Pos22, metric = "area")
# peptide.aa.props$Pos23.area <- mean.metric.calculator(peptide.seqs$Pos23, metric = "area")
# peptide.aa.props$Pos24.area <- mean.metric.calculator(peptide.seqs$Pos24, metric = "area")
# peptide.aa.props$Pos25.area <- mean.metric.calculator(peptide.seqs$Pos25, metric = "area")
# peptide.aa.props$Pos26.area <- mean.metric.calculator(peptide.seqs$Pos26, metric = "area")
# peptide.aa.props$Pos27.area <- mean.metric.calculator(peptide.seqs$Pos27, metric = "area")
# peptide.aa.props$Pos28.area <- mean.metric.calculator(peptide.seqs$Pos28, metric = "area")
# peptide.aa.props$Pos29.area <- mean.metric.calculator(peptide.seqs$Pos29, metric = "area")
# peptide.aa.props$Pos30.area <- mean.metric.calculator(peptide.seqs$Pos30, metric = "area")
# peptide.aa.props$Pos31.area <- mean.metric.calculator(peptide.seqs$Pos31, metric = "area")
# peptide.aa.props$Pos32.area <- mean.metric.calculator(peptide.seqs$Pos32, metric = "area")
# peptide.aa.props$Pos33.area <- mean.metric.calculator(peptide.seqs$Pos33, metric = "area")
# peptide.aa.props$Pos34.area <- mean.metric.calculator(peptide.seqs$Pos34, metric = "area")
# peptide.aa.props$Pos35.area <- mean.metric.calculator(peptide.seqs$Pos35, metric = "area")
# peptide.aa.props$Pos36.area <- mean.metric.calculator(peptide.seqs$Pos36, metric = "area")
# peptide.aa.props$Pos37.area <- mean.metric.calculator(peptide.seqs$Pos37, metric = "area")
# peptide.aa.props$Pos38.area <- mean.metric.calculator(peptide.seqs$Pos38, metric = "area")
# peptide.aa.props$Pos39.area <- mean.metric.calculator(peptide.seqs$Pos39, metric = "area")
# peptide.aa.props$Pos40.area <- mean.metric.calculator(peptide.seqs$Pos40, metric = "area")
# peptide.aa.props$Pos41.area <- mean.metric.calculator(peptide.seqs$Pos41, metric = "area")
# peptide.aa.props$Pos42.area <- mean.metric.calculator(peptide.seqs$Pos42, metric = "area")
# peptide.aa.props$Pos43.area <- mean.metric.calculator(peptide.seqs$Pos43, metric = "area")
# peptide.aa.props$Pos44.area <- mean.metric.calculator(peptide.seqs$Pos44, metric = "area")
# peptide.aa.props$Pos45.area <- mean.metric.calculator(peptide.seqs$Pos45, metric = "area")
# peptide.aa.props$Pos46.area <- mean.metric.calculator(peptide.seqs$Pos46, metric = "area")
# peptide.aa.props$Pos47.area <- mean.metric.calculator(peptide.seqs$Pos47, metric = "area")
# peptide.aa.props$Pos48.area <- mean.metric.calculator(peptide.seqs$Pos48, metric = "area")
# peptide.aa.props$Pos49.area <- mean.metric.calculator(peptide.seqs$Pos49, metric = "area")
# peptide.aa.props$Pos50.area <- mean.metric.calculator(peptide.seqs$Pos50, metric = "area")
# 
# peptide.aa.props$Pos1.weight <- mean.metric.calculator(peptide.seqs$Pos1, metric = "weight")
# peptide.aa.props$Pos2.weight <- mean.metric.calculator(peptide.seqs$Pos2, metric = "weight")
# peptide.aa.props$Pos3.weight <- mean.metric.calculator(peptide.seqs$Pos3, metric = "weight")
# peptide.aa.props$Pos4.weight <- mean.metric.calculator(peptide.seqs$Pos4, metric = "weight")
# peptide.aa.props$Pos5.weight <- mean.metric.calculator(peptide.seqs$Pos5, metric = "weight")
# peptide.aa.props$Pos6.weight <- mean.metric.calculator(peptide.seqs$Pos6, metric = "weight")
# peptide.aa.props$Pos7.weight <- mean.metric.calculator(peptide.seqs$Pos7, metric = "weight")
# peptide.aa.props$Pos8.weight <- mean.metric.calculator(peptide.seqs$Pos8, metric = "weight")
# peptide.aa.props$Pos9.weight <- mean.metric.calculator(peptide.seqs$Pos9, metric = "weight")
# peptide.aa.props$Pos10.weight <- mean.metric.calculator(peptide.seqs$Pos10, metric = "weight")
# peptide.aa.props$Pos11.weight <- mean.metric.calculator(peptide.seqs$Pos11, metric = "weight")
# peptide.aa.props$Pos12.weight <- mean.metric.calculator(peptide.seqs$Pos12, metric = "weight")
# peptide.aa.props$Pos13.weight <- mean.metric.calculator(peptide.seqs$Pos13, metric = "weight")
# peptide.aa.props$Pos14.weight <- mean.metric.calculator(peptide.seqs$Pos14, metric = "weight")
# peptide.aa.props$Pos15.weight <- mean.metric.calculator(peptide.seqs$Pos15, metric = "weight")
# peptide.aa.props$Pos16.weight <- mean.metric.calculator(peptide.seqs$Pos16, metric = "weight")
# peptide.aa.props$Pos17.weight <- mean.metric.calculator(peptide.seqs$Pos17, metric = "weight")
# peptide.aa.props$Pos18.weight <- mean.metric.calculator(peptide.seqs$Pos18, metric = "weight")
# peptide.aa.props$Pos19.weight <- mean.metric.calculator(peptide.seqs$Pos19, metric = "weight")
# peptide.aa.props$Pos20.weight <- mean.metric.calculator(peptide.seqs$Pos20, metric = "weight")
# peptide.aa.props$Pos21.weight <- mean.metric.calculator(peptide.seqs$Pos21, metric = "weight")
# peptide.aa.props$Pos22.weight <- mean.metric.calculator(peptide.seqs$Pos22, metric = "weight")
# peptide.aa.props$Pos23.weight <- mean.metric.calculator(peptide.seqs$Pos23, metric = "weight")
# peptide.aa.props$Pos24.weight <- mean.metric.calculator(peptide.seqs$Pos24, metric = "weight")
# peptide.aa.props$Pos25.weight <- mean.metric.calculator(peptide.seqs$Pos25, metric = "weight")
# peptide.aa.props$Pos26.weight <- mean.metric.calculator(peptide.seqs$Pos26, metric = "weight")
# peptide.aa.props$Pos27.weight <- mean.metric.calculator(peptide.seqs$Pos27, metric = "weight")
# peptide.aa.props$Pos28.weight <- mean.metric.calculator(peptide.seqs$Pos28, metric = "weight")
# peptide.aa.props$Pos29.weight <- mean.metric.calculator(peptide.seqs$Pos29, metric = "weight")
# peptide.aa.props$Pos30.weight <- mean.metric.calculator(peptide.seqs$Pos30, metric = "weight")
# peptide.aa.props$Pos31.weight <- mean.metric.calculator(peptide.seqs$Pos31, metric = "weight")
# peptide.aa.props$Pos32.weight <- mean.metric.calculator(peptide.seqs$Pos32, metric = "weight")
# peptide.aa.props$Pos33.weight <- mean.metric.calculator(peptide.seqs$Pos33, metric = "weight")
# peptide.aa.props$Pos34.weight <- mean.metric.calculator(peptide.seqs$Pos34, metric = "weight")
# peptide.aa.props$Pos35.weight <- mean.metric.calculator(peptide.seqs$Pos35, metric = "weight")
# peptide.aa.props$Pos36.weight <- mean.metric.calculator(peptide.seqs$Pos36, metric = "weight")
# peptide.aa.props$Pos37.weight <- mean.metric.calculator(peptide.seqs$Pos37, metric = "weight")
# peptide.aa.props$Pos38.weight <- mean.metric.calculator(peptide.seqs$Pos38, metric = "weight")
# peptide.aa.props$Pos39.weight <- mean.metric.calculator(peptide.seqs$Pos39, metric = "weight")
# peptide.aa.props$Pos40.weight <- mean.metric.calculator(peptide.seqs$Pos40, metric = "weight")
# peptide.aa.props$Pos41.weight <- mean.metric.calculator(peptide.seqs$Pos41, metric = "weight")
# peptide.aa.props$Pos42.weight <- mean.metric.calculator(peptide.seqs$Pos42, metric = "weight")
# peptide.aa.props$Pos43.weight <- mean.metric.calculator(peptide.seqs$Pos43, metric = "weight")
# peptide.aa.props$Pos44.weight <- mean.metric.calculator(peptide.seqs$Pos44, metric = "weight")
# peptide.aa.props$Pos45.weight <- mean.metric.calculator(peptide.seqs$Pos45, metric = "weight")
# peptide.aa.props$Pos46.weight <- mean.metric.calculator(peptide.seqs$Pos46, metric = "weight")
# peptide.aa.props$Pos47.weight <- mean.metric.calculator(peptide.seqs$Pos47, metric = "weight")
# peptide.aa.props$Pos48.weight <- mean.metric.calculator(peptide.seqs$Pos48, metric = "weight")
# peptide.aa.props$Pos49.weight <- mean.metric.calculator(peptide.seqs$Pos49, metric = "weight")
# peptide.aa.props$Pos50.weight <- mean.metric.calculator(peptide.seqs$Pos50, metric = "weight")
# 
# peptide.aa.props$Pos1.pI <- mean.metric.calculator(peptide.seqs$Pos1, metric = "pI")
# peptide.aa.props$Pos2.pI <- mean.metric.calculator(peptide.seqs$Pos2, metric = "pI")
# peptide.aa.props$Pos3.pI <- mean.metric.calculator(peptide.seqs$Pos3, metric = "pI")
# peptide.aa.props$Pos4.pI <- mean.metric.calculator(peptide.seqs$Pos4, metric = "pI")
# peptide.aa.props$Pos5.pI <- mean.metric.calculator(peptide.seqs$Pos5, metric = "pI")
# peptide.aa.props$Pos6.pI <- mean.metric.calculator(peptide.seqs$Pos6, metric = "pI")
# peptide.aa.props$Pos7.pI <- mean.metric.calculator(peptide.seqs$Pos7, metric = "pI")
# peptide.aa.props$Pos8.pI <- mean.metric.calculator(peptide.seqs$Pos8, metric = "pI")
# peptide.aa.props$Pos9.pI <- mean.metric.calculator(peptide.seqs$Pos9, metric = "pI")
# peptide.aa.props$Pos10.pI <- mean.metric.calculator(peptide.seqs$Pos10, metric = "pI")
# peptide.aa.props$Pos11.pI <- mean.metric.calculator(peptide.seqs$Pos11, metric = "pI")
# peptide.aa.props$Pos12.pI <- mean.metric.calculator(peptide.seqs$Pos12, metric = "pI")
# peptide.aa.props$Pos13.pI <- mean.metric.calculator(peptide.seqs$Pos13, metric = "pI")
# peptide.aa.props$Pos14.pI <- mean.metric.calculator(peptide.seqs$Pos14, metric = "pI")
# peptide.aa.props$Pos15.pI <- mean.metric.calculator(peptide.seqs$Pos15, metric = "pI")
# peptide.aa.props$Pos16.pI <- mean.metric.calculator(peptide.seqs$Pos16, metric = "pI")
# peptide.aa.props$Pos17.pI <- mean.metric.calculator(peptide.seqs$Pos17, metric = "pI")
# peptide.aa.props$Pos18.pI <- mean.metric.calculator(peptide.seqs$Pos18, metric = "pI")
# peptide.aa.props$Pos19.pI <- mean.metric.calculator(peptide.seqs$Pos19, metric = "pI")
# peptide.aa.props$Pos20.pI <- mean.metric.calculator(peptide.seqs$Pos20, metric = "pI")
# peptide.aa.props$Pos21.pI <- mean.metric.calculator(peptide.seqs$Pos21, metric = "pI")
# peptide.aa.props$Pos22.pI <- mean.metric.calculator(peptide.seqs$Pos22, metric = "pI")
# peptide.aa.props$Pos23.pI <- mean.metric.calculator(peptide.seqs$Pos23, metric = "pI")
# peptide.aa.props$Pos24.pI <- mean.metric.calculator(peptide.seqs$Pos24, metric = "pI")
# peptide.aa.props$Pos25.pI <- mean.metric.calculator(peptide.seqs$Pos25, metric = "pI")
# peptide.aa.props$Pos26.pI <- mean.metric.calculator(peptide.seqs$Pos26, metric = "pI")
# peptide.aa.props$Pos27.pI <- mean.metric.calculator(peptide.seqs$Pos27, metric = "pI")
# peptide.aa.props$Pos28.pI <- mean.metric.calculator(peptide.seqs$Pos28, metric = "pI")
# peptide.aa.props$Pos29.pI <- mean.metric.calculator(peptide.seqs$Pos29, metric = "pI")
# peptide.aa.props$Pos30.pI <- mean.metric.calculator(peptide.seqs$Pos30, metric = "pI")
# peptide.aa.props$Pos31.pI <- mean.metric.calculator(peptide.seqs$Pos31, metric = "pI")
# peptide.aa.props$Pos32.pI <- mean.metric.calculator(peptide.seqs$Pos32, metric = "pI")
# peptide.aa.props$Pos33.pI <- mean.metric.calculator(peptide.seqs$Pos33, metric = "pI")
# peptide.aa.props$Pos34.pI <- mean.metric.calculator(peptide.seqs$Pos34, metric = "pI")
# peptide.aa.props$Pos35.pI <- mean.metric.calculator(peptide.seqs$Pos35, metric = "pI")
# peptide.aa.props$Pos36.pI <- mean.metric.calculator(peptide.seqs$Pos36, metric = "pI")
# peptide.aa.props$Pos37.pI <- mean.metric.calculator(peptide.seqs$Pos37, metric = "pI")
# peptide.aa.props$Pos38.pI <- mean.metric.calculator(peptide.seqs$Pos38, metric = "pI")
# peptide.aa.props$Pos39.pI <- mean.metric.calculator(peptide.seqs$Pos39, metric = "pI")
# peptide.aa.props$Pos40.pI <- mean.metric.calculator(peptide.seqs$Pos40, metric = "pI")
# peptide.aa.props$Pos41.pI <- mean.metric.calculator(peptide.seqs$Pos41, metric = "pI")
# peptide.aa.props$Pos42.pI <- mean.metric.calculator(peptide.seqs$Pos42, metric = "pI")
# peptide.aa.props$Pos43.pI <- mean.metric.calculator(peptide.seqs$Pos43, metric = "pI")
# peptide.aa.props$Pos44.pI <- mean.metric.calculator(peptide.seqs$Pos44, metric = "pI")
# peptide.aa.props$Pos45.pI <- mean.metric.calculator(peptide.seqs$Pos45, metric = "pI")
# peptide.aa.props$Pos46.pI <- mean.metric.calculator(peptide.seqs$Pos46, metric = "pI")
# peptide.aa.props$Pos47.pI <- mean.metric.calculator(peptide.seqs$Pos47, metric = "pI")
# peptide.aa.props$Pos48.pI <- mean.metric.calculator(peptide.seqs$Pos48, metric = "pI")
# peptide.aa.props$Pos49.pI <- mean.metric.calculator(peptide.seqs$Pos49, metric = "pI")
# peptide.aa.props$Pos50.pI <- mean.metric.calculator(peptide.seqs$Pos50, metric = "pI")
# 
# peptide.aa.props$Pos1.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos1, metric = "cost-ecoli")
# peptide.aa.props$Pos2.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos2, metric = "cost-ecoli")
# peptide.aa.props$Pos3.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos3, metric = "cost-ecoli")
# peptide.aa.props$Pos4.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos4, metric = "cost-ecoli")
# peptide.aa.props$Pos5.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos5, metric = "cost-ecoli")
# peptide.aa.props$Pos6.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos6, metric = "cost-ecoli")
# peptide.aa.props$Pos7.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos7, metric = "cost-ecoli")
# peptide.aa.props$Pos8.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos8, metric = "cost-ecoli")
# peptide.aa.props$Pos9.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos9, metric = "cost-ecoli")
# peptide.aa.props$Pos10.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos10, metric = "cost-ecoli")
# peptide.aa.props$Pos11.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos11, metric = "cost-ecoli")
# peptide.aa.props$Pos12.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos12, metric = "cost-ecoli")
# peptide.aa.props$Pos13.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos13, metric = "cost-ecoli")
# peptide.aa.props$Pos14.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos14, metric = "cost-ecoli")
# peptide.aa.props$Pos15.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos15, metric = "cost-ecoli")
# peptide.aa.props$Pos16.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos16, metric = "cost-ecoli")
# peptide.aa.props$Pos17.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos17, metric = "cost-ecoli")
# peptide.aa.props$Pos18.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos18, metric = "cost-ecoli")
# peptide.aa.props$Pos19.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos19, metric = "cost-ecoli")
# peptide.aa.props$Pos20.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos20, metric = "cost-ecoli")
# peptide.aa.props$Pos21.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos21, metric = "cost-ecoli")
# peptide.aa.props$Pos22.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos22, metric = "cost-ecoli")
# peptide.aa.props$Pos23.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos23, metric = "cost-ecoli")
# peptide.aa.props$Pos24.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos24, metric = "cost-ecoli")
# peptide.aa.props$Pos25.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos25, metric = "cost-ecoli")
# peptide.aa.props$Pos26.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos26, metric = "cost-ecoli")
# peptide.aa.props$Pos27.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos27, metric = "cost-ecoli")
# peptide.aa.props$Pos28.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos28, metric = "cost-ecoli")
# peptide.aa.props$Pos29.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos29, metric = "cost-ecoli")
# peptide.aa.props$Pos30.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos30, metric = "cost-ecoli")
# peptide.aa.props$Pos31.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos31, metric = "cost-ecoli")
# peptide.aa.props$Pos32.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos32, metric = "cost-ecoli")
# peptide.aa.props$Pos33.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos33, metric = "cost-ecoli")
# peptide.aa.props$Pos34.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos34, metric = "cost-ecoli")
# peptide.aa.props$Pos35.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos35, metric = "cost-ecoli")
# peptide.aa.props$Pos36.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos36, metric = "cost-ecoli")
# peptide.aa.props$Pos37.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos37, metric = "cost-ecoli")
# peptide.aa.props$Pos38.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos38, metric = "cost-ecoli")
# peptide.aa.props$Pos39.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos39, metric = "cost-ecoli")
# peptide.aa.props$Pos40.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos40, metric = "cost-ecoli")
# peptide.aa.props$Pos41.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos41, metric = "cost-ecoli")
# peptide.aa.props$Pos42.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos42, metric = "cost-ecoli")
# peptide.aa.props$Pos43.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos43, metric = "cost-ecoli")
# peptide.aa.props$Pos44.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos44, metric = "cost-ecoli")
# peptide.aa.props$Pos45.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos45, metric = "cost-ecoli")
# peptide.aa.props$Pos46.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos46, metric = "cost-ecoli")
# peptide.aa.props$Pos47.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos47, metric = "cost-ecoli")
# peptide.aa.props$Pos48.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos48, metric = "cost-ecoli")
# peptide.aa.props$Pos49.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos49, metric = "cost-ecoli")
# peptide.aa.props$Pos50.cost.ecoli <- mean.metric.calculator(peptide.seqs$Pos50, metric = "cost-ecoli")
# 
# aa.properties.df <- peptide.aa.props[, c(2, 6:255)]

# Pre-processing. Normalize all variable before they go into the model. Center them on some median
# and then scale them to account for their distribution. We also want to split the data into
# training and test data sets.
set.seed(42)

# Partition for splitting data into training and test sets.
trainindex <- createDataPartition(
  seq.dmy.catonly.trimmed$Fitness.nb, p = 0.8, list = FALSE, times = 1
)

# Split data into training and test sets.
# peptide.train <- peptide.counts[trainindex, ]
# peptide.test <- peptide.counts[-trainindex, ]
seq.dmy.train <- seq.dmy.catonly.trimmed[trainindex,]
seq.dmy.test <- seq.dmy.catonly.trimmed[-trainindex,]
# aa.props.train <- aa.properties.df[trainindex,]
# aa.props.test <- aa.properties.df[-trainindex,]

# Pre-process our data, but not the column we are trying to predict. Standardizes all variable
# except Fitness.nb.
# pp <- preProcess(peptide.train[, -1],
#                  method = c("center", "scale"), # Can also add impute here by adding "knnImpute" or "bagImpute".
#                  outcome = peptide.train$Fitness.nb)
pp.dmy <- preProcess(seq.dmy.train,
                     method = c("center", "scale"),
                     outcome = seq.dmy.train$Fitness.nb)
# pp.aa.props <- preProcess(aa.props.train,
#                           method = c("center", "scale"),
#                           outcome = aa.props.train$Fitness.nb)
# peptide.train.pp <- predict(pp, peptide.train)
# peptide.test.pp <- predict(pp, peptide.test)
# glimpse(peptide.train.pp)
seq.dmy.train.pp <- predict(pp.dmy, seq.dmy.train)
seq.dmy.test.pp <- predict(pp.dmy, seq.dmy.test)
glimpse(seq.dmy.train.pp)
# aa.props.train.pp <- predict(pp.aa.props, aa.props.train)
# aa.props.test.pp <- predict(pp.aa.props, aa.props.test)

# Train using caret.
# peptide.ffnn <- train(
#   Fitness.nb ~ .,
#   data = seq.dmy.train.pp[,-c(647,649)],
#   method = "neuralnet"#,
#   #weights = seq.dmy.train.pp$Weight.nb,
#   #MaxNWts = 5000
# )
# peptide.ffnn
# fit.ffnn.pred <- predict(peptide.ffnn, seq.dmy.test.pp)
# fit.ffnn.pred
# cor.test(fit.ffnn.pred, seq.dmy.test.pp$Fitness.nb, method = "spearman")

# Train using neuralnet.
peptide.nn <- neuralnet(
  Fitness.nb ~ .,
  data = seq.dmy.train.pp[,-c(647,649)],
  startweights = seq.dmy.train$Weight.nb,
  hidden = 1
)
peptide.nn
fit.nn.pred <- predict(peptide.nn, seq.dmy.test.pp)
fit.nn.pred
cor.test(fit.nn.pred, seq.dmy.test.pp$Fitness.nb, method = "pearson")
summary(lm(formula = seq.dmy.test.pp$Fitness.nb ~ fit.nn.pred, weights = seq.dmy.test$Weight.nb))

# Train using nnet.
# peptide.ffnn.weights <- nnet(
#   Fitness.nb ~ .,
#   data = seq.dmy.train.pp[,-c(647,649)],
#   #weights = seq.dmy.train$Weight.nb,
#   size = 1
# )
# peptide.ffnn.weights
# predict(peptide.ffnn.weights, seq.dmy.test.pp)

# Nueral network test.
# peptide.ffnn <- neuralnet(
#   Fitness.nb ~ .,
#   data = seq.dmy.train.pp[,c(-2,-3)],
#   hidden = 1
# )
# peptide.ffnn
# plot(peptide.ffnn)

# Make prediction.
# peptide.pred <- compute(peptide.ffnn.weights, seq.dmy.test.pp)
# #fitness.pred <- peptide.pred$net.result
# fitness.pred <- predict(peptide.ffnn.weights, seq.dmy.test.pp)
# cor.test(seq.dmy.test.pp$Fitness.nb, fitness.pred, method = "pearson")
# peptide.test.pp$Pred.nn <- fitness.pred
# peptide.test$Pred.nn <- fitness.pred
# ggplot(
#   data = peptide.test.pp,
#   aes(
#     x = Pred.nn,
#     y = Fitness.nb
#   )
# ) +
#   geom_point()

# Comparing against my non-machine learning model.
trainindex.nocluster <- createDataPartition(
  peptide.data$Fitness.nb, p = 0.8, list = FALSE, times = 1
)
peptide.mixed.nb.lm <- lmer(
  data = peptide.data[trainindex.nocluster,],
  formula = log(Fitness.nb) ~
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys + Arg +
    #Clustering.Six +
    #WaltzBinary +
    #net.charge +
    (1|Cluster) +
    0
)
lmer.pred <- predict(peptide.mixed.nb.lm, newdata = peptide.data, re.form = NA, type = "response")
pred.lme.df <- tibble(
  "fitness" = peptide.data$Fitness.nb,
  "predicted" = lmer.pred,
  "weight" = peptide.data$Weight.nb,
  "cluster" = peptide.data$Cluster
)
pred.lme.df
pred.lme.bycluster <-
  pred.lme.df %>%
  group_by(cluster) %>%
  filter(weight == max(weight))
pred.lme.bycluster
cor.test(log(pred.lme.bycluster[-trainindex,]$fitness),
         pred.lme.bycluster[-trainindex,]$predicted, method = "pearson")
summary(lm(data = pred.lme.bycluster[-trainindex,],
           formula = log(fitness) ~ predicted))
# Seems good, but I can't actually compare this. The training index does not match the test index.

# Retrying with a fixed effects model.
peptide.all.maxweight <- 
  peptide.data %>%
  group_by(Cluster) %>%
  filter(Weight.nb == max(Weight.nb))
peptide.maxweight.lm <- lm(
  data = peptide.all.maxweight[trainindex,],
  formula = log(Fitness.nb) ~
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys + Arg +
    #Clustering.Six +
    #WaltzBinary +
    #net.charge +
    0
)
summary(peptide.maxweight.lm)
fit.lm.pred <- predict(peptide.maxweight.lm, newdata = peptide.all.maxweight[-trainindex,])
fit.lm.pred
cor.test(fit.lm.pred, log(peptide.all.maxweight[-trainindex,]$Fitness.nb), method = "pearson")
summary(lm(formula = log(peptide.all.maxweight[-trainindex,]$Fitness.nb) ~ fit.lm.pred,
           weights = peptide.all.maxweight[-trainindex,]$Weight.nb))


# I think the problem is that there is overfitting going on.
# Let's see if fewer predictors results in better predictions.
source(file = "Scripts/aa_comp_metrics.R")
peptide.maxweight.rsa <- peptide.clusters.maxweight
peptide.maxweight.rsa


peptide.maxweight.rsa$Pos1.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos1, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos2.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos2, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos3.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos3, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos4.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos4, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos5.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos5, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos6.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos6, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos7.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos7, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos8.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos8, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos9.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos9, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos10.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos10, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos11.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos11, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos12.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos12, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos13.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos13, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos14.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos14, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos15.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos15, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos16.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos16, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos17.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos17, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos18.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos18, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos19.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos19, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos20.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos20, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos21.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos21, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos22.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos22, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos23.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos23, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos24.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos24, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos25.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos25, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos26.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos26, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos27.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos27, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos28.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos28, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos29.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos29, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos30.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos30, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos31.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos31, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos32.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos32, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos33.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos33, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos34.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos34, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos35.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos35, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos36.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos36, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos37.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos37, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos38.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos38, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos39.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos39, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos40.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos40, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos41.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos41, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos42.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos42, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos43.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos43, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos44.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos44, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos45.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos45, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos46.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos46, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos47.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos47, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos48.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos48, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos49.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos49, metric = "RSA-hydrophilicity")
peptide.maxweight.rsa$Pos50.RSA <- mean.metric.calculator(peptide.maxweight.rsa$Pos50, metric = "RSA-hydrophilicity")

peptide.maxweight.rsa[,c(1:5, 56:105)]
peptide.maxweight.rsa <- peptide.maxweight.rsa[,c(1:5, 56:105)]

# New training index.
# trainindex.rsa <- createDataPartition(
#   peptide.maxweight.rsa$Fitness.nb, p = 0.8, list = FALSE, times = 1
# )

# Splitting into test and training.
rsa.train <- peptide.maxweight.rsa[trainindex,]
rsa.test <- peptide.maxweight.rsa[-trainindex,]

# Pre-processing
pp.rsa <- preProcess(rsa.train,
                     method = c("center", "scale"),
                     outcome = rsa.train$Fitness.nb)

rsa.train.pp <- predict(pp.rsa, rsa.train)
rsa.test.pp <- predict(pp.rsa, rsa.test)
glimpse(rsa.train.pp)

# Train using neuralnet. Should be the same as above, hopefully!
rsa.nn <- neuralnet(
  Fitness.nb ~ .,
  data = rsa.train.pp[,-c(1,3:5)],
  startweights = rsa.train$Weight.nb,
  hidden = 1
)
rsa.nn
rsa.nn.pred <- predict(rsa.nn, rsa.test.pp)
rsa.nn.pred
cor.test(rsa.nn.pred, rsa.test.pp$Fitness.nb, method = "pearson")
summary(lm(formula = rsa.test.pp$Fitness.nb ~ rsa.nn.pred, weights = rsa.test$Weight.nb))

# Train using just the 20 AA counts.
peptide.aa.counts <- dplyr::select(peptide.data,
                                   PeptideID,
                                   Fitness.nb,
                                   Weight.nb,
                                   Cluster,
                                   AASeq,
                                   Leu, Pro, Met, Trp, Ala,
                                   Val, Phe, Ile, Gly, Ser,
                                   Thr, Cys, Asn, Gln, Tyr,
                                   His, Asp, Glu, Lys, Arg,
                                   Clustering.Six, WaltzBinary, net.charge)
peptide.clusters.aacounts <- 
  peptide.aa.counts %>%
  group_by(Cluster) %>%
  filter(Weight.nb == max(Weight.nb))

peptide.clusters.aacounts

# Test and training.
aa.counts.train <- peptide.clusters.aacounts[trainindex,]
aa.counts.test <- peptide.clusters.aacounts[-trainindex,]

# Pre-processing.
pp.aa.counts <- preProcess(aa.counts.train,
                     method = c("center", "scale"),
                     outcome = aa.counts.train$Fitness.nb)

aa.counts.train.pp <- predict(pp.aa.counts, aa.counts.train)
aa.counts.test.pp <- predict(pp.aa.counts, aa.counts.test)
glimpse(aa.counts.train.pp)

start.weights <- coef(peptide.maxweight.lm)
as.vector(start.weights)

aa.counts.nn <- neuralnet(
  Fitness.nb ~ Leu + Pro + Met + Trp + Ala + Val + Phe + Ile + Gly + Ser +
               Thr + Cys + Asn + Gln + Tyr + His + Asp + Glu + Lys + Arg
               #+ Clustering.Six + WaltzBinary + net.charge
  ,
  data = aa.counts.train.pp,
  startweights = as.vector(start.weights),
  #startweights = NULL,
  rep = 1,
  hidden = c(10,2),
  lifesign = "full",
  #algorithm = "backprop",
  #learningrate = 0.1
  stepmax = 1e+05,
  threshold = 0.04
)
#plot(aa.counts.nn)
aa.counts.nn.pred <- predict(aa.counts.nn, aa.counts.test.pp)
aa.counts.nn.pred
cor.test(aa.counts.nn.pred, aa.counts.test.pp$Fitness.nb, method = "pearson")
summary(lm(formula = aa.counts.test.pp$Fitness.nb ~ aa.counts.nn.pred))
