# Model script only.

# Packages.
library(lme4)
library(stringr)

# Load peptide data.
peptide.data <- read.table(file = "Data/peptide_data_clusters_7-20-19.tsv", header = T, stringsAsFactors = F)
peptide.data

# Building the model.
# First, the model with Clustering and Waltz scores.
peptide.mixed.nb.freq.lm <- lmer(
  data = peptide.data,
  formula = log(Fitness.nb) ~
    Leu.freq + Pro.freq + Met.freq + Trp.freq + Ala.freq +
    Val.freq + Phe.freq + Ile.freq + Gly.freq + Ser.freq +
    Thr.freq + Cys.freq + Asn.freq + Gln.freq + Tyr.freq +
    His.freq + Asp.freq + Glu.freq + Lys.freq + Arg.freq +
    Clustering.Six +
    #WaltzBinary +
    #sqrt(ISD) +
    #TangoBinary +
    #CamSol.avg +
    (1|Cluster) +
    0,
  weights = Weight.nb
)
summary(peptide.mixed.nb.freq.lm)
drop1(peptide.mixed.nb.freq.lm, test = "Chisq")

# And now the model with amino acid frequencies as the only predictors.
peptide.mixed.nb.freq.aaonly.lm <- lmer(
  data = peptide.data,
  formula = log(Fitness.nb) ~
    Leu.freq + Pro.freq + Met.freq + Trp.freq + Ala.freq +
    Val.freq + Phe.freq + Ile.freq + Gly.freq + Ser.freq +
    Thr.freq + Cys.freq + Asn.freq + Gln.freq + Tyr.freq +
    His.freq + Asp.freq + Glu.freq + Lys.freq + Arg.freq +
    (1|Cluster) +
    0,
  weights = Weight.nb
)
summary(peptide.mixed.nb.freq.aaonly.lm)