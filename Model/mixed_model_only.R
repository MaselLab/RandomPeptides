# Model script only.

# Packages.
library(lme4)
library(stringr)
library(Hmisc)

# Load peptide data.
peptide.data <- read.table(file = "Scripts/RandomPeptides/Data/supplemental_table_1.tsv", header = T, stringsAsFactors = F)
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
    WaltzBinary +
    net.charge +
    #sqrt(ISD.iupred2) +
    #TangoBinary +
    #CamSol.avg +
    #PredHel +
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
drop1(peptide.mixed.nb.freq.aaonly.lm, test = "Chisq")

# Non-freq models.
peptide.mixed.nb.aaonly.lm <- lmer(
  data = peptide.data,
  formula = log(Fitness.nb) ~
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys + Arg +
    (1|Cluster) +
    0,
  weights = Weight.nb
)
summary(peptide.mixed.nb.aaonly.lm)
drop1(peptide.mixed.nb.aaonly.lm, test = "Chisq")

peptide.mixed.nb.lm <- lmer(
  data = peptide.data,
  formula = log(Fitness.nb) ~
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys + Arg +
    Clustering.Six +
    WaltzBinary +
    net.charge +
    (1|Cluster) +
    0,
  weights = Weight.nb
)
summary(peptide.mixed.nb.lm)
drop1(peptide.mixed.nb.lm, test = "Chisq")

peptide.mixed.other.predictors.nb.lm <- lmer(
  data = peptide.data,
  formula = log(Fitness.nb) ~
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys + Arg +
    Clustering.Six +
    WaltzBinary +
    net.charge +
    (1|Cluster) +
    0,
  weights = Weight.nb
)
summary(peptide.mixed.other.predictors.nb.lm)
anova(peptide.mixed.nb.aaonly.lm, peptide.mixed.other.predictors.nb.lm, test = "LRT")
AIC(peptide.mixed.nb.aaonly.lm)
AIC(peptide.mixed.other.predictors.nb.lm)