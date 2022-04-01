# Model script only.

# Packages.
library(lme4)
library(stringr)
library(Hmisc)
library(tidyverse)

# Load peptide data.
peptide.data <- read.table(file = "Scripts/RandomPeptides/Data/supplemental_dataset_1.tsv", header = T, stringsAsFactors = F)
peptide.data

# Calculating disorder propensity.
source("~/MaselLab/RandomPeptides/Scripts/RandomPeptides/Metrics/aa_comp_metrics.R")
peptide.data$disorder <- mean.metric.calculator(
  aa.sequence = peptide.data$AASeq,
  metric = "disorder"
)
hist(peptide.data$disorder)

# Collapsing clusters into pseudo-datapoints.
# First, doing the mean.
cluster_means <- peptide.data %>% group_by(Cluster) %>%
  summarise(Weight.sum = sum(WEIGHT),
            ISD.iupred2 = wtd.mean(ISD.iupred2, weights = WEIGHT), FITNESS = wtd.mean(FITNESS, weights = WEIGHT),
            Leu = wtd.mean(Leu, weights = WEIGHT), Phe = wtd.mean(Phe, weights = WEIGHT),
            Met = wtd.mean(Met, weights = WEIGHT), Val = wtd.mean(Val, weights = WEIGHT),
            Ile = wtd.mean(Ile, weights = WEIGHT), Lys = wtd.mean(Lys, weights = WEIGHT),
            His = wtd.mean(His, weights = WEIGHT), Arg = wtd.mean(Arg, weights = WEIGHT),
            Glu = wtd.mean(Glu, weights = WEIGHT), Asp = wtd.mean(Asp, weights = WEIGHT),
            Gln = wtd.mean(Gln, weights = WEIGHT), Asn = wtd.mean(Asn, weights = WEIGHT),
            Gly = wtd.mean(Gly, weights = WEIGHT), Ala = wtd.mean(Ala, weights = WEIGHT),
            Pro = wtd.mean(Pro, weights = WEIGHT), Ser = wtd.mean(Ser, weights = WEIGHT),
            Trp = wtd.mean(Trp, weights = WEIGHT), Tyr = wtd.mean(Tyr, weights = WEIGHT),
            Thr = wtd.mean(Thr, weights = WEIGHT), Cys = wtd.mean(Cys, weights = WEIGHT),
            Clustering.Six = wtd.mean(Clustering.Six, weights = WEIGHT),
            TangoAAsInAPRs = wtd.mean(TangoAAsInAPRs, weights = WEIGHT),
            CamSol.avg = wtd.mean(CamSol.avg, weights = WEIGHT),
            charge.pos = wtd.mean(charge.pos, weights = WEIGHT), charge.neg = wtd.mean(charge.neg),
            net.charge = wtd.mean(net.charge, weights = WEIGHT),
            disorder = wtd.mean(disorder, weights = WEIGHT),
            ISD.delta = wtd.mean(ISD.delta, weights = WEIGHT),
            PredHel = wtd.mean(PredHel, weights = WEIGHT),
            GC.avg = wtd.mean(GC.avg, weights = WEIGHT),
            pI = wtd.mean(pI, weights = WEIGHT),
            WaltzAAsInAPRs = wtd.mean(WaltzAAsInAPRs, weights = WEIGHT)
  )

# Now picking the max weight peptide.
cluster_max <- peptide.data %>% group_by(Cluster) %>% filter(WEIGHT == max(WEIGHT))

# Non-freq models.
# Clusters collapsed using the weighted mean.
aa_only_means_lm <- lm(
  data = cluster_means,
  formula = FITNESS ~
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys + Arg +
    0,
  weights = Weight.sum
)
summary(aa_only_means_lm)
aa_summary_means <- summary(aa_only_means_lm)
drop1(aa_only_means_lm, test = "Chisq")

# Clusters collapsed using max weight.
aa_only_max_lm <- lm(
  data = cluster_max,
  formula = FITNESS ~
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys + Arg +
    0,
  weights = WEIGHT
)
summary(aa_only_max_lm)
aa_summary_max <- summary(aa_only_max_lm)
drop1(aa_only_max_lm, test = "Chisq")

# Comparing the two models.
cbind(aa_summary_means$coefficients, aa_summary_max$coefficients)

# Moving forward with the max model because it gaurantees the number of AAs sum to 50.

# Building a frequency model.
aa_freq_max_lm <- lm(
  data = cluster_max,
  formula = FITNESS ~
    Leu.freq + Pro.freq + Met.freq + Trp.freq + Ala.freq +
    Val.freq + Phe.freq + Ile.freq + Gly.freq + Ser.freq +
    Thr.freq + Cys.freq + Asn.freq + Gln.freq + Tyr.freq +
    His.freq + Asp.freq + Glu.freq + Lys.freq + Arg.freq +
    0,
  weights = WEIGHT
)
summary(aa_only_max_lm)

# Writing the beta coefficients and standard errors to file.
sorted_betas <- aa_summary_max$coefficients[base::order(aa_summary$coefficients[,1], decreasing = F),]
sorted_betas
write.csv(sorted_betas, file = paste0(
  "aa_betas_sorted_", Sys.Date(), ".csv"
))

# Checking a TMH only model.
predhel_lm <- lm(
  data = cluster_max,
  formula = FITNESS ~
    PredHel
    ,
  weights = WEIGHT
)
summary(predhel_lm)
drop1(predhel_lm, test = "Chisq")

aa_other_predictors_lm <- lm(
  data = cluster_max,
  formula = FITNESS ~
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys + Arg +
    #sqrt(ISD.iupred2) +
    #pI +
    #CamSol.avg +
    Clustering.Six +
    #TangoAAsInAPRs +
    #ISD.delta +
    #PredHel +
    #WaltzAAsInAPRs +
    #net.charge +
    #GC.avg +
    #charge.pos +
    #charge.neg +
    0
  ,
  weights = WEIGHT
)
summary(aa_other_predictors_lm)
drop1(aa_other_predictors_lm, test = "Chisq")
anova(aa_only_max_lm, aa_other_predictors_lm, test = "LRT")
AIC(aa_only_max_lm)
AIC(aa_other_predictors_lm)
