# Compaing the predictive power of GC content to that of amino acid composition.

# Packages.
library(lme4)
library(tidyverse)
library(stringr)
library(Hmisc)

# Load peptide data.
peptide.data <- read.table(file = "Scripts/RandomPeptides/Data/supplemental_table_1.tsv", header = T, stringsAsFactors = F)
peptide.data

# Comparing GC content as a predictor to amino acid composition.
peptide.mixed.intercept <- lmer(
  data = peptide.data,
  formula = Fitness.nb ~ 1 + (1|Cluster),
  weights = Weight.nb.5.7
)

peptide.mixed.gc <- lmer(
  data = peptide.data[!is.na(peptide.data$GC.avg),],
  formula = Fitness.nb ~ GC.avg +
    (1|Cluster),
  weights = Weight.nb.5.7
)
summary(peptide.mixed.gc)
peptide.mixed.aaonly.gc <- lmer(
  data = peptide.data[!is.na(peptide.data$GC.avg),],
  # contrasts = c("Leu" = -0.045, "Pro" = 0.068, "Met" = -0.084, "Trp" = -0.0013, "Ala" = 0.059,
  #           "Val" = -0.051, "Phe" = -0.11, "Ile" = -0.16, "Gly" = 0.017, "Ser" = 0.035,
  #           "Thr" = -0.0094, "Cys" = -0.024, "Asn" = -0.047, "Gln" = 0.0053, "Tyr" = -0.077,
  #           "His" = -0.075, "Asp" = 0.0093, "Glu" = -0.029, "Lys" = -0.09, "Arg" = -0.043,
  #           "GC.avg" = 10),
  formula = Fitness.nb ~
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys + Arg +
    GC.avg +
    (1|Cluster) +
    0,
  weights = Weight.nb.5.7
)
summary(peptide.mixed.aaonly.gc)
AIC(peptide.mixed.gc)
AIC(peptide.mixed.aaonly.gc)
anova(peptide.mixed.aaonly.gc, peptide.mixed.gc, test = "LRT")
drop1(peptide.mixed.aaonly.gc, test = "Chisq")

peptide.mixed.nb.aaonly.lm <- lmer(
  data = peptide.data[!is.na(peptide.data$GC.avg),],
  formula = Fitness.nb ~
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys + Arg +
    (1|Cluster) +
    0,
  weights = Weight.nb.5.7
)
summary(peptide.mixed.nb.aaonly.lm)
drop1(peptide.mixed.nb.aaonly.lm, test = "Chisq")
anova(peptide.mixed.aaonly.gc, peptide.mixed.nb.aaonly.lm, test = "LRT")
#anova(peptide.mixed.nb.aaonly.lm, peptide.mixed.intercept, test = "LRT")

# Checking LOGOCV.
source(file = "Scripts/RandomPeptides/Metrics/logocv.R")

wmse.gc <- logo_cv(
  data.df = peptide.data,
  predictors = "GC.avg",
  dependent.variable = "Fitness.nb",
  group.col = "Cluster",
  weight.col = "Weight.nb.5.7",
  re.form = NA, random.only = F, type = "response"
)
wmse.gc

wmse.aa_plus_gc <- logo_cv(
  data.df = peptide.data,
  predictors = c("Leu", "Pro", "Met", "Trp", "Ala",
                 "Val", "Phe", "Ile", "Gly", "Ser",
                 "Thr", "Cys", "Asn", "Gln", "Tyr",
                 "His", "Asp", "Glu", "Lys",
                 "GC.avg"),
  dependent.variable = "Fitness.nb",
  group.col = "Cluster",
  weight.col = "Weight.nb.5.7",
  re.form = NA, random.only = F, type = "response"
)
wmse.aa_plus_gc

# Calculating disorder propensity for a disorder + GC content model.
# Calculating disorder propensity.
source("~/MaselLab/RandomPeptides/Scripts/RandomPeptides/Metrics/aa_comp_metrics.R")
peptide.data$disorder <- mean.metric.calculator(
  aa.sequence = peptide.data$AASeq,
  metric = "disorder"
)

# Checking the model.
peptide.mixed.gc.disorder <- lmer(
  data = peptide.data[!is.na(peptide.data$GC.avg),],
  formula = Fitness.nb ~ GC.avg + disorder +
    (1|Cluster),
  weights = Weight.nb.5.7
)
summary(peptide.mixed.gc.disorder)
drop1(peptide.mixed.gc.disorder, test = "Chisq")

wmse.disprop_plus_gc <- logo_cv(
  data.df = peptide.data,
  predictors = c("disorder", "GC.avg"),
  dependent.variable = "Fitness.nb",
  group.col = "Cluster",
  weight.col = "Weight.nb.5.7",
  re.form = NA, random.only = F, type = "response"
)
wmse.disprop_plus_gc

# Checking lmer starting values.
aaonly.lm.ff <- lFormula(data = peptide.data[!is.na(peptide.data$GC.avg),], Fitness.nb ~
                           Leu + Pro + Met + Trp + Ala +
                           Val + Phe + Ile + Gly + Ser +
                           Thr + Cys + Asn + Gln + Tyr +
                           His + Asp + Glu + Lys + Arg +
                           (1|Cluster) +
                           0)
ifelse(aaonly.lm.ff$reTrms$lower == 0, 1, 0)

# Do single amino acids improve a GC content only model?
peptide.mixed.aa.gc <- lmer(
  data = peptide.data[!is.na(peptide.data$GC.avg),],
  formula = Fitness.nb ~ GC.avg + Pro +
    (1|Cluster),
  weights = Weight.nb.5.7
)
summary(peptide.mixed.aa.gc)
drop1(peptide.mixed.aa.gc, test = "Chisq")

# And now checking weighted R^2 for fitness ~ predicted fitness.
peptide.data$fit.aa <- predict(peptide.mixed.nb.aaonly.lm,
                               newdata = peptide.data,
                               type = "response",
                               random.only = F,
                               re.form = NA)
peptide.data$fit.gc <- predict(peptide.mixed.gc,
                               newdata = peptide.data,
                               type = "response",
                               random.only = F,
                               re.form = NA)
peptide.data$fit.aa.gc <- predict(peptide.mixed.aa.gc,
                                  newdata = peptide.data,
                                  type = "response",
                                  random.only = F,
                                  re.form = NA)
peptide.data[is.na(peptide.data$fit.gc), "PeptideID"]
peptide.data[is.na(peptide.data$fit.aa.gc), "PeptideID"]

peptide.cluster <- 
  peptide.data %>%
  group_by(Cluster) %>%
  summarise(Weight.nb.sum = sum(Weight.nb.5.7), Fitness.nb = wtd.mean(Fitness.nb, weights = Weight.nb.5.7),
            fit.aa = wtd.mean(fit.aa, weights = Weight.nb.5.7), fit.gc = wtd.mean(fit.gc, weights = Weight.nb.5.7, na.rm = T),
            fit.aa.gc = wtd.mean(fit.aa.gc, weights = Weight.nb.5.7, na.rm = T),
            Leu = wtd.mean(Leu, weights = Weight.nb.5.7), Phe = wtd.mean(Phe, weights = Weight.nb.5.7),
            Met = wtd.mean(Met, weights = Weight.nb.5.7), Val = wtd.mean(Val, weights = Weight.nb.5.7),
            Ile = wtd.mean(Ile, weights = Weight.nb.5.7), Lys = wtd.mean(Lys, weights = Weight.nb.5.7),
            His = wtd.mean(His, weights = Weight.nb.5.7), Arg = wtd.mean(Arg, weights = Weight.nb.5.7),
            Glu = wtd.mean(Glu, weights = Weight.nb.5.7), Asp = wtd.mean(Asp, weights = Weight.nb.5.7),
            Gln = wtd.mean(Gln, weights = Weight.nb.5.7), Asn = wtd.mean(Asn, weights = Weight.nb.5.7),
            Gly = wtd.mean(Gly, weights = Weight.nb.5.7), Ala = wtd.mean(Ala, weights = Weight.nb.5.7),
            Pro = wtd.mean(Pro, weights = Weight.nb.5.7), Ser = wtd.mean(Ser, weights = Weight.nb.5.7),
            Trp = wtd.mean(Trp, weights = Weight.nb.5.7), Tyr = wtd.mean(Tyr, weights = Weight.nb.5.7),
            Thr = wtd.mean(Thr, weights = Weight.nb.5.7), Cys = wtd.mean(Cys, weights = Weight.nb.5.7),
            GC.avg = wtd.mean(GC.avg, weights = Weight.nb.5.7, na.rm = T),
            Clustering.Six = wtd.mean(Clustering.Six, weights = Weight.nb.5.7),
            net.charge = wtd.mean(net.charge, weights = Weight.nb.5.7),
            TangoAAsInAPRs = wtd.mean(TangoAAsInAPRs, weights = Weight.nb.5.7))
peptide.cluster

# Predicting fitness with predicted fitness to get R^2 values.
predfit.aa <- lm(
  data = peptide.cluster,
  formula = Fitness.nb ~ fit.aa,
  weights = Weight.nb.sum
)
summary(predfit.aa)

predfit.gc <- lm(
  data = peptide.cluster,
  formula = Fitness.nb ~ fit.gc,
  weights = Weight.nb.sum
)
summary(predfit.gc)

predfit.aa.gc <- lm(
  data = peptide.cluster,
  formula = Fitness.nb ~ fit.aa.gc,
  weights = Weight.nb.sum
)
summary(predfit.aa.gc)

# Checking just fixed effects models. This makes the comparison much simpler, and I don't need to
# worry about the fact that the mixed models are fitted numerically and might be stuck in some
# local optimum.
aa.lm <- lm(
  data = peptide.cluster,
  formula = Fitness.nb ~ 
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys #+ Arg +
  #0
  ,
  weights = Weight.nb.sum
)
aa.lm.summary <- summary(aa.lm)
aa.lm.summary

gc.lm <- lm(
  data = peptide.cluster,
  formula = Fitness.nb ~ GC.avg,
  weights = Weight.nb.sum
)
gc.lm.summary <- summary(gc.lm)
gc.lm.summary

aa.gc.lm <- lm(
  data = peptide.cluster,
  formula = Fitness.nb ~ GC.avg +
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys #+ Arg +
  #0
  ,
  weights = Weight.nb.sum
)
aa.gc.lm.summary <- summary(aa.gc.lm)
aa.gc.lm.summary

# full.lm <- lm(
#   data = peptide.cluster,
#   formula = Fitness.nb ~ 
#     Leu + Pro + Met + Trp + Ala +
#     Val + Phe + Ile + Gly + Ser +
#     Thr + Cys + Asn + Gln + Tyr +
#     His + Asp + Glu + Lys #+ Arg +
#   #0
#   + TangoAAsInAPRs + net.charge + Clustering.Six
#   ,
#   weights = Weight.nb.sum
# )
# full.lm.summary <- summary(full.lm)
# full.lm.summary

# full.gc.lm <- lm(
#   data = peptide.cluster,
#   formula = Fitness.nb ~ 
#     Leu + Pro + Met + Trp + Ala +
#     Val + Phe + Ile + Gly + Ser +
#     Thr + Cys + Asn + Gln + Tyr +
#     His + Asp + Glu + Lys #+ Arg +
#   #0
#   + TangoAAsInAPRs + net.charge + Clustering.Six
#   + GC.avg
#   ,
#   weights = Weight.nb.sum
# )
# full.gc.lm.summary <- summary(full.gc.lm)
# full.gc.lm.summary
# drop1(full.gc.lm, test = "Chisq")

# Comparing the GC content only, aa comp only, and GC content + aa comp models to intercept only.
intercept.lm <- lm(
  data = peptide.cluster,
  formula = Fitness.nb ~ 1,
  weights = Weight.nb.sum
)
summary(intercept.lm)

aa.lrt <- anova(intercept.lm, aa.lm, test = "LRT")
aa.lrt
aa.lrt$`Pr(>Chi)`
gc.lrt <- anova(intercept.lm, gc.lm, test = "LRT")
gc.lrt$`Pr(>Chi)`
aa.gc.lrt <- anova(intercept.lm, aa.gc.lm, test = "LRT")
aa.gc.lrt$`Pr(>Chi)`

# Seeing if GC content improves the aa comp only model, and vice versa.
aa.gc.gc.lrt <- anova(gc.lm, aa.gc.lm, test = "LRT")
aa.gc.gc.lrt
aa.gc.aa.lrt <- anova(aa.lm, aa.gc.lm, test = "LRT")
aa.gc.aa.lrt
