# Test vs training for various models.

# Packages.
library(lme4)
library(stringr)
library(Hmisc)
library(tidyverse)

# Load peptide data.
peptide.data <- read.table(file = "Scripts/RandomPeptides/Data/supplemental_table_1.tsv", header = T, stringsAsFactors = F)
peptide.data

# Calculating disorder propensity.
source("~/MaselLab/RandomPeptides/Scripts/RandomPeptides/Metrics/aa_comp_metrics.R")
peptide.data$disorder <- mean.metric.calculator(
  aa.sequence = peptide.data$AASeq,
  metric = "disorder"
)
hist(peptide.data$disorder)

# Set seed.
set.seed(55588)

## 75% of the sample size
sample_size <- floor(0.90 * max(peptide.data$Cluster))

## set the seed to make your partition reproducible
train_ind <- sample(seq_len(max(peptide.data$Cluster)), size = sample_size)

train <- peptide.data[peptide.data$Cluster %in% train_ind, ]
test <- peptide.data[!(peptide.data$Cluster %in% train_ind), ]

# Checking results per cluster, not per data point.
test.maxw <- test %>% group_by(Cluster) %>% filter(Weight.nb.5.7 == max(Weight.nb.5.7))

# Non-freq models.
peptide.lmm <- lmer(
  data = train,
  formula = Fitness.nb ~
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys + Arg +
    (1|Cluster) +
    0,
  weights = Weight.nb.5.7
)
summary(peptide.lmm)
AIC(peptide.lmm)
drop1(peptide.lmm, test = "Chisq")

pred.lmm.peptide <- predict(peptide.lmm, test.maxw,
                        re.form = NA, random.only = F, type = "response")

# ISD only.
# Non-freq models.
isd.lmm <- lmer(
  data = train,
  formula = Fitness.nb ~
    sqrt(ISD.iupred2) +
    (1|Cluster),
  weights = Weight.nb.5.7
)
summary(isd.lmm)
AIC(isd.lmm)
drop1(isd.lmm, test = "Chisq")

pred.lmm.isd <- predict(isd.lmm, test.maxw,
                        re.form = NA, random.only = F, type = "response")

# Disorder only.
disorder.lmm <- lmer(
  data = train,
  formula = Fitness.nb ~
    disorder +
    (1|Cluster),
  weights = Weight.nb.5.7
)
summary(disorder.lmm)
AIC(disorder.lmm)
drop1(disorder.lmm, test = "Chisq")

pred.lmm.disorder <- predict(disorder.lmm, test.maxw,
                        re.form = NA, random.only = F, type = "response")

# CamSol only.
# Non-freq models.
camsol.lmm <- lmer(
  data = train,
  formula = Fitness.nb ~
    CamSol.avg +
    (1|Cluster),
  weights = Weight.nb.5.7
)
summary(camsol.lmm)
AIC(camsol.lmm)
drop1(camsol.lmm, test = "Chisq")

pred.lmm.camsol <- predict(camsol.lmm, test.maxw,
                    re.form = NA, random.only = F, type = "response")

# clustering only.
# Non-freq models.
clustering.lmm <- lmer(
  data = train,
  formula = Fitness.nb ~
    Clustering.Six +
    (1|Cluster),
  weights = Weight.nb.5.7
)
summary(clustering.lmm)
AIC(clustering.lmm)
drop1(clustering.lmm, test = "Chisq")

pred.lmm.clustering <- predict(clustering.lmm, test.maxw,
                       re.form = NA, random.only = F, type = "response")

# Tango only.
# Non-freq models.
Tango.lmm <- lmer(
  data = train,
  formula = Fitness.nb ~
    TangoAAsInAPRs +
    (1|Cluster),
  weights = Weight.nb.5.7
)
summary(Tango.lmm)
AIC(Tango.lmm)
drop1(Tango.lmm, test = "Chisq")

pred.lmm.Tango <- predict(Tango.lmm, test.maxw,
                           re.form = NA, random.only = F, type = "response")

# Results
cor(pred.lmm.peptide, test.maxw$Fitness.nb)
cor(pred.lmm.isd, test.maxw$Fitness.nb)
cor(pred.lmm.disorder, test.maxw$Fitness.nb)
cor(pred.lmm.camsol, test.maxw$Fitness.nb)
cor(pred.lmm.clustering, test.maxw$Fitness.nb)
cor(pred.lmm.Tango, test.maxw$Fitness.nb)

# Squared correlation coefficient (not really R^2, but something of an approximation.)
cor(pred.lmm.peptide, test.maxw$Fitness.nb) ^ 2
cor(pred.lmm.isd, test.maxw$Fitness.nb) ^ 2
cor(pred.lmm.disorder, test.maxw$Fitness.nb) ^ 2
cor(pred.lmm.camsol, test.maxw$Fitness.nb) ^ 2
cor(pred.lmm.clustering, test.maxw$Fitness.nb) ^ 2
cor(pred.lmm.Tango, test.maxw$Fitness.nb) ^ 2

# Mean squared error.
(1/nrow(test.maxw)) * sum((test.maxw$Fitness.nb - pred.lmm.peptide)^2)
(1/nrow(test.maxw)) * sum((test.maxw$Fitness.nb - pred.lmm.isd)^2)
(1/nrow(test.maxw)) * sum((test.maxw$Fitness.nb - pred.lmm.disorder)^2)
(1/nrow(test.maxw)) * sum((test.maxw$Fitness.nb - pred.lmm.camsol)^2)
(1/nrow(test.maxw)) * sum((test.maxw$Fitness.nb - pred.lmm.clustering)^2)
(1/nrow(test.maxw)) * sum((test.maxw$Fitness.nb - pred.lmm.Tango)^2)

# Leave one group out cross validation.
source("~/MaselLab/RandomPeptides/Scripts/logocv.R")
wmse.aa <- logo_cv(
  data.df = peptide.data,
  predictors = c("Leu", "Pro", "Met", "Trp", "Ala",
                   "Val", "Phe", "Ile", "Gly", "Ser",
                   "Thr", "Cys", "Asn", "Gln", "Tyr",
                   "His", "Asp", "Glu", "Lys"),
  dependent.variable = "Fitness.nb",
  group.col = "Cluster",
  weight.col = "Weight.nb.5.7",
  re.form = NA, random.only = F, type = "response"
)

# Checking the rest of the predictors.
peptide.data$ISD.sqrt <- sqrt(peptide.data$ISD.iupred2)
pred.list <- c("ISD.sqrt", "disorder", "CamSol.avg", "TangoAAsInAPRs", "Clustering.Six")
wmse.list <-
  lapply(pred.list, logo_cv,
         data.df = peptide.data,
         dependent.variable = "Fitness.nb",
         group.col = "Cluster",
         weight.col = "Weight.nb.5.7",
         re.form = NA, random.only = F, type = "response")

tibble(
  "Predictor" = c("AA", pred.list),
  "WMSE" = as.numeric(c(wmse.aa, wmse.list))
)
