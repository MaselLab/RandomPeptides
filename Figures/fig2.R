# Figure 2 Script.

# Packages.
library(lme4)
library(stringr)
library(tidyverse)
library(wCorr)
library(Hmisc)

# Load peptide data.
peptide.data <- read.table(file = "Scripts/RandomPeptides/Data/supplemental_table_1.tsv", header = T, stringsAsFactors = F)

# Full model.
fitness.nb.full.lm <- lmer(data = peptide.data,
                           formula = log(Fitness.nb) ~
                             Leu + Pro + Met + Trp + Ala +
                             Val + Phe + Ile + Gly + Ser +
                             Thr + Cys + Asn + Gln + Tyr +
                             His + Asp + Glu + Lys + Arg +
                             Clustering.Six +
                             WaltzBinary +
                             net.charge +
                             (1|Cluster) + 0,
                           weights = Weight.nb)
fitness.nb.full.summary <- summary(fitness.nb.full.lm)

peptide.data$fit.full <- 
  predict(
    fitness.nb.full.lm,
    newdata = peptide.data,
    re.form = NA
  )

# Predicted estimated fitness with predicted fitness from the full model.
fitness.pred.fit.cluster.lm <- lmer(
  data = peptide.data,
  formula = log(Fitness.nb) ~ fit.full + (1|Cluster),
  weights = Weight.nb
)
summary(fitness.pred.fit.cluster.lm)
#pred.full.summary <- summary(fitness.pred.fit.cluster.lm)
#pred.full.summary$coefficients

# Amino acid composition only.
fitness.nb.aa.lm <- lmer(data = peptide.data,
                           formula = log(Fitness.nb) ~
                             Leu + Pro + Met + Trp + Ala +
                             Val + Phe + Ile + Gly + Ser +
                             Thr + Cys + Asn + Gln + Tyr +
                             His + Asp + Glu + Lys + Arg +
                             (1|Cluster) + 0,
                           weights = Weight.nb)
fitness.nb.aa.summary <- summary(fitness.nb.aa.lm)

peptide.data$fit.aa <- 
  predict(
    fitness.nb.aa.lm,
    newdata = peptide.data,
    re.form = NA
  )

# Predicted estimated fitness with predicted fitness from the aa comp only model.
fitness.pred.aa.fit.cluster.lm <- lmer(
  data = peptide.data,
  formula = log(Fitness.nb) ~ predict(fitness.nb.aa.lm, newdata = peptide.data, re.form = NA) + (1|Cluster),
  weights = Weight.nb
)
summary(fitness.pred.aa.fit.cluster.lm)
#pred.aa.summary <- summary(fitness.pred.aa.fit.cluster.lm)
#pred.aa.summary$coefficients

# Combining the data by cluster for plotting.
by_cluster <-
  peptide.data %>% 
  group_by(Cluster) %>%
  summarise(Weight.nb.sum = sum(Weight.nb), 
            ISD.iupred2 = wtd.mean(ISD.iupred2, weights = Weight.nb), Fitness.nb = wtd.mean(Fitness.nb, weights = Weight.nb),
            Leu = wtd.mean(Leu, weights = Weight.nb), Phe = wtd.mean(Phe, weights = Weight.nb),
            Met = wtd.mean(Met, weights = Weight.nb), Val = wtd.mean(Val, weights = Weight.nb),
            Ile = wtd.mean(Ile, weights = Weight.nb), Lys = wtd.mean(Lys, weights = Weight.nb),
            His = wtd.mean(His, weights = Weight.nb), Arg = wtd.mean(Arg, weights = Weight.nb),
            Glu = wtd.mean(Glu, weights = Weight.nb), Asp = wtd.mean(Asp, weights = Weight.nb),
            Gln = wtd.mean(Gln, weights = Weight.nb), Asn = wtd.mean(Asn, weights = Weight.nb),
            Gly = wtd.mean(Gly, weights = Weight.nb), Ala = wtd.mean(Ala, weights = Weight.nb),
            Pro = wtd.mean(Pro, weights = Weight.nb), Ser = wtd.mean(Ser, weights = Weight.nb),
            Trp = wtd.mean(Trp, weights = Weight.nb), Tyr = wtd.mean(Tyr, weights = Weight.nb),
            Thr = wtd.mean(Thr, weights = Weight.nb), Cys = wtd.mean(Cys, weights = Weight.nb),
            Clustering.Six = wtd.mean(Clustering.Six, weights = Weight.nb),
            WaltzBinary = wtd.mean(WaltzBinary, weights = Weight.nb),
            CamSol.avg = wtd.mean(CamSol.avg, weights = Weight.nb),
            charge.pos = wtd.mean(charge.pos, weights = Weight.nb), charge.neg = wtd.mean(charge.neg),
            net.charge = wtd.mean(net.charge, weights = Weight.nb),
            fit.full = wtd.mean(fit.full, weights = Weight.nb),
            fit.aa = wtd.mean(fit.aa, weights = Weight.nb)
  )
# by_cluster$Fitness.nb.weighted <- rep(NA, length(by_cluster$Fitness.nb))
# for (i in 1:length(by_cluster$Fitness.nb)) {
#   by_cluster[by_cluster$Cluster == i, "Fitness.nb.weighted"] <- 
#     weighted.mean(peptide.data[peptide.data$Cluster == i, "Fitness.nb"],
#                   peptide.data[peptide.data$Cluster == i, "Weight.nb"])
# }
# by_cluster$fit.full.weighted <- rep(NA, length(by_cluster$CamSol.avg))
# by_cluster$fit.aa.weighted <- rep(NA, length(by_cluster$WaltzBinary))
# for (i in 1:length(by_cluster$Fitness.nb)) {
#   by_cluster[by_cluster$Cluster == i, "fit.full.weighted"] <- 
#     weighted.mean(predict(fitness.nb.full.lm,
#                           newdata = peptide.data[peptide.data$Cluster == i, ],
#                           re.form = NA),
#                   peptide.data[peptide.data$Cluster == i, "Weight.nb"])
#   by_cluster[by_cluster$Cluster == i, "fit.aa.weighted"] <- 
#     weighted.mean(predict(fitness.nb.aa.lm,
#                           newdata = peptide.data[peptide.data$Cluster == i, ],
#                           re.form = NA),
#                   peptide.data[peptide.data$Cluster == i, "Weight.nb"])
# }

# Correlation tests.
fit.full.lm <- lm(data = by_cluster,
                  formula = log(Fitness.nb) ~ fit.full,
                  weights = Weight.nb.sum)
summary(fit.full.lm)
pred.full.summary <- summary(fit.full.lm)
pred.full.summary$coefficients

fit.aa.lm <- lm(data = by_cluster,
                formula = log(Fitness.nb) ~ fit.aa,
                weights = Weight.nb.sum)
summary(fit.aa.lm)
pred.aa.summary <- summary(fit.aa.lm)
pred.aa.summary$coefficients

with(by_cluster, weightedCorr(log(Fitness.nb), fit.full, method = "Pearson", weights = Weight.nb.sum))
with(by_cluster, weightedCorr(log(Fitness.nb), fit.aa, method = "Pearson", weights = Weight.nb.sum))

# Trying to flip the regression to X~Y, seeing what happens.
flipped.fit.full.lm <- lm(data = by_cluster,
                  formula = fit.full ~ log(Fitness.nb),
                  weights = Weight.nb.sum)
summary(flipped.fit.full.lm)
flipped.pred.full.summary <- summary(flipped.fit.full.lm)
flipped.pred.full.summary$coefficients
flipped.full.int <- -flipped.pred.full.summary$coefficients[1,1] / flipped.pred.full.summary$coefficients[2,1]
flipped.full.beta <- 1/flipped.pred.full.summary$coefficients[2,1]

flipped.fit.aa.lm <- lm(data = by_cluster,
                formula = fit.aa ~ log(Fitness.nb),
                weights = Weight.nb.sum)
summary(flipped.fit.aa.lm)
flipped.pred.aa.summary <- summary(flipped.fit.aa.lm)
flipped.pred.aa.summary$coefficients
flipped.aa.int <- -flipped.pred.aa.summary$coefficients[1,1] / flipped.pred.aa.summary$coefficients[2,1]
flipped.aa.beta <- 1/flipped.pred.aa.summary$coefficients[2,1]

# Plotting part A.
todays.date <- "4-9-20"
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
png(filename = paste("Scripts/Figures/fitness_pred_full_", todays.date, ".png", sep = ""),
    height = 500, width = 500)
ggplot(data = by_cluster,
       aes(x = fit.full,
           y = log(Fitness.nb),
           size = Weight.nb.sum,
           weight = Weight.nb.sum)
) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = cbbPalette[2], size = 1.5) +
  geom_smooth(method = "lm", color = cbbPalette[6], size = 1.5, se = F) +
  #stat_function(fun = function(x)flipped.full.int + flipped.full.beta*x,
  #              geom = "line", color = cbbPalette[4], size = 1.5) +
  ylab("Fitness") +
  xlab("Predicted fitness") +
  scale_y_continuous(breaks = log(c(0.2, 0.5, 1, 2, 10)),
                     labels = c(0.2, 0.5, 1, 2, 10)) +
  scale_x_continuous(breaks = log(c(0.2, 0.5, 1)),
                     labels = c(0.2, 0.5, 1)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

# Plotting part B.
png(filename = paste("Scripts/Figures/fitness_pred_aacomp_", todays.date, ".png", sep = ""),
    height = 500, width = 500)
ggplot(data = by_cluster,
       aes(y = fit.aa,
           x = log(Fitness.nb),
           size = Weight.nb.sum,
           weight = Weight.nb.sum)
) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = cbbPalette[2], size = 1.5) +
  stat_function(fun = function(x)pred.aa.summary$coefficients[1,1] + pred.aa.summary$coefficients[2,1]*x,
                geom = "line", color = cbbPalette[6], size = 1.5) +
  ylab("Fitness") +
  xlab("Predicted fitness") +
  scale_y_continuous(breaks = log(c(0.2, 0.5, 1, 2, 10)),
                     labels = c(0.2, 0.5, 1, 2, 10)) +
  scale_x_continuous(breaks = log(c(0.2, 0.5, 1)),
                     labels = c(0.2, 0.5, 1)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()
