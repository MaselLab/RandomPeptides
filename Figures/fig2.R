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
                           formula = Fitness.nb ~
                             Leu + Pro + Met + Trp + Ala +
                             Val + Phe + Ile + Gly + Ser +
                             Thr + Cys + Asn + Gln + Tyr +
                             His + Asp + Glu + Lys + Arg +
                             #sqrt(ISD.iupred2) +
                             #Clustering.Six +
                             #TangoAAsInAPRs +
                             #net.charge +
                             #charge.pos +
                             #charge.neg +
                             #pI +
                             #GC.avg +
                             (1|Cluster) + 0,
                           weights = Weight.nb.5.7)
fitness.nb.full.summary <- summary(fitness.nb.full.lm)
drop1(fitness.nb.full.lm, test = "Chisq")
AIC(fitness.nb.full.lm)

peptide.data$fit.full <- 
  predict(
    fitness.nb.full.lm,
    newdata = peptide.data,
    re.form = NA
  )

# Checking effect size.
quantile(peptide.data$Fitness.nb)
quantile(peptide.data$fit.full, probs = seq(0, 1, by = 0.1))

# Predicted estimated fitness with predicted fitness from the full model.
fitness.pred.fit.cluster.lm <- lmer(
  data = peptide.data,
  formula = Fitness.nb ~ fit.full + (1|Cluster),
  weights = Weight.nb.5.7
)
summary(fitness.pred.fit.cluster.lm)
#pred.full.summary <- summary(fitness.pred.fit.cluster.lm)
#pred.full.summary$coefficients

# Amino acid composition only.
fitness.nb.aa.lm <- lmer(data = peptide.data,
                           formula = Fitness.nb ~
                             Leu + Pro + Met + Trp + Ala +
                             Val + Phe + Ile + Gly + Ser +
                             Thr + Cys + Asn + Gln + Tyr +
                             His + Asp + Glu + Lys + Arg +
                             (1|Cluster) + 0,
                           weights = Weight.nb.5.7)
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
  formula = Fitness.nb ~ predict(fitness.nb.aa.lm, newdata = peptide.data, re.form = NA) + (1|Cluster),
  weights = Weight.nb.5.7
)
summary(fitness.pred.aa.fit.cluster.lm)
#pred.aa.summary <- summary(fitness.pred.aa.fit.cluster.lm)
#pred.aa.summary$coefficients

intercept.only.lm <- lmer(
  data = peptide.data,
  formula = Fitness.nb ~ (1|Cluster),
  weights = Weight.nb.5.7
)
anova(fitness.nb.aa.lm, intercept.only.lm, test = "lrt")

# Combining the data by cluster for plotting.
by_cluster <-
  peptide.data %>% 
  group_by(Cluster) %>%
  summarise(Weight.nb.sum = sum(Weight.nb.5.7), 
            ISD.iupred2 = wtd.mean(ISD.iupred2, weights = Weight.nb.5.7), Fitness.nb = wtd.mean(Fitness.nb, weights = Weight.nb.5.7),
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
            Clustering.Six = wtd.mean(Clustering.Six, weights = Weight.nb.5.7),
            TangoAAsInAPRs = wtd.mean(TangoAAsInAPRs, weights = Weight.nb.5.7),
            CamSol.avg = wtd.mean(CamSol.avg, weights = Weight.nb.5.7),
            charge.pos = wtd.mean(charge.pos, weights = Weight.nb.5.7), charge.neg = wtd.mean(charge.neg, weights = Weight.nb.5.7),
            net.charge = wtd.mean(net.charge, weights = Weight.nb.5.7), pI = wtd.mean(pI, weights = Weight.nb.5.7),
            fit.full = wtd.mean(fit.full, weights = Weight.nb.5.7),
            fit.aa = wtd.mean(fit.aa, weights = Weight.nb.5.7)
  )
# by_cluster$Fitness.nb.weighted <- rep(NA, length(by_cluster$Fitness.nb))
# for (i in 1:length(by_cluster$Fitness.nb)) {
#   by_cluster[by_cluster$Cluster == i, "Fitness.nb.weighted"] <- 
#     weighted.mean(peptide.data[peptide.data$Cluster == i, "Fitness.nb"],
#                   peptide.data[peptide.data$Cluster == i, "Weight.nb"])
# }
# by_cluster$fit.full.weighted <- rep(NA, length(by_cluster$CamSol.avg))
# by_cluster$fit.aa.weighted <- rep(NA, length(by_cluster$TangoAAsInAPRs))
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
                  formula = Fitness.nb ~ fit.full,
                  weights = Weight.nb.sum)
summary(fit.full.lm)
pred.full.summary <- summary(fit.full.lm)
pred.full.summary$coefficients

fit.aa.lm <- lm(data = by_cluster,
                formula = Fitness.nb ~ fit.aa,
                weights = Weight.nb.sum)
summary(fit.aa.lm)
pred.aa.summary <- summary(fit.aa.lm)
pred.aa.summary$coefficients

with(by_cluster, weightedCorr(Fitness.nb, fit.full, method = "Pearson", weights = Weight.nb.sum))
with(by_cluster, weightedCorr(Fitness.nb, fit.aa, method = "Pearson", weights = Weight.nb.sum))

# Trying to flip the regression to X~Y, seeing what happens.
flipped.fit.full.lm <- lm(data = by_cluster,
                  formula = fit.full ~ Fitness.nb,
                  weights = Weight.nb.sum)
summary(flipped.fit.full.lm)
flipped.pred.full.summary <- summary(flipped.fit.full.lm)
flipped.pred.full.summary$coefficients
flipped.full.int <- -flipped.pred.full.summary$coefficients[1,1] / flipped.pred.full.summary$coefficients[2,1]
flipped.full.beta <- 1/flipped.pred.full.summary$coefficients[2,1]

flipped.fit.aa.lm <- lm(data = by_cluster,
                formula = fit.aa ~ Fitness.nb,
                weights = Weight.nb.sum)
summary(flipped.fit.aa.lm)
flipped.pred.aa.summary <- summary(flipped.fit.aa.lm)
flipped.pred.aa.summary$coefficients
flipped.aa.int <- -flipped.pred.aa.summary$coefficients[1,1] / flipped.pred.aa.summary$coefficients[2,1]
flipped.aa.beta <- 1/flipped.pred.aa.summary$coefficients[2,1]

# Plotting part A.
todays.date <- "5-24-20"
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
png(filename = paste("Scripts/Figures/fitness_pred_full_", todays.date, ".png", sep = ""),
    height = 490, width = 490)
ggplot(data = by_cluster,
       aes(x = fit.full,
           y = Fitness.nb,
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
  scale_y_continuous(limits = c(0, 2)) +
  #scale_y_continuous(breaks = sqrt(c(0.2, 0.5, 1, 2, 5)),
  #                   labels = c(0.2, 0.5, 1, 2, 5)) +
  #scale_x_continuous(breaks = sqrt(c(0.2, 0.5, 1)),
  #                   labels = c(0.2, 0.5, 1)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

# Plotting part B.
png(filename = paste("Scripts/Figures/fitness_pred_aacomp_", todays.date, ".png", sep = ""),
    height = 490, width = 490)
ggplot(data = by_cluster,
       aes(y = Fitness.nb,
           x = fit.aa,
           size = Weight.nb.sum,
           weight = Weight.nb.sum)
) +
  geom_point(alpha = 0.4) +
  #geom_abline(slope = 1, intercept = 0, color = cbbPalette[2], size = 1.5) +
  #stat_function(fun = function(x)pred.aa.summary$coefficients[1,1] + pred.aa.summary$coefficients[2,1]*x,
  #              geom = "line", color = cbbPalette[6], size = 1.5) +
  geom_smooth(method = "lm", color = cbbPalette[6], size = 1.5, se = F) +
  #geom_smooth(method = "loess", color = cbbPalette[2], size = 1.5, se = F) +
  ylab("Fitness") +
  xlab("AA-predicted fitness") +
  scale_y_continuous(limits = c(0, 2)) +
  #scale_y_continuous(breaks = sqrt(c(0.2, 0.5, 1, 2, 5)),
  #                   labels = c(0.2, 0.5, 1, 2, 5)) +
  #scale_x_continuous(breaks = sqrt(c(0.2, 0.5, 1)),
  #                   labels = c(0.2, 0.5, 1)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

# Zoomed out version.
png(filename = paste("Scripts/Figures/fitness_pred_aacomp_zoomout_", todays.date, ".png", sep = ""),
    height = 490, width = 490)
ggplot(data = by_cluster,
       aes(y = Fitness.nb,
           x = fit.aa,
           size = Weight.nb.sum,
           weight = Weight.nb.sum)
) +
  geom_point(alpha = 0.4) +
  #geom_abline(slope = 1, intercept = 0, color = cbbPalette[2], size = 1.5) +
  #stat_function(fun = function(x)pred.aa.summary$coefficients[1,1] + pred.aa.summary$coefficients[2,1]*x,
  #              geom = "line", color = cbbPalette[6], size = 1.5) +
  geom_smooth(method = "lm", color = cbbPalette[6], size = 1.5, se = F) +
  #geom_smooth(method = "loess", color = cbbPalette[2], size = 1.5, se = F) +
  ylab("Fitness") +
  xlab("AA-predicted fitness") +
  #scale_y_continuous(breaks = sqrt(c(0.2, 0.5, 1, 2, 5)),
  #                   labels = c(0.2, 0.5, 1, 2, 5)) +
  #scale_x_continuous(breaks = sqrt(c(0.2, 0.5, 1)),
  #                   labels = c(0.2, 0.5, 1)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()
