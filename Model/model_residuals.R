# QQ plot and residuals plot for our mixed model.

# Loading packages.
library(tidyverse)
library(lme4)
library(Hmisc)
library(MASS)
library(stats)

# Loading the data.
peptide.data <- read.table(file = "Scripts/RandomPeptides/Data/supplemental_table_1.tsv", header = T, stringsAsFactors = F)

# Messing with the weights.
peptide.data$fit.var <- 1 / peptide.data$WEIGHT
peptide.data$fit.var.log <- ((1 / peptide.data$FITNESS) ^ 2) * peptide.data$fit.var
peptide.data$weight.log <- 1 / peptide.data$fit.var.log
peptide.data$fit.var.sqrt <- (1 / (4 * peptide.data$FITNESS)) * peptide.data$fit.var
peptide.data[1:10, c("PeptideID", "WEIGHT", "Weight.nb.log", "weight.log")]
peptide.data[1:10, c("PeptideID", "fit.var", "fit.var.log", "fit.var.sqrt")]

# Building the model.
peptide.mixed.nb.log.lm <-
  lmer(
    data = peptide.data,
    formula = log(FITNESS) ~
      Leu + Pro + Met + Trp + Ala +
      Val + Phe + Ile + Gly + Ser +
      Thr + Cys + Asn + Gln + Tyr +
      His + Asp + Glu + Lys + Arg +
      #Cys.squared +
      #Clustering.Six +
      #CamSol.avg +
      #sqrt(ISD.iupred2) +
      #TangoAAsInAPRs +
      #WaltzAAsInAPRs +
      #WaltzNumAPRs +
      #AnchorAvg +
      #TangoBinary +
      #TangoAAsInAPRs +
      #Tango.delta +
      #Waltz.delta +
      #mean.run.norm +
      #max.run.length +
      #ISD.delta +
      #PredHel +
      #ExpAA +
      #net.charge +
      #I(abs(net.charge)) +
      #net.pos +
      #net.neg +
      #net.neither +
      #pI +
      #pI.ecoli +
      #AA.cost.ecoli +
      (1|Cluster)
    + 0
    # ,
    # weights = WEIGHT
  )
summary(peptide.mixed.nb.log.lm)
drop1(peptide.mixed.nb.log.lm, test = "Chisq")
peptide.mixed.log.lm.lrt <- drop1(peptide.mixed.nb.log.lm, test = "Chisq")
peptide.mixed.log.lm.lrt$`Pr(Chi)`
AIC(peptide.mixed.nb.log.lm)
peptide.mixed.nb.log.summary <- summary(peptide.mixed.nb.log.lm)
peptide.mixed.nb.log.summary

# Predicting fitness with predicted fitness.
peptide.data$full.predict <- predict(peptide.mixed.nb.log.lm,
                                     newdata = peptide.data,
                                     type = "response",
                                     random.only = F,
                                     re.form = NA)

# Combining data by cluster.
# Merging the data by cluster.
by_cluster <-
  peptide.data %>% 
  group_by(Cluster) %>%
  summarise(Weight.nb.sum = sum(WEIGHT), ISD.iupred2 = wtd.mean(ISD.iupred2, weights = WEIGHT),
            FITNESS = wtd.mean(FITNESS, weights = WEIGHT), Weight.nb.log.sum = sum(Weight.nb.log),
            FITNESS.log = wtd.mean(log(FITNESS), weights = Weight.nb.log),
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
            ISD.delta = wtd.mean(ISD.delta, weights = WEIGHT),
            net.charge = wtd.mean(net.charge, weights = WEIGHT),
            charge.pos = wtd.mean(charge.pos, weights = WEIGHT), charge.neg = wtd.mean(charge.neg, weights = WEIGHT),
            full.predict = wtd.mean(full.predict, weights = WEIGHT),
            Trp.unweighted = mean(Trp), Arg.unweighted = mean(Arg)
  )

max.weight.cluster <- peptide.data %>% group_by(Cluster) %>% filter(WEIGHT == max(WEIGHT))
max.logweight.cluster <- peptide.data %>% group_by(Cluster) %>% filter(Weight.nb.log == max(Weight.nb.log))

# Assessing linearity of dependent and independent variables.
ggplot(data = by_cluster,
       aes(x = Leu, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Val, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Phe, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Met, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Ile, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Gly, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Ala, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Ser, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Cys, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Pro, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

png(filename = "Scripts/Figures/trp_vs_fitness.png")
ggplot(data = by_cluster,
       aes(x = Trp, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "gam") +
  geom_smooth(method = "lm", color = "red") +
  ylab("Fitness") +
  xlab("Tryptophan residues") +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

png(filename = "Scripts/Figures/arg_vs_trp.png", height = 600, width = 600)
ggplot(data = by_cluster,
       aes(x = Trp.unweighted, y = Arg.unweighted)) +
  geom_jitter(width = 0.15, height = 0.15) +
  geom_smooth(method = "loess") +
  ylab("#arg") +
  xlab("#trp") +
  theme_bw(base_size = 28)
dev.off()

png(filename = "Scripts/Figures/trp_vs_arg.png", height = 600, width = 600)
ggplot(data = by_cluster,
       aes(x = Arg.unweighted, y = Trp.unweighted)) +
  geom_jitter(width = 0.15, height = 0.15) +
  geom_smooth(method = "loess") +
  xlab("#arg") +
  ylab("#trp") +
  theme_bw(base_size = 28)
dev.off()

#cor.test(by_cluster$Trp, by_cluster$Arg, method = "spearman")

png(filename = "Scripts/Figures/trp_histogram.png")
ggplot(data = by_cluster,
       aes(x = Trp)) +
  geom_histogram(bins = 7) +
  theme_bw(base_size = 28)
dev.off()

wtd.mean(by_cluster$Trp, weights = by_cluster$Weight.nb.sum)

ggplot(data = by_cluster,
       aes(x = His, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Lys, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Arg, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Glu, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Asp, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Gln, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Asn, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Tyr, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Thr, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "gam") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

# Fitness vs predicted fitness.
ggplot(data = by_cluster,
       aes(x = full.predict, y = FITNESS, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

# Looking at fitness vs variance.
ggplot(data = peptide.data,
       aes(y = sqrt(fit.var), x = FITNESS)) +
  ylab("Estimated standard error") +
  xlab("Fitness") +
  geom_point() +
  theme_bw(base_size = 28)

ggplot(data = peptide.data,
       aes(x = WEIGHT, y = log(FITNESS))) +
  geom_point()

ggplot(data = peptide.data,
       aes(x = fit.var.log, y = FITNESS)) +
  geom_point()

ggplot(data = peptide.data,
       aes(x = fit.var.sqrt, y = sqrt(FITNESS))) +
  geom_point()

# Looking at weight vs fitness plots
todays.date <- "5-24-2020"

png(filename = paste("Scripts/Figures/logfit_normweights_", todays.date, ".png", sep = ""))
ggplot(data = peptide.data,
       aes(y = WEIGHT, x = log(FITNESS))) +
  geom_point() +
  geom_smooth() +
  xlab("log(Fitness)") +
  ylab("Weight") +
  theme_bw(base_size = 28)
dev.off()

png(filename = paste("Scripts/Figures/logfit_logweights_", todays.date, ".png", sep = ""))
ggplot(data = max.logweight.cluster,
       aes(y = Weight.nb.log, x = log(FITNESS))) +
  geom_point() +
  geom_smooth() +
  xlab("log(Fitness)") +
  ylab("Weight") +
  theme_bw(base_size = 28)
dev.off()

png(filename = paste("Scripts/Figures/normfit_normweights_", todays.date, ".png", sep = ""))
ggplot(data = max.weight.cluster,
       aes(y = WEIGHT, x = FITNESS)) +
  geom_point() +
  geom_smooth() +
  xlab("Fitness") +
  ylab("Weight") +
  theme_bw(base_size = 28)
dev.off()

# Transform choice via Box-Cox method.
boxcox(peptide.data$FITNESS ~ 1, lambda = seq(-0.5, 0.5, by = 0.01))
# Optimal values include zero, so a log transform is chosen.

# Looking at the data.
png(filename = paste("Scripts/Figures/fitness_histogram_all_", todays.date, ".png", sep = ""),
    width = 500, height = 500)
ggplot(
  data = peptide.data,
  aes(
    x = FITNESS
  )
) +
  geom_histogram() +
  theme_bw(base_size = 28)
dev.off()

png(filename = paste("Scripts/Figures/fitness_histogram_cluster_", todays.date, ".png", sep = ""),
    width = 500, height = 500)
ggplot(
  data = by_cluster,
  aes(
    x = FITNESS
  )
) +
  geom_histogram() +
  theme_bw(base_size = 28)
dev.off()

ggplot(
  data = by_cluster,
  aes(
    x = (FITNESS)^-0.2
  )
) +
  geom_histogram() +
  theme_bw(base_size = 28)

# Diagnostics
plot(peptide.mixed.nb.log.lm)
ggplot(
  data = data.frame("residuals.model" = peptide.mixed.nb.log.summary$residuals),
  aes(
    x = residuals.model
  )
) +
  geom_histogram() +
  xlab("residuals") +
  theme_bw(base_size = 28)
qqmath(peptide.mixed.nb.log.lm)

# Regular linear model with weighted means from clusters.
full.lm <-
  lm(
    data = by_cluster,
    formula = FITNESS ~
      Leu + Pro + Met + Trp + Ala +
      Val + Phe + Ile + Gly + Ser +
      Thr + Cys + Asn + Gln + Tyr +
      His + Asp + Glu + Lys + Arg +
      #Cys.squared +
      #Clustering.Six +
      #CamSol.avg +
      #sqrt(ISD.iupred2) +
      #TangoAAsInAPRs +
      #WaltzAAsInAPRs +
      #WaltzNumAPRs +
      #AnchorAvg +
      #TangoBinary +
      #TangoAAsInAPRs +
      #Tango.delta +
      #Waltz.delta +
      #mean.run.norm +
      #max.run.length +
      #ISD.delta +
      #PredHel +
      #ExpAA +
      #net.charge +
      #I(abs(net.charge)) +
      #net.pos +
      #net.neg +
      #net.neither +
      #abs(charge.pos) +
      #abs(charge.neg) +
      #pI +
      #pI.ecoli +
      #AA.cost.ecoli
      + 0
    ,
    weights = Weight.nb.sum
  )
full.lm.summary <- summary(full.lm)
full.lm.summary
plot(full.lm)
hist(weighted.residuals(full.lm))
qqnorm(weighted.residuals(full.lm))
qqline(weighted.residuals(full.lm))
# There's a mutlicollinearity problem here where the algorithm is having trouble inverting
# the matrix because net charge and the charged amino acids are highly correlated.

# Calculating residuals using only the fixed effects from the mixed model.
by_cluster$full.resid <- by_cluster$FITNESS - by_cluster$full.predict
# Calculating weighted residuals by multiplying raw residuals by the square root of the weight.
by_cluster$full.resid.weighted <- by_cluster$full.resid * sqrt(by_cluster$Weight.nb.sum)
plot(by_cluster$full.predict, by_cluster$full.resid.weighted)
plot(by_cluster$full.predict, by_cluster$full.resid)

# Plotting residuals vs fitted.
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
png(filename = paste("Scripts/Figures/residuals_weighted_fitted_", todays.date, ".png", sep = ""),
    width = 500, height = 500)
ggplot(
  data = by_cluster,
  aes(
    x = full.predict,
    y = full.resid.weighted
  )
) +
  geom_point() +
  geom_smooth(color = cbbPalette[6]) +
  xlab("Fitted values") +
  ylab("Weighted residuals") +
  #scale_x_continuous(breaks = log(c(0.2, 0.5, 1)),
  #                   labels = c(0.2, 0.5, 1)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

# Those residuals are not around zero.

qqnorm(residuals(full.lm, type = "response"))
qqnorm(residuals(full.lm, type = "pearson"))
qqnorm(weighted.residuals(full.lm))
qqline(weighted.residuals(full.lm))

# Looking at the same plots of residuals, but for an amino acid only model.
# Regular linear model with weighted means from clusters.
aa.lm <-
  lm(
    data = by_cluster,
    formula = log(FITNESS) ~
      Leu + Pro + Met + Trp + Ala +
      Val + Phe + Ile + Gly + Ser +
      Thr + Cys + Asn + Gln + Tyr +
      His + Asp + Glu + Lys + Arg +
      + 0
    # ,
    # weights = Weight.nb.log.sum
  )
aa.lm.summary <- summary(aa.lm)
aa.lm.summary
plot(aa.lm)
hist(weighted.residuals(aa.lm))
qqnorm(weighted.residuals(aa.lm))
qqline(weighted.residuals(aa.lm))

# And now residuals vs fitted.
png(filename = paste("Scripts/Figures/aa_only_residuals_unweighted_logfitted_", todays.date, ".png", sep = ""),
    width = 500, height = 500)
ggplot(
  data = data.frame(
    "fitted.y" = predict(aa.lm, by_cluster),
    "residuals.model" = weighted.residuals(aa.lm)
  ),
  aes(
    x = fitted.y,
    y = residuals.model
  )
) +
  geom_point() +
  geom_smooth(color = cbbPalette[6], lwd = 1.2) +
  xlab("Fitted values") +
  ylab("Weighted residuals") +
  # scale_x_continuous(breaks = log(c(0.2, 0.5, 1, 1.5)),
  #                    labels = c(0.2, 0.5, 1, 1.5)) +
  #scale_x_continuous(breaks = log(c(0.75, 1, 1.5)),
  #                   labels = c(0.75, 1, 1.5)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

png(filename = paste("Scripts/Figures/aa_only_residuals_logfitted_qq_", todays.date, ".png", sep = ""),
    width = 500, height = 500)
ggplot(
  data = data.frame(
    "fitted.y" = predict(aa.lm, by_cluster),
    "residuals.model" = weighted.residuals(aa.lm)
  ),
  aes(
    sample = residuals.model
  )
) +
  geom_qq() +
  geom_qq_line(color = cbbPalette[6], lwd = 1.2) +
  xlab("Theoretical quantiles") +
  ylab("Weighted residuals' quantiles") +
  #scale_x_continuous(breaks = log(c(0.2, 0.5, 1, 1.5)),
  #                   labels = c(0.2, 0.5, 1, 1.5)) +
  theme_bw(base_size = 28)
dev.off()
