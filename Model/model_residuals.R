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
peptide.data$fit.var <- 1 / peptide.data$Weight.nb.5.7
peptide.data$fit.var.log <- ((1 / peptide.data$Fitness.nb) ^ 2) * peptide.data$fit.var
peptide.data$weight.log <- 1 / peptide.data$fit.var.log
peptide.data[1:10, c("PeptideID", "Weight.nb.5.7", "Weight.nb.log", "weight.log")]

# Building the model.
peptide.mixed.nb.log.lm <-
  lmer(
    data = peptide.data,
    formula = Fitness.nb ~
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
      #WaltzBinary +
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
    ,
    weights = Weight.nb.5.7
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
  summarise(Weight.nb.sum = sum(Weight.nb.5.7), ISD.iupred2 = wtd.mean(ISD.iupred2, weights = Weight.nb.5.7),
            Fitness.nb = wtd.mean(Fitness.nb, weights = Weight.nb.5.7),
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
            WaltzBinary = wtd.mean(WaltzBinary, weights = Weight.nb.5.7),
            CamSol.avg = wtd.mean(CamSol.avg, weights = Weight.nb.5.7),
            ISD.delta = wtd.mean(ISD.delta, weights = Weight.nb.5.7),
            net.charge = wtd.mean(net.charge, weights = Weight.nb.5.7),
            charge.pos = wtd.mean(charge.pos, weights = Weight.nb.5.7), charge.neg = wtd.mean(charge.neg, weights = Weight.nb.5.7),
            full.predict = wtd.mean(full.predict, weights = Weight.nb.5.7)
  )

# Assessing linearity of dependent and independent variables.
ggplot(data = by_cluster,
       aes(x = Leu, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Val, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Phe, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Met, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Ile, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Gly, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Ala, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Ser, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Cys, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Pro, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Trp, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = His, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Lys, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Arg, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Glu, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Asp, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Gln, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Asn, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Tyr, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

ggplot(data = by_cluster,
       aes(x = Thr, y = Fitness.nb, size = Weight.nb.sum, weight = Weight.nb.sum)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "gam") +
  geom_smooth(method = "lm", color = "red") +
  theme(legend.position = "none")

# Transform choice via Box-Cox method.
boxcox(peptide.data$Fitness.nb ~ 1, lambda = seq(-0.5, 0.5, by = 0.01))
# Optimal values include zero, so a log transform is chosen.

# Looking at the data.
todays.date <- "5-7-2020"
png(filename = paste("Scripts/Figures/fitness_histogram_all_", todays.date, ".png", sep = ""),
    width = 500, height = 500)
ggplot(
  data = peptide.data,
  aes(
    x = log(Fitness.nb)
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
    x = log(Fitness.nb)
  )
) +
  geom_histogram() +
  theme_bw(base_size = 28)
dev.off()

ggplot(
  data = by_cluster,
  aes(
    x = (Fitness.nb)^-0.2
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
    formula = log(Fitness.nb) ~
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
      #WaltzBinary +
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
    weights = Weight.nb.log.sum
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
by_cluster$full.resid <- by_cluster$Fitness.nb - by_cluster$full.predict
# Calculating weighted residuals by multiplying raw residuals by the square root of the weight.
by_cluster$full.resid.weighted <- by_cluster$full.resid * sqrt(by_cluster$Weight.nb.log.sum)
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
  scale_x_continuous(breaks = log(c(0.2, 0.5, 1)),
                     labels = c(0.2, 0.5, 1)) +
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
    formula = Fitness.nb ~
      Leu + Pro + Met + Trp + Ala +
      Val + Phe + Ile + Gly + Ser +
      Thr + Cys + Asn + Gln + Tyr +
      His + Asp + Glu + Lys + Arg +
      + 0
    ,
    weights = Weight.nb.sum
  )
aa.lm.summary <- summary(aa.lm)
aa.lm.summary
plot(aa.lm)
hist(weighted.residuals(aa.lm))
qqnorm(weighted.residuals(aa.lm))
qqline(weighted.residuals(aa.lm))

# And now residuals vs fitted.
png(filename = paste("Scripts/Figures/aa_only_residuals_weighted_fitted_", todays.date, ".png", sep = ""),
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

png(filename = paste("Scripts/Figures/aa_only_residuals_weighted_qq_", todays.date, ".png", sep = ""),
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
