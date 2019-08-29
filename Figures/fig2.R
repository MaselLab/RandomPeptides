# Figure 2 Script.

# Packages.
library(lme4)
library(stringr)
library(tidyverse)

# Load peptide data.
peptide.data <- read.table(file = "Data/peptide_data_clusters_7-20-19.tsv", header = T, stringsAsFactors = F)

# Full model.
fitness.nb.full.lm <- lmer(data = peptide.data,
                           formula = log(Fitness.nb) ~
                             Leu + Pro + Met + Trp + Ala +
                             Val + Phe + Ile + Gly + Ser +
                             Thr + Cys + Asn + Gln + Tyr +
                             His + Asp + Glu + Lys + Arg +
                             Clustering.Six +
                             WaltzBinary +
                             (1|Cluster) + 0,
                           weights = Weight.nb)
summary(fitness.nb.full.lm)

# Predicted estimated fitness with predicted fitness from the full model.
fitness.pred.fit.cluster.lm <- lmer(
  data = peptide.data,
  formula = log(Fitness.nb) ~ predict(fitness.nb.full.lm, newdata = peptide.data, re.form = NA) + (1|Cluster),
  weights = Weight.nb
)
fitness.pred.fit.cluster.lm

# Amino acid composition only.
fitness.nb.aa.lm <- lmer(data = peptide.data,
                           formula = log(Fitness.nb) ~
                             Leu + Pro + Met + Trp + Ala +
                             Val + Phe + Ile + Gly + Ser +
                             Thr + Cys + Asn + Gln + Tyr +
                             His + Asp + Glu + Lys + Arg +
                             (1|Cluster) + 0,
                           weights = Weight.nb)
summary(fitness.nb.aa.lm)

# Predicted estimated fitness with predicted fitness from the aa comp only model.
fitness.pred.aa.fit.cluster.lm <- lmer(
  data = peptide.data,
  formula = log(Fitness.nb) ~ predict(fitness.nb.aa.lm, newdata = peptide.data, re.form = NA) + (1|Cluster),
  weights = Weight.nb
)
fitness.pred.aa.fit.cluster.lm

# Combining the data by cluster for plotting.
by_cluster <-
  peptide.data %>% 
  group_by(Cluster) %>%
  summarise(Weight.nb = sum(Weight.nb), ISD = mean(ISD),
            Fitness.nb = mean(Fitness.nb),
            Leu = mean(Leu), Phe = mean(Phe), Met = mean(Met), Val = mean(Val), Ile = mean(Ile),
            Lys = mean(Lys), His = mean(His), Arg = mean(Arg), Glu = mean(Glu), Asp = mean(Asp),
            Gln = mean(Gln), Asn = mean(Asn), Gly = mean(Gly), Ala = mean(Ala), Pro = mean(Pro),
            Ser = mean(Ser), Trp = mean(Trp), Tyr = mean(Tyr), Thr = mean(Thr), Cys = mean(Cys),
            Clustering.Six = mean(Clustering.Six), WaltzBinary = mean(WaltzBinary), WaltzAAsInAPRs = mean(WaltzAAsInAPRs),
            Waltz.delta = mean(Waltz.delta), TangoBinary = mean(TangoBinary), TangoAAsInAPRs = mean(TangoAAsInAPRs),
            Tango.delta = mean(Tango.delta), Cys.squared = mean(Cys.squared), CamSol.avg = mean(CamSol.avg),
            AnchorAvg = mean(AnchorAvg), ISD.delta = mean(ISD.delta), mean.run.norm = mean(mean.run.norm),
            max.run.length = mean(max.run.length), Weight.nb = sum(Weight.nb))
by_cluster$Fitness.nb.weighted <- rep(NA, length(by_cluster$Fitness.nb))
for (i in 1:length(by_cluster$Fitness.nb)) {
  by_cluster[by_cluster$Cluster == i, "Fitness.nb.weighted"] <- 
    weighted.mean(peptide.data[peptide.data$Cluster == i, "Fitness.nb"],
                  peptide.data[peptide.data$Cluster == i, "Weight.nb"])
}
by_cluster$fit.full.weighted <- rep(NA, length(by_cluster$CamSol.avg))
by_cluster$fit.aa.weighted <- rep(NA, length(by_cluster$WaltzBinary))
for (i in 1:length(by_cluster$Fitness.nb)) {
  by_cluster[by_cluster$Cluster == i, "fit.full.weighted"] <- 
    weighted.mean(predict(fitness.nb.full.lm,
                          newdata = peptide.data[peptide.data$Cluster == i, ],
                          re.form = NA),
                  peptide.data[peptide.data$Cluster == i, "Weight.nb"])
  by_cluster[by_cluster$Cluster == i, "fit.aa.weighted"] <- 
    weighted.mean(predict(fitness.nb.aa.lm,
                          newdata = peptide.data[peptide.data$Cluster == i, ],
                          re.form = NA),
                  peptide.data[peptide.data$Cluster == i, "Weight.nb"])
}

# Plotting part A.
todays.date <- "8-26-19"
png(filename = paste("Scripts/Figures/fitness_pred_full_", todays.date, ".png", sep = ""),
    height = 500, width = 500)
ggplot(data = by_cluster,
       aes(x = fit.full.weighted,
           y = log(Fitness.nb.weighted),
           size = Weight.nb,
           weight = Weight.nb)
) +
  geom_point(alpha = 0.4) +
  stat_function(fun = function(x)-0.0008267 + 0.9991354*x, geom = "line", color = "blue", size = 1.5) +
  ylab("Fitness") +
  xlab("Predicted fitness") +
  scale_y_continuous(breaks = log(c(0.05, 0.5, 1, 2, 10)),
                     labels = c(0.05, 0.5, 1, 2, 10)) +
  scale_x_continuous(breaks = log(c(0.2, 0.5, 1)),
                     labels = c(0.2, 0.5, 1)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

# Plotting part B.
png(filename = paste("Scripts/Figures/fitness_pred_aacomp_", todays.date, ".png", sep = ""),
    height = 500, width = 500)
ggplot(data = by_cluster,
       aes(x = fit.aa.weighted,
           y = log(Fitness.nb.weighted),
           size = Weight.nb,
           weight = Weight.nb)
) +
  geom_point(alpha = 0.4) +
  stat_function(fun = function(x)-0.000727 + 0.999252*x, geom = "line", color = "blue", size = 1.5) +
  ylab("Fitness") +
  xlab("Predicted fitness") +
  scale_y_continuous(breaks = log(c(0.05, 0.5, 1, 2, 10)),
                     labels = c(0.05, 0.5, 1, 2, 10)) +
  scale_x_continuous(breaks = log(c(0.2, 0.5, 1)),
                     labels = c(0.2, 0.5, 1)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()
