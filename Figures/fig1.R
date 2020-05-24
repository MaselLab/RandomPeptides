#Figure 1 script.

# Packages.
library(tidyverse)
library(lme4)
library(stringr)
library(Hmisc)

# Functions
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Today's date
todays.date <- "5-23-2020"

# Load peptide data.
peptide.data <- read.table(file = "Scripts/RandomPeptides/Data/supplemental_table_1.tsv", header = T, stringsAsFactors = F)

# Building the models of fitness ~ ISD, Clustering, Tango, and CamSol.
isd.lme <-
  lmer(
    data = peptide.data,
    formula = Fitness.nb ~ sqrt(ISD.iupred2) + (1|Cluster),
    weights = Weight.nb.5.7
  )
drop1(isd.lme, test = "Chisq")
isd.summary <- summary(isd.lme)

clustering.lme <-
  lmer(
    data = peptide.data,
    formula = Fitness.nb ~ Clustering.Six + (1|Cluster),
    weights = Weight.nb.5.7
  )
drop1(clustering.lme, test = "Chisq")
clustering.summary <- summary(clustering.lme)

camsol.lme <-
  lmer(
    data = peptide.data,
    formula = Fitness.nb ~ CamSol.avg + (1|Cluster),
    weights = Weight.nb.5.7
  )
drop1(camsol.lme, test = "Chisq")
camsol.summary <- summary(camsol.lme)

tango.lme <-
  lmer(
    data = peptide.data,
    formula = Fitness.nb ~ TangoAAsInAPRs + (1|Cluster),
    weights = Weight.nb.5.7
  )
drop1(tango.lme, test = "Chisq")
tango.summary <- summary(tango.lme)

charge.lme <-
  lmer(
    data = peptide.data,
    formula = Fitness.nb ~ net.charge + (1|Cluster),
    weights = Weight.nb.5.7
  )
drop1(charge.lme, test = "Chisq")
charge.summary <- summary(charge.lme)

# Checking a two-parameter net charge model.
charge.2param.lme <-
  lmer(
    data = peptide.data,
    formula = Fitness.nb ~ abs(charge.neg) + abs(charge.pos) + (1|Cluster),
    weights = Weight.nb.5.7
  )
drop1(charge.2param.lme, test = "Chisq")
charge.2param.summary <- summary(charge.lme)
AIC(charge.lme)
AIC(charge.2param.lme)
# Single parameter model has lower AIC, so we go with the simpler model.
# But neither are significant.

# Calculating R-squared from predicted fitness vs estimated fitness from each model.
peptide.data$ISD.fit <- predict(isd.lme,
                                newdata = peptide.data,
                                type = "response",
                                random.only = F,
                                re.form = NA)
peptide.data$Clustering.fit <- predict(clustering.lme,
                                       newdata = peptide.data,
                                       type = "response",
                                       random.only = F,
                                       re.form = NA)
peptide.data$camsol.fit <- predict(camsol.lme,
                                   newdata = peptide.data,
                                   type = "response",
                                   random.only = F,
                                   re.form = NA)
peptide.data$Tango.fit <- predict(tango.lme,
                                  newdata = peptide.data,
                                  type = "response",
                                  random.only = F,
                                  re.form = NA)
peptide.data$charge.fit <- predict(charge.lme,
                                   newdata = peptide.data,
                                   type = "response",
                                   random.only = F,
                                   re.form = NA)

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
            charge.pos = wtd.mean(charge.pos, weights = Weight.nb.5.7), charge.neg = wtd.mean(charge.neg),
            net.charge = wtd.mean(net.charge, weights = Weight.nb.5.7),
            ISD.fit = wtd.mean(ISD.fit, weights = Weight.nb.5.7), Clustering.fit = wtd.mean(Clustering.fit, weights = Weight.nb.5.7),
            camsol.fit = wtd.mean(camsol.fit, weights = Weight.nb.5.7), Tango.fit = wtd.mean(Tango.fit, weights = Weight.nb.5.7),
            charge.fit = wtd.mean(charge.fit, weights = Weight.nb.5.7)
  )
# by_cluster$Fitness.nb.weighted <- rep(NA, length(by_cluster$Fitness.nb))
# for (i in 1:length(by_cluster$Fitness.nb)) {
#   by_cluster[by_cluster$Cluster == i, "Fitness.nb.weighted"] <- 
#     weighted.mean(peptide.data[peptide.data$Cluster == i, "Fitness.nb"],
#                   peptide.data[peptide.data$Cluster == i, "Weight.nb"])
# }
# by_cluster$ISD.fit <- rep(NA, length(by_cluster$Cluster))
# by_cluster$Clustering.fit <- rep(NA, length(by_cluster$Clustering.Six))
# by_cluster$CamSol.fit <- rep(NA, length(by_cluster$CamSol.avg))
# by_cluster$Tango.fit <- rep(NA, length(by_cluster$TangoAAsInAPRs))
# by_cluster$ISD.weighted <- rep(NA, length(by_cluster$ISD))
# by_cluster$Clustering.Six.weighted <- rep(NA, length(by_cluster$Clustering.Six))
# by_cluster$CamSol.weighted <- rep(NA, length(by_cluster$CamSol.avg))
# by_cluster$TangoAAsInAPRs.mode <- rep(NA, length(by_cluster$TangoAAsInAPRs))
# for (i in 1:length(by_cluster$Fitness.nb)) {
#   by_cluster[by_cluster$Cluster == i, "ISD.fit"] <- 
#     weighted.mean(predict(isd.lme,
#                           newdata = peptide.data[peptide.data$Cluster == i, ],
#                           re.form = NA),
#                   peptide.data[peptide.data$Cluster == i, "Weight.nb"])
#   by_cluster[by_cluster$Cluster == i, "Clustering.fit"] <- 
#     weighted.mean(predict(clustering.lme,
#                           newdata = peptide.data[peptide.data$Cluster == i, ],
#                           re.form = NA),
#                   peptide.data[peptide.data$Cluster == i, "Weight.nb"])
#   by_cluster[by_cluster$Cluster == i, "CamSol.fit"] <- 
#     weighted.mean(predict(camsol.lme,
#                           newdata = peptide.data[peptide.data$Cluster == i, ],
#                           re.form = NA),
#                   peptide.data[peptide.data$Cluster == i, "Weight.nb"])
#   by_cluster[by_cluster$Cluster == i, "Tango.fit"] <- 
#     weighted.mean(predict(tango.lme,
#                           newdata = peptide.data[peptide.data$Cluster == i, ],
#                           re.form = NA),
#                   peptide.data[peptide.data$Cluster == i, "Weight.nb"])
#   by_cluster[by_cluster$Cluster == i, "ISD.weighted"] <- 
#     weighted.mean(peptide.data[peptide.data$Cluster == i, "ISD"],
#                   peptide.data[peptide.data$Cluster == i, "Weight.nb"])
#   by_cluster[by_cluster$Cluster == i, "Clustering.Six.weighted"] <- 
#     weighted.mean(peptide.data[peptide.data$Cluster == i, "Clustering.Six"],
#                   peptide.data[peptide.data$Cluster == i, "Weight.nb"])
#   by_cluster[by_cluster$Cluster == i, "CamSol.weighted"] <- 
#     weighted.mean(peptide.data[peptide.data$Cluster == i, "CamSol.avg"],
#                   peptide.data[peptide.data$Cluster == i, "Weight.nb"])
#   by_cluster[by_cluster$Cluster == i, "TangoAAsInAPRs.mode"] <- 
#     getmode(peptide.data[peptide.data$Cluster == i, "TangoAAsInAPRs"])
# }

# And now for the R-squared.
isd.fit.lm <- lm(
  data = by_cluster,
  formula = Fitness.nb ~ ISD.fit,
  weights = Weight.nb.sum
)
summary(isd.fit.lm)
clustering.fit.lm <- lm(
  data = by_cluster,
  formula = Fitness.nb ~ Clustering.fit,
  weights = Weight.nb.sum
)
summary(clustering.fit.lm)
camsol.fit.lm <- lm(
  data = by_cluster,
  formula = Fitness.nb ~ camsol.fit,
  weights = Weight.nb.sum
)
summary(camsol.fit.lm)
Tango.fit.lm <- lm(
  data = by_cluster,
  formula = Fitness.nb ~ Tango.fit,
  weights = Weight.nb.sum
)
summary(Tango.fit.lm)
charge.fit.lm <- lm(
  data = by_cluster,
  formula = Fitness.nb ~ charge.fit,
  weights = Weight.nb.sum
)
summary(charge.fit.lm)

# Making the plots and exporting.
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
png(filename = paste("Scripts/Figures/fitness_isd_", todays.date, ".png", sep = ""),
    height = 500, width = 500)
ggplot(
  data = by_cluster,
  aes(
    y = Fitness.nb,
    x = sqrt(ISD.iupred2),
    size = Weight.nb.sum,
    weight = Weight.nb.sum
  )
) +
  geom_point(alpha = 0.4) +
  stat_function(fun = function(x)isd.summary$coefficients[1,1]+isd.summary$coefficients[2,1]*x,
                geom = "line", color = cbbPalette[6], size = 1.5) +
  #geom_smooth(method = "loess") +
  ylab("Fitness") +
  #geom_abline(slope = 1, intercept = 0, color = cbbPalette[2], size = 1.5) +
  xlab("ISD") +
  scale_y_continuous(limits = c(0, 2)) +
  #scale_y_continuous(breaks = log(c(0.05, 0.5, 1, 2)),
  #                   labels = c(0.05, 0.5, 1, 2),
  #                   limits = log(c(0.04, 6))) +
  scale_x_continuous(breaks = sqrt(c(0.04, 0.16, 0.36, 0.64)),
                     labels = c(0.04, 0.16, 0.36, 0.64)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()    

png(filename = paste("Scripts/Figures/fitness_clustering_", todays.date, ".png", sep = ""),
    height = 500, width = 500)
ggplot(
  data = by_cluster,
  aes(
    y = Fitness.nb,
    x = Clustering.Six,
    size = Weight.nb.sum,
    weight = Weight.nb.sum
  )
) +
  geom_point(alpha = 0.4) +
  #geom_abline(slope = 1, intercept = 0, color = cbbPalette[2], size = 1.5) +
  stat_function(fun = function(x)clustering.summary$coefficients[1,1]+clustering.summary$coefficients[2,1]*x,
                geom = "line", color = cbbPalette[6], size = 1.5) +
  geom_smooth(method = "loess") +
  ylab("Fitness") +
  xlab("Clustering") +
  scale_y_continuous(limits = c(0, 2)) +
  #scale_y_continuous(breaks = log(c(0.05, 0.5, 1, 2)),
  #                   labels = c(0.05, 0.5, 1, 2),
  #                   limits = log(c(0.04, 6))) +
  #scale_x_continuous(breaks = log(c(0.2, 0.5, 1)),
  #                   labels = c(0.2, 0.5, 1)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

png(filename = paste("Scripts/Figures/fitness_camsol_", todays.date, ".png", sep = ""),
    height = 500, width = 500)
ggplot(
  data = by_cluster,
  aes(
    y = Fitness.nb,
    x = CamSol.avg,
    size = Weight.nb.sum,
    weight = Weight.nb.sum
  )
) +
  geom_point(alpha = 0.4) +
  #geom_abline(slope = 1, intercept = 0, color = cbbPalette[2], size = 1.5) +
  stat_function(fun = function(x)camsol.summary$coefficients[1,1]+camsol.summary$coefficients[2,1]*x,
                geom = "line", color = cbbPalette[6], size = 1.5) +
  #geom_smooth(method = "loess") +
  ylab("Fitness") +
  xlab("CamSol") +
  scale_y_continuous(limits = c(0, 2)) +
  #scale_y_continuous(breaks = log(c(0.05, 0.5, 1, 2)),
  #                   labels = c(0.05, 0.5, 1, 2),
  #                   limits = log(c(0.04, 6))) +
  #scale_x_continuous(breaks = log(c(0.2, 0.5, 1)),
  #                   labels = c(0.2, 0.5, 1)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

# Getting the boxplot quantiles for graphing whiskers with 9% and 91% quantiles.
boxplot.quantiles <- function(x){
  qntls <- quantile(x, probs = c(0.09, 0.25, 0.5, 0.75, 0.91))
  names(qntls) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(qntls)
}

png(filename = paste("Scripts/Figures/fitness_tango_", todays.date, ".png", sep = ""),
    height = 500, width = 500)
ggplot(
  data = by_cluster,
  aes(
    y = Fitness.nb,
    x = round(TangoAAsInAPRs),
    size = Weight.nb.sum,
    weight = Weight.nb.sum
  )
) +
  geom_point(alpha = 0.4) +
  #stat_summary(fun.data = boxplot.quantiles, geom = "boxplot") +
  stat_function(fun = function(x)tango.summary$coefficients[1,1]+tango.summary$coefficients[2,1]*x,
                geom = "line", color = cbbPalette[6], size = 1.5) +
  #geom_smooth(method = "loess") +
  #geom_smooth(method = "lm") +
  ylab("Fitness") +
  xlab("AAs in Tango predicted APRs") +
  scale_y_continuous(limits = c(0, 2)) +
  #scale_y_continuous(breaks = log(c(0.05, 0.5, 1, 2)),
  #                   labels = c(0.05, 0.5, 1, 2),
  #                   limits = log(c(0.04, 6))) +
  #scale_x_continuous(breaks = c(0,1),
  #                   labels = c("None", "1+")) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

# Plot net charge against fitness.
png(filename = paste("Scripts/Figures/fitness_charge_", todays.date, ".png", sep = ""),
                     height = 500, width = 500)
ggplot(
  data = by_cluster,
  aes(
    y = Fitness.nb,
    x = net.charge,
    size = Weight.nb.sum,
    weight = Weight.nb.sum
  )
) +
  geom_point(alpha = 0.5) +
  #geom_abline(slope = 1, intercept = 0, color = cbbPalette[2], size = 1.5) +
  #geom_smooth(method = "lm") +
  #geom_smooth(method = "loess") +
  ylab("Fitness") +
  xlab("Expected net charge") +
  scale_y_continuous(limits = c(0, 2)) +
  #scale_y_continuous(breaks = log(c(0.05, 0.5, 1, 2)),
  #                   labels = c(0.05, 0.5, 1, 2),
  #                   limits = log(c(0.04, 6))) +
  stat_function(fun = function(x)charge.summary$coefficients[1,1] + charge.summary$coefficients[2,1]*x,
                geom = "line", color = cbbPalette[6], size = 1.5) +
  #stat_function(fun = function(x)charge.summary$coefficients[1,1] + charge.summary$coefficients[3,1]*x,
  #              geom = "line", color = cbbPalette[6], size = 1.5, xlim = c(-6, 0)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

# Zoomed out plots.
png(filename = paste("Scripts/Figures/fitness_tango_zoomout_", todays.date, ".png", sep = ""),
    height = 500, width = 500)
ggplot(
  data = by_cluster,
  aes(
    y = Fitness.nb,
    x = round(TangoAAsInAPRs),
    size = Weight.nb.sum,
    weight = Weight.nb.sum
  )
) +
  geom_point(alpha = 0.4) +
  #stat_summary(fun.data = boxplot.quantiles, geom = "boxplot") +
  stat_function(fun = function(x)tango.summary$coefficients[1,1]+tango.summary$coefficients[2,1]*x,
                geom = "line", color = cbbPalette[6], size = 1.5) +
  #geom_smooth(method = "loess") +
  #geom_smooth(method = "lm") +
  ylab("Fitness") +
  xlab("AAs in Tango predicted APRs") +
  #scale_y_continuous(limits = c(0, 2)) +
  #scale_y_continuous(breaks = log(c(0.05, 0.5, 1, 2)),
  #                   labels = c(0.05, 0.5, 1, 2),
  #                   limits = log(c(0.04, 6))) +
  #scale_x_continuous(breaks = c(0,1),
  #                   labels = c("None", "1+")) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

png(filename = paste("Scripts/Figures/fitness_camsol_zoomout_", todays.date, ".png", sep = ""),
    height = 500, width = 500)
ggplot(
  data = by_cluster,
  aes(
    y = Fitness.nb,
    x = CamSol.avg,
    size = Weight.nb.sum,
    weight = Weight.nb.sum
  )
) +
  geom_point(alpha = 0.4) +
  #geom_abline(slope = 1, intercept = 0, color = cbbPalette[2], size = 1.5) +
  stat_function(fun = function(x)camsol.summary$coefficients[1,1]+camsol.summary$coefficients[2,1]*x,
                geom = "line", color = cbbPalette[6], size = 1.5) +
  #geom_smooth(method = "loess") +
  ylab("Fitness") +
  xlab("CamSol") +
  #scale_y_continuous(breaks = log(c(0.05, 0.5, 1, 2)),
  #                   labels = c(0.05, 0.5, 1, 2),
  #                   limits = log(c(0.04, 6))) +
  #scale_x_continuous(breaks = log(c(0.2, 0.5, 1)),
  #                   labels = c(0.2, 0.5, 1)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

png(filename = paste("Scripts/Figures/fitness_isd_zoomout_", todays.date, ".png", sep = ""),
    height = 500, width = 500)
ggplot(
  data = by_cluster,
  aes(
    y = Fitness.nb,
    x = sqrt(ISD.iupred2),
    size = Weight.nb.sum,
    weight = Weight.nb.sum
  )
) +
  geom_point(alpha = 0.4) +
  stat_function(fun = function(x)isd.summary$coefficients[1,1]+isd.summary$coefficients[2,1]*x,
                geom = "line", color = cbbPalette[6], size = 1.5) +
  #geom_smooth(method = "loess") +
  ylab("Fitness") +
  #geom_abline(slope = 1, intercept = 0, color = cbbPalette[2], size = 1.5) +
  xlab("ISD") +
  #scale_y_continuous(breaks = log(c(0.05, 0.5, 1, 2)),
  #                   labels = c(0.05, 0.5, 1, 2),
  #                   limits = log(c(0.04, 6))) +
  scale_x_continuous(breaks = sqrt(c(0.04, 0.16, 0.36, 0.64)),
                     labels = c(0.04, 0.16, 0.36, 0.64)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()
