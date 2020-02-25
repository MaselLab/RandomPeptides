#Figure 1 script.

# Packages.
library(lme4)
library(stringr)
library(Hmisc)

# Functions
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Load peptide data.
peptide.data <- read.table(file = "Data/peptide_data_clusters_2-14-20.tsv", header = T, stringsAsFactors = F)

# Building the models of fitness ~ ISD, Clustering, Waltz, and CamSol.
isd.lme <-
  lmer(
    data = peptide.data,
    formula = log(Fitness.nb) ~ sqrt(ISD) + (1|Cluster),
    weights = Weight.nb
  )
drop1(isd.lme, test = "Chisq")
isd.summary <- summary(isd.lme)

clustering.lme <-
  lmer(
    data = peptide.data,
    formula = log(Fitness.nb) ~ Clustering.Six + (1|Cluster),
    weights = Weight.nb
  )
drop1(clustering.lme, test = "Chisq")
clustering.summary <- summary(clustering.lme)

camsol.lme <-
  lmer(
    data = peptide.data,
    formula = log(Fitness.nb) ~ CamSol.avg + (1|Cluster),
    weights = Weight.nb
  )
drop1(camsol.lme, test = "Chisq")
camsol.summary <- summary(camsol.lme)

waltz.lme <-
  lmer(
    data = peptide.data,
    formula = log(Fitness.nb) ~ WaltzBinary + (1|Cluster),
    weights = Weight.nb
  )
drop1(waltz.lme, test = "Chisq")

charge.lme <-
  lmer(
    data = peptide.data,
    formula = log(Fitness.nb) ~ charge.pos + charge.neg + (1|Cluster),
    weights = Weight.nb
  )
drop1(charge.lme, test = "Chisq")
charge.summary <- summary(charge.lme)

# Combining the data by cluster for plotting.
by_cluster <-
  peptide.data %>% 
  group_by(Cluster) %>%
  summarise(Weight.nb.sum = sum(Weight.nb), ISD = wtd.mean(ISD, weights = Weight.nb),
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
            net.charge = wtd.mean(net.charge, weights = Weight.nb)
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
# by_cluster$Waltz.fit <- rep(NA, length(by_cluster$WaltzBinary))
# by_cluster$ISD.weighted <- rep(NA, length(by_cluster$ISD))
# by_cluster$Clustering.Six.weighted <- rep(NA, length(by_cluster$Clustering.Six))
# by_cluster$CamSol.weighted <- rep(NA, length(by_cluster$CamSol.avg))
# by_cluster$WaltzBinary.mode <- rep(NA, length(by_cluster$WaltzBinary))
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
#   by_cluster[by_cluster$Cluster == i, "Waltz.fit"] <- 
#     weighted.mean(predict(waltz.lme,
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
#   by_cluster[by_cluster$Cluster == i, "WaltzBinary.mode"] <- 
#     getmode(peptide.data[peptide.data$Cluster == i, "WaltzBinary"])
# }

# Making the plots and exporting.
todays.date <- "2-24-20"
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
png(filename = paste("Scripts/Figures/fitness_isd_", todays.date, ".png", sep = ""),
    height = 500, width = 500)
ggplot(
  data = by_cluster,
  aes(
    y = log(Fitness.nb),
    x = sqrt(ISD),
    size = Weight.nb.sum,
    weight = Weight.nb.sum
  )
) +
  geom_point(alpha = 0.4) +
  stat_function(fun = function(x)isd.summary$coefficients[1,1]+isd.summary$coefficients[2,1]*x,
                geom = "line", color = cbbPalette[6], size = 1.5) +
  ylab("Fitness") +
  xlab(expression(sqrt("ISD"))) +
  scale_y_continuous(breaks = log(c(0.05, 0.5, 1, 2)),
                     labels = c(0.05, 0.5, 1, 2),
                     limits = log(c(0.04, 6))) +
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
    y = log(Fitness.nb),
    x = Clustering.Six,
    size = Weight.nb.sum,
    weight = Weight.nb.sum
  )
) +
  geom_point(alpha = 0.4) +
  stat_function(fun = function(x)clustering.summary$coefficients[1,1]+clustering.summary$coefficients[2,1]*x,
                geom = "line", color = cbbPalette[6], size = 1.5) +
  ylab("Fitness") +
  xlab("Clustering") +
  scale_y_continuous(breaks = log(c(0.05, 0.5, 1, 2)),
                     labels = c(0.05, 0.5, 1, 2),
                     limits = log(c(0.04, 6))) +
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
    y = log(Fitness.nb),
    x = CamSol.avg,
    size = Weight.nb.sum,
    weight = Weight.nb.sum
  )
) +
  geom_point(alpha = 0.4) +
  stat_function(fun = function(x)camsol.summary$coefficients[1,1]+camsol.summary$coefficients[2,1]*x,
                geom = "line", color = cbbPalette[6], size = 1.5) +
  ylab("Fitness") +
  xlab("CamSol") +
  scale_y_continuous(breaks = log(c(0.05, 0.5, 1, 2)),
                     labels = c(0.05, 0.5, 1, 2),
                     limits = log(c(0.04, 6))) +
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

png(filename = paste("Scripts/Figures/fitness_waltz_", todays.date, ".png", sep = ""),
    height = 500, width = 500)
ggplot(
  data = by_cluster,
  aes(
    y = log(Fitness.nb),
    x = round(WaltzBinary),
    group = factor(round(WaltzBinary)),
    weight = Weight.nb.sum
  )
) +
  stat_summary(fun.data = boxplot.quantiles, geom = "boxplot") +
  #stat_function(fun = function(x)-0.8099 - 0.1367*x, geom = "line", color = cbbPalette[6], size = 1.5) +
  ylab("Fitness") +
  xlab("Waltz predicted APRs") +
  scale_y_continuous(breaks = log(c(0.05, 0.5, 1, 2)),
                     labels = c(0.05, 0.5, 1, 2),
                     limits = log(c(0.04, 6))) +
  scale_x_continuous(breaks = c(0,1),
                     labels = c("None", "1+")) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

# Plot net charge against fitness.
png(filename = paste("Scripts/Figures/fitness_charge_", todays.date, ".png", sep = ""),
                     height = 500, width = 500)
ggplot(
  data = by_cluster,
  aes(
    y = log(Fitness.nb),
    x = net.charge,
    size = Weight.nb.sum,
    weight = Weight.nb.sum
  )
) +
  geom_point(alpha = 0.5) +
  #geom_smooth(method = "lm") +
  ylab("Fitness") +
  xlab("Expected net charge") +
  scale_y_continuous(breaks = log(c(0.05, 0.5, 1, 2)),
                     labels = c(0.05, 0.5, 1, 2),
                     limits = log(c(0.04, 6))) +
  stat_function(fun = function(x)charge.summary$coefficients[1,1] + charge.summary$coefficients[2,1]*x,
                geom = "line", color = cbbPalette[2], size = 1.5, xlim = c(0,11)) +
  stat_function(fun = function(x)charge.summary$coefficients[1,1] + charge.summary$coefficients[3,1]*x,
                geom = "line", color = cbbPalette[6], size = 1.5, xlim = c(-6, 0)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()