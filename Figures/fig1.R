#Figure 1 script.

# Packages.
library(lme4)
library(stringr)

# Functions
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Load peptide data.
peptide.data <- read.table(file = "Data/peptide_data_clusters_7-20-19.tsv", header = T, stringsAsFactors = F)

# Building the models of fitness ~ ISD, Clustering, Waltz, and CamSol.
isd.lme <-
  lmer(
    data = peptide.data,
    formula = log(Fitness.nb) ~ sqrt(ISD) + (1|Cluster),
    weights = Weight.nb
  )
drop1(isd.lme, test = "Chisq")

clustering.lme <-
  lmer(
    data = peptide.data,
    formula = log(Fitness.nb) ~ Clustering.Six + (1|Cluster),
    weights = Weight.nb
  )
drop1(clustering.lme, test = "Chisq")

camsol.lme <-
  lmer(
    data = peptide.data,
    formula = log(Fitness.nb) ~ CamSol.avg + (1|Cluster),
    weights = Weight.nb
  )
drop1(camsol.lme, test = "Chisq")

waltz.lme <-
  lmer(
    data = peptide.data,
    formula = log(Fitness.nb) ~ WaltzBinary + (1|Cluster),
    weights = Weight.nb
  )
drop1(waltz.lme, test = "Chisq")

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
by_cluster$ISD.fit <- rep(NA, length(by_cluster$Cluster))
by_cluster$Clustering.fit <- rep(NA, length(by_cluster$Clustering.Six))
by_cluster$CamSol.fit <- rep(NA, length(by_cluster$CamSol.avg))
by_cluster$Waltz.fit <- rep(NA, length(by_cluster$WaltzBinary))
by_cluster$ISD.weighted <- rep(NA, length(by_cluster$ISD))
by_cluster$Clustering.Six.weighted <- rep(NA, length(by_cluster$Clustering.Six))
by_cluster$CamSol.weighted <- rep(NA, length(by_cluster$CamSol.avg))
by_cluster$WaltzBinary.mode <- rep(NA, length(by_cluster$WaltzBinary))
for (i in 1:length(by_cluster$Fitness.nb)) {
  by_cluster[by_cluster$Cluster == i, "ISD.fit"] <- 
    weighted.mean(predict(isd.lme,
                          newdata = peptide.data[peptide.data$Cluster == i, ],
                          re.form = NA),
                  peptide.data[peptide.data$Cluster == i, "Weight.nb"])
  by_cluster[by_cluster$Cluster == i, "Clustering.fit"] <- 
    weighted.mean(predict(clustering.lme,
                          newdata = peptide.data[peptide.data$Cluster == i, ],
                          re.form = NA),
                  peptide.data[peptide.data$Cluster == i, "Weight.nb"])
  by_cluster[by_cluster$Cluster == i, "CamSol.fit"] <- 
    weighted.mean(predict(camsol.lme,
                          newdata = peptide.data[peptide.data$Cluster == i, ],
                          re.form = NA),
                  peptide.data[peptide.data$Cluster == i, "Weight.nb"])
  by_cluster[by_cluster$Cluster == i, "Waltz.fit"] <- 
    weighted.mean(predict(waltz.lme,
                          newdata = peptide.data[peptide.data$Cluster == i, ],
                          re.form = NA),
                  peptide.data[peptide.data$Cluster == i, "Weight.nb"])
  by_cluster[by_cluster$Cluster == i, "ISD.weighted"] <- 
    weighted.mean(peptide.data[peptide.data$Cluster == i, "ISD"],
                  peptide.data[peptide.data$Cluster == i, "Weight.nb"])
  by_cluster[by_cluster$Cluster == i, "Clustering.Six.weighted"] <- 
    weighted.mean(peptide.data[peptide.data$Cluster == i, "Clustering.Six"],
                  peptide.data[peptide.data$Cluster == i, "Weight.nb"])
  by_cluster[by_cluster$Cluster == i, "CamSol.weighted"] <- 
    weighted.mean(peptide.data[peptide.data$Cluster == i, "CamSol.avg"],
                  peptide.data[peptide.data$Cluster == i, "Weight.nb"])
  by_cluster[by_cluster$Cluster == i, "WaltzBinary.mode"] <- 
    getmode(peptide.data[peptide.data$Cluster == i, "WaltzBinary"])
}

# Making the plots and exporting.
todays.date <- "8-7-19"
png(filename = paste("Scripts/Figures/fitness_isd_", todays.date, ".png", sep = ""),
    height = 500, width = 500)
ggplot(
  data = by_cluster,
  aes(
    y = log(Fitness.nb.weighted),
    x = sqrt(ISD.weighted),
    size = Weight.nb,
    weight = Weight.nb
  )
) +
  geom_point(alpha = 0.4) +
  stat_function(fun = function(x)-1.342 + 1.093*x, geom = "line", color = "blue", size = 1.5) +
  ylab("Fitness") +
  xlab("ISD") +
  scale_y_continuous(breaks = log(c(0.05, 0.5, 1, 2, 10)),
                     labels = c(0.05, 0.5, 1, 2, 10)) +
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
    y = log(Fitness.nb.weighted),
    x = Clustering.Six.weighted,
    size = Weight.nb,
    weight = Weight.nb
  )
) +
  geom_point(alpha = 0.4) +
  stat_function(fun = function(x)-0.6385 - 0.2138*x, geom = "line", color = "blue", size = 1.5) +
  ylab("Fitness") +
  xlab("Clustering") +
  scale_y_continuous(breaks = log(c(0.05, 0.5, 1, 2, 10)),
                     labels = c(0.05, 0.5, 1, 2, 10)) +
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
    y = log(Fitness.nb.weighted),
    x = CamSol.weighted,
    size = Weight.nb,
    weight = Weight.nb
  )
) +
  geom_point(alpha = 0.4) +
  stat_function(fun = function(x)-1.0870 + 0.1697*x, geom = "line", color = "blue", size = 1.5) +
  ylab("Fitness") +
  xlab("CamSol") +
  scale_y_continuous(breaks = log(c(0.05, 0.5, 1, 2, 10)),
                     labels = c(0.05, 0.5, 1, 2, 10)) +
  #scale_x_continuous(breaks = log(c(0.2, 0.5, 1)),
  #                   labels = c(0.2, 0.5, 1)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

png(filename = paste("Scripts/Figures/fitness_waltz_", todays.date, ".png", sep = ""),
    height = 500, width = 500)
ggplot(
  data = by_cluster,
  aes(
    y = log(Fitness.nb.weighted),
    x = WaltzBinary.mode,
    size = Weight.nb,
    weight = Weight.nb
  )
) +
  geom_point(alpha = 0.4) +
  stat_function(fun = function(x)-0.8099 - 0.1367*x, geom = "line", color = "blue", size = 1.5) +
  ylab("Fitness") +
  xlab("Waltz predicted APRs") +
  scale_y_continuous(breaks = log(c(0.05, 0.5, 1, 2, 10)),
                     labels = c(0.05, 0.5, 1, 2, 10)) +
  scale_x_continuous(breaks = c(0,1),
                     labels = c("None", "1+")) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()
