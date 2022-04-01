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
todays.date <- Sys.Date()

# Load peptide data.
peptide.data <- read.table(file = "Scripts/RandomPeptides/Data/supplemental_dataset_1.tsv", header = T, stringsAsFactors = F)

# Calculating disorder propensity.
source("~/MaselLab/RandomPeptides/Scripts/RandomPeptides/Metrics/aa_comp_metrics.R")
peptide.data$disorder <- mean.metric.calculator(
  aa.sequence = peptide.data$AASeq,
  metric = "disorder"
)
hist(peptide.data$disorder)

# Combining the data by cluster using the max weight peptide from each cluster.
cluster_max <- peptide.data %>% group_by(Cluster) %>% filter(WEIGHT == max(WEIGHT))

# Building the models of fitness ~ ISD, Clustering, Tango, and CamSol.
isd.lm <-
  lm(
    data = cluster_max,
    formula = FITNESS ~ sqrt(ISD.iupred2) ,
    weights = WEIGHT
  )
drop1(isd.lm, test = "Chisq")
isd.summary <- summary(isd.lm)
isd.summary

clustering.lm <-
  lm(
    data = cluster_max,
    formula = FITNESS ~ Clustering.Six ,
    weights = WEIGHT
  )
drop1(clustering.lm, test = "Chisq")
clustering.summary <- summary(clustering.lm)
clustering.summary

camsol.lm <-
  lm(
    data = cluster_max,
    formula = FITNESS ~ CamSol.avg ,
    weights = WEIGHT
  )
drop1(camsol.lm, test = "Chisq")
camsol.summary <- summary(camsol.lm)
camsol.summary

tango.lm <-
  lm(
    data = cluster_max,
    formula = FITNESS ~ TangoAAsInAPRs ,
    weights = WEIGHT
  )
drop1(tango.lm, test = "Chisq")
tango.summary <- summary(tango.lm)
tango.summary

waltz.lm <-
  lm(
    data = cluster_max,
    formula = FITNESS ~ WaltzAAsInAPRs ,
    weights = WEIGHT
  )
drop1(waltz.lm, test = "Chisq")
waltz.summary <- summary(waltz.lm)
waltz.summary

disorder.lm <-
  lm(
    data = cluster_max,
    formula = FITNESS ~ disorder ,
    weights = WEIGHT
  )
drop1(disorder.lm, test = "Chisq")
disorder.summary <- summary(disorder.lm)
disorder.summary

FITNESS.aa.lm <- lm(data = cluster_max,
                         formula = FITNESS ~
                           Leu + Pro + Met + Trp + Ala +
                           Val + Phe + Ile + Gly + Ser +
                           Thr + Cys + Asn + Gln + Tyr +
                           His + Asp + Glu + Lys + Arg +
                           + 0,
                         weights = WEIGHT)
drop1(FITNESS.aa.lm, test = "Chisq")
aa.summary <- summary(FITNESS.aa.lm)
aa.summary

# Checking whether delta ISD improves the model.
delta.isd.lm <- 
  lm(
    data = cluster_max,
    formula = FITNESS ~ 
      Leu + Pro + Met + Trp + Ala +
      Val + Phe + Ile + Gly + Ser +
      Thr + Cys + Asn + Gln + Tyr +
      His + Asp + Glu + Lys + Arg +
      ISD.delta  + 0,
    weights = WEIGHT
  )
drop1(delta.isd.lm, test = "Chisq")

# Calculating predicted fitness each model.
cluster_max$ISD.fit <- predict(isd.lm,
                                newdata = cluster_max,
                                type = "response")
cluster_max$Clustering.fit <- predict(clustering.lm,
                                       newdata = cluster_max,
                                       type = "response")
cluster_max$CamSol.fit <- predict(camsol.lm,
                                   newdata = cluster_max,
                                   type = "response")
cluster_max$Tango.fit <- predict(tango.lm,
                                  newdata = cluster_max,
                                  type = "response")
cluster_max$disorder.fit <- predict(disorder.lm,
                                     newdata = cluster_max,
                                     type = "response")
cluster_max$fit.aa <- 
  predict(
    FITNESS.aa.lm,
    newdata = cluster_max,
    type = "response"
  )

# Checking significance for AA model.
intercept.only.lm <- lm(
  data = cluster_max,
  formula = FITNESS ~ 1,
  weights = WEIGHT
)
lrt_aa_int <- anova(FITNESS.aa.lm, intercept.only.lm, test = "LRT")
lrt_aa_int$`Pr(>Chi)`

# Looking at effect size.
quantile(cluster_max$FITNESS)
quantile(cluster_max$ISD.fit, probs = seq(0, 1, by = 0.1))
quantile(cluster_max$Clustering.fit, probs = seq(0, 1, by = 0.1))
quantile(cluster_max$CamSol.fit, probs = seq(0, 1, by = 0.1))
quantile(cluster_max$Tango.fit, probs = seq(0, 1, by = 0.1))
quantile(cluster_max$disorder.fit, probs = seq(0, 1, by = 0.1))
quantile(cluster_max$fit.aa, probs = seq(0, 1, by = 0.1))

# And now for the R-squared.
isd.summary
clustering.summary
camsol.summary
tango.summary
disorder.summary
fit.aa.lm <- lm(data = cluster_max,
                formula = FITNESS ~ Leu + Pro + Met + Trp + Ala +
                  Val + Phe + Ile + Gly + Ser +
                  Thr + Cys + Asn + Gln + Tyr +
                  His + Asp + Glu + Lys #+ Arg
                ,
                weights = WEIGHT)
summary(fit.aa.lm)

delta.isd.lm <- lm(
  data = cluster_max,
  formula = FITNESS ~ ISD.delta,
  weights = WEIGHT
)
summary(delta.isd.lm)

## Why is the fixed effect model for clustering significant but the mixed effect model is not?
## Maybe the within-cluster trend is contrary to the between-cluster trend.
math.sign <- function(x){
  if (is.na(x)){
    return(NA_character_)
  } else if (!is.numeric(x)){
    return("Not numeric")
  } else if (x > 0){
    return("+")
  } else if (x < 0){
    return("-")
  } else if (x == 0){
    return("0")
  } else {
    return("ERROR")
  }
}

peptide.data %>% group_by(Cluster) %>% summarise(
  Slope = lm(formula = FITNESS ~ Clustering.Six, weights = WEIGHT)$coefficients[2],
  Sign = math.sign(Slope)
) %>% group_by(Sign) %>% summarise(N = n())

# Plotting the data.
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cluster_max <- cluster_max %>% mutate(
  "IUPred (square root)" = sqrt(ISD.iupred2),
  "CamSol" = CamSol.avg,
  "#Aggregation-Prone AAs (Tango)" = TangoAAsInAPRs,
  "Disorder propensity" = disorder,
  "AA-predicted fitness" = fit.aa,
  "Hydrophobic clustering" = Clustering.Six
)

fit.list <- c(
  "IUPred (square root)",
  "CamSol",
  "#Aggregation-Prone AAs (Tango)",
  "Disorder propensity",
  "AA-predicted fitness",
  "Hydrophobic clustering"
)

fit.plotting <- function(.x, y.limit = c(0,2)){
  require(stringr)
  require(ggplot2)
  
  # Making a new DF for plotting.
  plot.df <- data.frame(
    "Fitness" = cluster_max[["FITNESS"]],
    "Fit" = cluster_max[[.x]],
    "WEIGHT" = cluster_max[["WEIGHT"]]
  )
  
  fit.plot <-
    ggplot(
      data = plot.df,
      aes(
        y = Fitness,
        x = Fit,
        size = WEIGHT,
        weight = WEIGHT
      )
    ) +
    ylab("Fitness") + xlab(.x) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "lm",
                se = F, size = 1.5, color = cbbPalette[6]) +
    scale_y_continuous(limits = y.limit) +
    #scale_x_continuous(limits = x.limit) +
    theme_bw(base_size = 28) +
    theme(legend.position = "none")
  
  return(fit.plot)
}

fit.plots <- lapply(fit.list, fit.plotting)

names(fit.plots) <- c("ISD", "CamSol", "Tango", "Disorder", "AA", "Clustering")

mapply(ggsave,
       file = paste0("Scripts/Figures/fitness_", names(fit.plots), "_", todays.date, ".png"),
       plot = fit.plots,
       height = 7, width = 7, units = "in", dpi = 300)


fit.zoomout.plots <- lapply(fit.list, fit.plotting, y.limit = c(0, 5.8))

names(fit.zoomout.plots) <- c("ISD", "CamSol", "Tango", "Disorder", "AA", "Clustering")

mapply(ggsave,
       file = paste0("Scripts/Figures/zoomout_fitness_", names(fit.zoomout.plots), "_", todays.date, ".png"),
       plot = fit.zoomout.plots,
       height = 7, width = 7, units = "in", dpi = 300)
