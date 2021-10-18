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
peptide.data <- read.table(file = "Scripts/RandomPeptides/Data/supplemental_table_1.tsv", header = T, stringsAsFactors = F)

# Calculating disorder propensity.
source("~/MaselLab/RandomPeptides/Scripts/RandomPeptides/Metrics/aa_comp_metrics.R")
peptide.data$disorder <- mean.metric.calculator(
  aa.sequence = peptide.data$AASeq,
  metric = "disorder"
)
hist(peptide.data$disorder)

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

disorder.lme <-
  lmer(
    data = peptide.data,
    formula = Fitness.nb ~ disorder + (1|Cluster),
    weights = Weight.nb.5.7
  )
drop1(disorder.lme, test = "Chisq")
disorder.summary <- summary(disorder.lme)

fitness.nb.aa.lm <- lmer(data = peptide.data,
                         formula = Fitness.nb ~
                           Leu + Pro + Met + Trp + Ala +
                           Val + Phe + Ile + Gly + Ser +
                           Thr + Cys + Asn + Gln + Tyr +
                           His + Asp + Glu + Lys + Arg +
                           (1|Cluster) + 0,
                         weights = Weight.nb.5.7)

# Calculating predicted fitness each model.
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
peptide.data$CamSol.fit <- predict(camsol.lme,
                                   newdata = peptide.data,
                                   type = "response",
                                   random.only = F,
                                   re.form = NA)
peptide.data$Tango.fit <- predict(tango.lme,
                                  newdata = peptide.data,
                                  type = "response",
                                  random.only = F,
                                  re.form = NA)
peptide.data$disorder.fit <- predict(disorder.lme,
                                     newdata = peptide.data,
                                     type = "response",
                                     random.only = F,
                                     re.form = NA)
peptide.data$fit.aa <- 
  predict(
    fitness.nb.aa.lm,
    newdata = peptide.data,
    type = "response",
    random.only = F,
    re.form = NA
  )

# Checking significance for AA model.
intercept.only.lm <- lmer(
  data = peptide.data,
  formula = Fitness.nb ~ (1|Cluster),
  weights = Weight.nb.5.7
)
anova(fitness.nb.aa.lm, intercept.only.lm, test = "lrt")

# Looking at effect size.
quantile(peptide.data$Fitness.nb)
quantile(peptide.data$ISD.fit, probs = seq(0, 1, by = 0.1))
quantile(peptide.data$Clustering.fit, probs = seq(0, 1, by = 0.1))
quantile(peptide.data$CamSol.fit, probs = seq(0, 1, by = 0.1))
quantile(peptide.data$Tango.fit, probs = seq(0, 1, by = 0.1))
quantile(peptide.data$disorder.fit, probs = seq(0, 1, by = 0.1))
quantile(peptide.data$fit.aa, probs = seq(0, 1, by = 0.1))

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
            disorder = wtd.mean(disorder, weights = Weight.nb.5.7),
            ISD.delta = wtd.mean(ISD.delta, weights = Weight.nb.5.7),
            ISD.fit = wtd.mean(ISD.fit, weights = Weight.nb.5.7), Clustering.fit = wtd.mean(Clustering.fit, weights = Weight.nb.5.7),
            CamSol.fit = wtd.mean(CamSol.fit, weights = Weight.nb.5.7), Tango.fit = wtd.mean(Tango.fit, weights = Weight.nb.5.7),
            disorder.fit = wtd.mean(disorder.fit, weights = Weight.nb.5.7),
            AA.fit = wtd.mean(fit.aa, weights = Weight.nb.5.7)
  )

# And now for the R-squared.
isd.fit.lm <- lm(
  data = by_cluster,
  formula = Fitness.nb ~ sqrt(ISD.iupred2),
  weights = Weight.nb.sum
)
summary(isd.fit.lm)
clustering.fit.lm <- lm(
  data = by_cluster,
  formula = Fitness.nb ~ Clustering.Six,
  weights = Weight.nb.sum
)
summary(clustering.fit.lm)
CamSol.fit.lm <- lm(
  data = by_cluster,
  formula = Fitness.nb ~ CamSol.avg,
  weights = Weight.nb.sum
)
summary(CamSol.fit.lm)
Tango.fit.lm <- lm(
  data = by_cluster,
  formula = Fitness.nb ~ TangoAAsInAPRs,
  weights = Weight.nb.sum
)
summary(Tango.fit.lm)
disorder.fit.lm <- lm(
  data = by_cluster,
  formula = Fitness.nb ~ disorder,
  weights = Weight.nb.sum
)
summary(disorder.fit.lm)
fit.aa.lm <- lm(data = by_cluster,
                formula = Fitness.nb ~ AA.fit,
                weights = Weight.nb.sum)
summary(fit.aa.lm)

delta.isd.lm <- lm(
  data = by_cluster,
  formula = Fitness.nb ~ ISD.delta,
  weights = Weight.nb.sum
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
  Slope = lm(formula = Fitness.nb ~ Clustering.Six, weights = Weight.nb.5.7)$coefficients[2],
  Sign = math.sign(Slope)
) %>% group_by(Sign) %>% summarise(N = n())

# Plotting the data.
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

by_cluster <- by_cluster %>% mutate(
  "IUPred (square root)" = sqrt(ISD.iupred2),
  "CamSol" = CamSol.avg,
  "#Aggregation-Prone AAs (Tango)" = TangoAAsInAPRs,
  "Disorder propensity" = disorder,
  "AA-predicted fitness" = AA.fit,
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
    "Fitness" = by_cluster[["Fitness.nb"]],
    "Fit" = by_cluster[[.x]],
    "Weight.nb.sum" = by_cluster[["Weight.nb.sum"]]
  )
  
  fit.plot <-
    ggplot(
      data = plot.df,
      aes(
        y = Fitness,
        x = Fit,
        size = Weight.nb.sum,
        weight = Weight.nb.sum
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


fit.zoomout.plots <- lapply(fit.list, fit.plotting, y.limit = c(0, 5.3))

names(fit.zoomout.plots) <- c("ISD", "CamSol", "Tango", "Disorder", "AA", "Clustering")

mapply(ggsave,
       file = paste0("Scripts/Figures/zoomout_fitness_", names(fit.zoomout.plots), "_", todays.date, ".png"),
       plot = fit.zoomout.plots,
       height = 7, width = 7, units = "in", dpi = 300)
