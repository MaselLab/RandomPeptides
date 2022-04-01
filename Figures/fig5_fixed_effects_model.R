# Compaing the predictive power of GC content to that of amino acid composition.

# Packages.
library(lme4)
library(tidyverse)
library(stringr)
library(Hmisc)

# Load peptide data.
peptide.data <- read.table(file = "Scripts/RandomPeptides/Data/supplemental_dataset_1.tsv", header = T, stringsAsFactors = F)
peptide.data

# Combining the data by cluster using the max weight peptide from each cluster.
cluster_max <- peptide.data %>% group_by(Cluster) %>% filter(WEIGHT == max(WEIGHT)) %>% ungroup()

# Comparing GC content as a predictor to amino acid composition.
peptide.mixed.intercept <- lm(
  data = cluster_max,
  formula = FITNESS ~ 1 ,
  weights = WEIGHT
)

peptide.mixed.gc <- lm(
  data = cluster_max[!is.na(cluster_max$GC.avg),],
  formula = FITNESS ~ GC.avg ,
  weights = WEIGHT
)
summary(peptide.mixed.gc)
drop1(peptide.mixed.gc, test = "Chisq")
peptide.mixed.aaonly.gc <- lm(
  data = cluster_max[!is.na(cluster_max$GC.avg),],
  # contrasts = c("Leu" = -0.045, "Pro" = 0.068, "Met" = -0.084, "Trp" = -0.0013, "Ala" = 0.059,
  #           "Val" = -0.051, "Phe" = -0.11, "Ile" = -0.16, "Gly" = 0.017, "Ser" = 0.035,
  #           "Thr" = -0.0094, "Cys" = -0.024, "Asn" = -0.047, "Gln" = 0.0053, "Tyr" = -0.077,
  #           "His" = -0.075, "Asp" = 0.0093, "Glu" = -0.029, "Lys" = -0.09, "Arg" = -0.043,
  #           "GC.avg" = 10),
  formula = FITNESS ~
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys + Arg +
    GC.avg +
    
    0,
  weights = WEIGHT
)
summary(peptide.mixed.aaonly.gc)
anova(peptide.mixed.aaonly.gc, peptide.mixed.gc, test = "LRT")
drop1(peptide.mixed.aaonly.gc, test = "Chisq")

peptide.mixed.nb.aaonly.lm <- lm(
  data = cluster_max[!is.na(cluster_max$GC.avg),],
  formula = FITNESS ~
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys + Arg +
    
    0,
  weights = WEIGHT
)
summary(peptide.mixed.nb.aaonly.lm)
drop1(peptide.mixed.nb.aaonly.lm, test = "Chisq")
anova(peptide.mixed.aaonly.gc, peptide.mixed.nb.aaonly.lm, test = "LRT")
#anova(peptide.mixed.nb.aaonly.lm, peptide.mixed.intercept, test = "LRT")

# Combined model.
peptide.mixed.aa.gc <- lm(
  data = cluster_max[!is.na(cluster_max$GC.avg),],
  formula = FITNESS ~ GC.avg +
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys + Arg +
     0,
  weights = WEIGHT
)
summary(peptide.mixed.aa.gc)
drop1(peptide.mixed.aa.gc, test = "Chisq")

# Checking just fixed effects models. This makes the comparison much simpler, and I don't need to
# worry about the fact that the mixed models are fitted numerically and might be stuck in some
# local optimum.
aa.lm <- lm(
  data = cluster_max,
  formula = FITNESS ~ 
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys #+ Arg +
  #0
  ,
  weights = WEIGHT
)
aa.lm.summary <- summary(aa.lm)
aa.lm.summary

gc.lm <- lm(
  data = cluster_max,
  formula = FITNESS ~ GC.avg,
  weights = WEIGHT
)
gc.lm.summary <- summary(gc.lm)
gc.lm.summary

aa.gc.lm <- lm(
  data = cluster_max,
  formula = FITNESS ~ GC.avg +
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys #+ Arg +
  #0
  ,
  weights = WEIGHT
)
aa.gc.lm.summary <- summary(aa.gc.lm)
aa.gc.lm.summary

# Comparing the GC content only, aa comp only, and GC content + aa comp models to intercept only.
intercept.lm <- lm(
  data = cluster_max,
  formula = FITNESS ~ 1,
  weights = WEIGHT
)
summary(intercept.lm)

aa.lrt <- anova(intercept.lm, aa.lm, test = "LRT")
aa.lrt
aa.lrt$`Pr(>Chi)`
gc.lrt <- anova(intercept.lm, gc.lm, test = "LRT")
gc.lrt$`Pr(>Chi)`
aa.gc.lrt <- anova(intercept.lm, aa.gc.lm, test = "LRT")
aa.gc.lrt$`Pr(>Chi)`

# Seeing if GC content improves the aa comp only model, and vice versa.
aa.gc.gc.lrt <- anova(gc.lm, aa.gc.lm, test = "LRT")
aa.gc.gc.lrt
aa.gc.aa.lrt <- anova(aa.lm, aa.gc.lm, test = "LRT")
aa.gc.aa.lrt

# Splitting the data into low and high GC content and repeating figure 1D.
hist(cluster_max$GC.avg)
summary(cluster_max$GC.avg)
gc.median <- median(cluster_max$GC.avg, na.rm = T)

high.gc.clusters <- cluster_max %>% filter(GC.avg > gc.median) %>% select(Cluster)

highgc.aaonly.lm <- lm(
  data = cluster_max %>% filter(Cluster %in% high.gc.clusters$Cluster),
  formula = FITNESS ~
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys,
  weights = WEIGHT
)
summary(highgc.aaonly.lm)
drop1(highgc.aaonly.lm, test = "Chisq")

# Testing p-value compared to intercept-only.
highgc.intercept.lm <- lm(
  data = cluster_max %>% filter(Cluster %in% high.gc.clusters$Cluster),
  formula = FITNESS ~
    1,
  weights = WEIGHT
)
anova(highgc.intercept.lm, highgc.aaonly.lm, test = "LRT")
# Checking variability in fitness for high gc-content.
cluster_max %>% filter(Cluster %in% high.gc.clusters$Cluster) %>%
  ggplot(
    data = .,
    aes(x = FITNESS)
  ) + geom_histogram() + theme_bw()

# Checking low GC content.
lowgc.aaonly.lm <- lm(
  data = cluster_max %>% filter(!(Cluster %in% high.gc.clusters$Cluster)),
  formula = FITNESS ~
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys,
  weights = WEIGHT
)
summary(lowgc.aaonly.lm)
drop1(lowgc.aaonly.lm, test = "Chisq")

lowgc.intercept.lm <- lm(
  data = cluster_max %>% filter(!(Cluster %in% high.gc.clusters$Cluster)),
  formula = FITNESS ~
    1,
  weights = WEIGHT
)
anova(lowgc.intercept.lm, lowgc.aaonly.lm, test = "LRT")
# Checking variability in fitness in low gc-content peptides.
cluster_max %>% filter(!(Cluster %in% high.gc.clusters$Cluster)) %>%
  ggplot(
    data = .,
    aes(x = FITNESS)
  ) +
  geom_histogram() + theme_bw()

# Predicting from the two models.
cluster_max <- cluster_max %>% mutate(
  "GC Type" = if_else(
    Cluster %in% high.gc.clusters$Cluster,
    "High GC (> 57.4%)", "Low GC (\U2264 57.4%)"
  )
)

gc.high <- cluster_max %>% filter(Cluster %in% high.gc.clusters$Cluster)
gc.low <- cluster_max %>% filter(!(Cluster %in% high.gc.clusters$Cluster))

gc.high <- gc.high %>% mutate(fit.highlow.gc = predict.lm(highgc.aaonly.lm, newdata = ., type = "response"))
gc.low <- gc.low %>% mutate(fit.highlow.gc = predict.lm(lowgc.aaonly.lm, newdata = ., type = "response"))

gc.cluster <- rbind(gc.high, gc.low)

cbbPalette <- c("#0072B2", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7")

highgc.lm <- lm(
  data = gc.cluster %>% filter(`GC Type` == "High GC (> 57.4%)"),
  formula = FITNESS ~ fit.highlow.gc,
  weights = WEIGHT
)
summary(highgc.lm)
drop1(highgc.lm, test = "Chisq")
lowgc.lm <- lm(
  data = gc.cluster %>% filter(`GC Type` != "High GC (> 57.4%)"),
  formula = FITNESS ~ fit.highlow.gc,
  weights = WEIGHT
)
summary(lowgc.lm)
drop1(lowgc.lm, test = "Chisq")

gc.cluster %>% ggplot(
  data = .,
  aes(
    x = fit.highlow.gc,
    y = FITNESS,
    size = WEIGHT,
    weight = WEIGHT,
    color = `GC Type`
  )
) +
  geom_smooth(method = "lm", lwd = 1.5, se = F) +
  geom_point(alpha = 0.4) + ylim(c(0,2)) +
  scale_x_continuous(breaks = c(0.25, 0.50, 0.75, 1.00)) +
  ylab("Fitness") + xlab("AA-predicted fitness") +
  scale_color_manual(values = cbbPalette) +
  theme_bw(base_size = 28)

ggsave(filename = paste0("Scripts/Figures/gc_high_low_", Sys.Date(), ".png"),
       units = "in", width = 10, height = 6)
