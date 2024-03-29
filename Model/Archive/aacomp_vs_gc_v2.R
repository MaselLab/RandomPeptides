# Compaing the predictive power of GC content to that of amino acid composition.

# Packages.
library(lme4)
library(tidyverse)
library(stringr)
library(Hmisc)

# Load peptide data.
peptide.data <- read.table(file = "RandomPeptides/Data/supplemental_table_1.tsv", header = T, stringsAsFactors = F)
peptide.data

# Comparing GC content as a predictor to amino acid composition.
peptide.mixed.intercept <- lmer(
  data = peptide.data,
  formula = FITNESS ~ 1 + (1|Cluster),
  weights = WEIGHT
)

peptide.mixed.gc <- lmer(
  data = peptide.data[!is.na(peptide.data$GC.avg),],
  formula = FITNESS ~ GC.avg +
    (1|Cluster),
  weights = WEIGHT
)
summary(peptide.mixed.gc)
peptide.mixed.aaonly.gc <- lmer(
  data = peptide.data[!is.na(peptide.data$GC.avg),],
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
    (1|Cluster) +
    0,
  weights = WEIGHT
)
summary(peptide.mixed.aaonly.gc)
anova(peptide.mixed.aaonly.gc, peptide.mixed.gc, test = "LRT")
drop1(peptide.mixed.aaonly.gc, test = "Chisq")

peptide.mixed.nb.aaonly.lm <- lmer(
  data = peptide.data[!is.na(peptide.data$GC.avg),],
  formula = FITNESS ~
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys + Arg +
    (1|Cluster) +
    0,
  weights = WEIGHT
)
summary(peptide.mixed.nb.aaonly.lm)
drop1(peptide.mixed.nb.aaonly.lm, test = "Chisq")
anova(peptide.mixed.aaonly.gc, peptide.mixed.nb.aaonly.lm, test = "LRT")
#anova(peptide.mixed.nb.aaonly.lm, peptide.mixed.intercept, test = "LRT")

# Combined model.
peptide.mixed.aa.gc <- lmer(
  data = peptide.data[!is.na(peptide.data$GC.avg),],
  formula = FITNESS ~ GC.avg +
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys + Arg +
    (1|Cluster) + 0,
  weights = WEIGHT
)
summary(peptide.mixed.aa.gc)
drop1(peptide.mixed.aa.gc, test = "Chisq")

# And now checking weighted R^2 for fitness ~ predicted fitness.
peptide.data$fit.aa <- predict(peptide.mixed.nb.aaonly.lm,
                               newdata = peptide.data,
                               type = "response",
                               random.only = F,
                               re.form = NA)
peptide.data$fit.gc <- predict(peptide.mixed.gc,
                               newdata = peptide.data,
                               type = "response",
                               random.only = F,
                               re.form = NA)
peptide.data$fit.aa.gc <- predict(peptide.mixed.aa.gc,
                                  newdata = peptide.data,
                                  type = "response",
                                  random.only = F,
                                  re.form = NA)
peptide.data[is.na(peptide.data$fit.gc), "PeptideID"]
peptide.data[is.na(peptide.data$fit.aa.gc), "PeptideID"]

peptide.cluster <- 
  peptide.data %>%
  group_by(Cluster) %>%
  summarise(Weight.nb.sum = sum(WEIGHT), FITNESS = wtd.mean(FITNESS, weights = WEIGHT),
            fit.aa = wtd.mean(fit.aa, weights = WEIGHT), fit.gc = wtd.mean(fit.gc, weights = WEIGHT, na.rm = T),
            fit.aa.gc = wtd.mean(fit.aa.gc, weights = WEIGHT, na.rm = T),
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
            GC.avg = wtd.mean(GC.avg, weights = WEIGHT, na.rm = T),
            Clustering.Six = wtd.mean(Clustering.Six, weights = WEIGHT),
            net.charge = wtd.mean(net.charge, weights = WEIGHT),
            TangoAAsInAPRs = wtd.mean(TangoAAsInAPRs, weights = WEIGHT))
peptide.cluster

# Predicting fitness with predicted fitness to get R^2 values.
predfit.aa <- lm(
  data = peptide.cluster,
  formula = FITNESS ~ fit.aa,
  weights = Weight.nb.sum
)
summary(predfit.aa)

predfit.gc <- lm(
  data = peptide.cluster,
  formula = FITNESS ~ fit.gc,
  weights = Weight.nb.sum
)
summary(predfit.gc)

predfit.aa.gc <- lm(
  data = peptide.cluster,
  formula = FITNESS ~ fit.aa.gc,
  weights = Weight.nb.sum
)
summary(predfit.aa.gc)

# Checking just fixed effects models. This makes the comparison much simpler, and I don't need to
# worry about the fact that the mixed models are fitted numerically and might be stuck in some
# local optimum.
aa.lm <- lm(
  data = peptide.cluster,
  formula = FITNESS ~ 
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys #+ Arg +
  #0
  ,
  weights = Weight.nb.sum
)
aa.lm.summary <- summary(aa.lm)
aa.lm.summary

gc.lm <- lm(
  data = peptide.cluster,
  formula = FITNESS ~ GC.avg,
  weights = Weight.nb.sum
)
gc.lm.summary <- summary(gc.lm)
gc.lm.summary

aa.gc.lm <- lm(
  data = peptide.cluster,
  formula = FITNESS ~ GC.avg +
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys #+ Arg +
  #0
  ,
  weights = Weight.nb.sum
)
aa.gc.lm.summary <- summary(aa.gc.lm)
aa.gc.lm.summary

# Comparing the GC content only, aa comp only, and GC content + aa comp models to intercept only.
intercept.lm <- lm(
  data = peptide.cluster,
  formula = FITNESS ~ 1,
  weights = Weight.nb.sum
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
hist(peptide.cluster$GC.avg)
summary(peptide.cluster$GC.avg)
gc.median <- median(peptide.cluster$GC.avg, na.rm = T)

high.gc.clusters <- peptide.cluster %>% filter(GC.avg > gc.median) %>% select(Cluster)

highgc.aaonly.lm <- lmer(
  data = peptide.data %>% filter(Cluster %in% high.gc.clusters$Cluster),
  formula = FITNESS ~
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys + Arg +
    (1|Cluster) +
    0,
  weights = WEIGHT
)
summary(highgc.aaonly.lm)
drop1(highgc.aaonly.lm, test = "Chisq")

lowgc.aaonly.lm <- lmer(
  data = peptide.data %>% filter(!(Cluster %in% high.gc.clusters$Cluster)),
  formula = FITNESS ~
    Leu + Pro + Met + Trp + Ala +
    Val + Phe + Ile + Gly + Ser +
    Thr + Cys + Asn + Gln + Tyr +
    His + Asp + Glu + Lys + Arg +
    (1|Cluster) +
    0,
  weights = WEIGHT
)
summary(lowgc.aaonly.lm)
drop1(lowgc.aaonly.lm, test = "Chisq")

peptide.data <- peptide.data %>% mutate(
  fit.highlow.gc = if_else(
    Cluster %in% high.gc.clusters$Cluster,
    predict(highgc.aaonly.lm,
            newdata = .,
            type = "response",
            random.only = F,
            re.form = NA),
    predict(lowgc.aaonly.lm,
            newdata = .,
            type = "response",
            random.only = F,
            re.form = NA)
  )
)

gc.cluster <- peptide.data %>%
  group_by(Cluster) %>%
  summarise(
    Weight.nb.sum = sum(WEIGHT), FITNESS = wtd.mean(FITNESS, weights = WEIGHT),
    fit.highlow.gc = wtd.mean(fit.highlow.gc, weights = WEIGHT)
  )

gc.cluster <- gc.cluster %>% mutate(
  "GC Type" = if_else(
    Cluster %in% high.gc.clusters$Cluster,
    "High GC (> 57.4%)", "Low GC (\U2264 57.4%)"
  )
)

cbbPalette <- c("#0072B2", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7")

highgc.lm <- lm(
  data = gc.cluster %>% filter(`GC Type` == "High GC (> 57.4%)"),
  formula = FITNESS ~ fit.highlow.gc,
  weights = Weight.nb.sum
)
summary(highgc.lm)
lowgc.lm <- lm(
  data = gc.cluster %>% filter(`GC Type` != "High GC (> 57.4%)"),
  formula = FITNESS ~ fit.highlow.gc,
  weights = Weight.nb.sum
)
summary(lowgc.lm)

gc.cluster %>% ggplot(
  data = .,
  aes(
    x = fit.highlow.gc,
    y = FITNESS,
    size = Weight.nb.sum,
    weight = Weight.nb.sum,
    color = `GC Type`
  )
) +
  geom_smooth(method = "lm", lwd = 1.5, se = F) +
  geom_point(alpha = 0.4) + ylim(c(0,2)) +
  scale_x_continuous(breaks = c(0.25, 0.50, 0.75, 1.00)) +
  ylab("Fitness") + xlab("AA-predicted fitness") +
  scale_color_manual(values = cbbPalette) +
  theme_bw(base_size = 28)

ggsave(filename = "Figures/gc_high_low.png",
       units = "in", width = 10, height = 6)
