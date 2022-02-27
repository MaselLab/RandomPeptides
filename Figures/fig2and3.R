# Figure 3 script.

# Load packages
library(tidyverse)
library(lme4)
library(Hmisc)
library(Biostrings)

# Load data.
peptide.data <- read.table(file = "Scripts/RandomPeptides/Data/supplemental_dataset_1.tsv", header = T, stringsAsFactors = F)

# Date of file generation.
today.date <- Sys.Date()

# Calculating the marginal effects.
# Amino acid frequency only model.
peptide.mixed.nb.log.aaonly.lm <-
  lmer(
    data = peptide.data,
    formula = FITNESS ~
      Leu + Pro + Met + Trp + Ala +
      Val + Phe + Ile + Gly + Ser +
      Thr + Cys + Asn + Gln + Tyr +
      His + Asp + Glu + Lys + Arg +
      (1|Cluster) +
      0,
    weights = WEIGHT
  )
log.aaonly.summary <- summary(peptide.mixed.nb.log.aaonly.lm)
log.aaonly.summary
drop1(peptide.mixed.nb.log.aaonly.lm, test = "Chisq")

# Estimating the marginal effects.
custom.margins.exp.var <- function(weight.vector, linear.model){
  # First, we calculate the expected value (the mean). Note that beta coefficients are unbiased estimators.
  model.summary <- summary(linear.model)
  #print(model.summary)
  model.coefs <- model.summary$coefficients[,1]
  #print(model.coefs)
  mean.margin <- sum(weight.vector * model.coefs)
  #print(mean.margin)
  
  # Next, we calculate the variance.
  weight.matrix <- outer(weight.vector, weight.vector)
  cov.matrix <- vcov(linear.model)
  full.matrix <- weight.matrix * cov.matrix
  var.margin <- sum(full.matrix)
  std.dv <- sqrt(var.margin)
  
  # Now, both are returned as a vector.
  return(c(mean.margin, std.dv))
}

# Writing a second function to increment the positions of the elements in a vector. So, element 1 should
# become element 2, 2 becomes 3, etc up to n becomes 1.
increment.element.position <- function(any.vector){
  vector.length <- length(any.vector)
  vector.matrix <- matrix(nrow = vector.length, ncol = vector.length)
  vector.matrix[1,] <- any.vector
  for (i in 2:vector.length) {
    vector.matrix[i, ] <- any.vector[c((vector.length + 2 - i):vector.length, 1:(vector.length + 1 - i))]
  }
  return(vector.matrix)
}

# Function for building a weights matrix where the weights for each row vary depending on the
# expected frequency/counts of the focal coefficient.
increment.position.freq <- function(counts.vector){
  total.counts <- sum(counts.vector)
  vector.length <- length(counts.vector)
  vector.matrix <- matrix(nrow = vector.length, ncol = vector.length)
  for (i in 1:vector.length) {
    for (j in 1:vector.length) {
      if (j == i) {
        vector.matrix[i, j] <- 1
      } else {
        vector.matrix[i, j] <- -counts.vector[j] / (total.counts - counts.vector[i])
      }
    }
  }
  return(vector.matrix)
}

# Function for ordering a data frame with column "AminoAcid" by the order of amino acids
# in a linear model.
aa.sort.lm <- function(linear.model, aa.df, aa.col = "AminoAcid") {
  lm.summary <- summary(linear.model)
  #print(lm.summary)
  aa.order <- names(lm.summary$coefficients[,1])
  #print(aa.order)
  aa.unordered <- as.list(aa.df[, aa.col])[[1]]
  #print(aa.unordered)
  new.order <- match(aa.order, aa.unordered)
  #print(new.order)
  aa.df.ordered <- aa.df[new.order,]
  return(aa.df.ordered)
}

# Writing a third function to get all necessary combinations of the weight vector, and use those to
# call the above function, all while keeping the coefficient names with the correct output values.
# Should return everything as a data frame with coefficient name, mean marginal effect, and variance
# as its three columns.
# Note: weight.df needs to have a column to sort the amino acids by!
custom.contrasts <- function(weight.df, linear.model, aa.col = "AminoAcid", weight.col = "Count"){
  weight.df.ordered <- aa.sort.lm(linear.model = linear.model, aa.df = weight.df, aa.col = aa.col)
  weight.vector <- as.list(weight.df.ordered[, weight.col])[[1]]
  coef.count <- length(weight.vector)
  model.summary <- summary(linear.model)
  names.vector <- names(model.summary$coefficients[,1])
  mean.var.df <- data.frame("Coef" = names.vector,
                            "Mean" = c(rep(NA, coef.count)),
                            "Std.Err" = c(rep(NA, coef.count)))
  full.weight.matrix <- increment.position.freq(weight.vector)
  for (i in 1:coef.count) {
    mean.var.df[i, c(2,3)] <- custom.margins.exp.var(full.weight.matrix[i,], linear.model)
  }
  return(mean.var.df)
}

# Grouping all data by cluster, and taking the weighted mean of each predicter for each cluster.
by_cluster <-
  peptide.data %>% 
  group_by(Cluster) %>%
  summarise(Weight.nb.sum = sum(WEIGHT), 
            ISD.iupred2 = wtd.mean(ISD.iupred2, weights = WEIGHT), FITNESS = wtd.mean(FITNESS, weights = WEIGHT),
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
            Trp.unweighted = mean(Trp), Arg.unweighted = mean(Arg)
  )

# Calculating frequency from cluster.
aa.freqs.wmax <-
  tibble(
    "AminoAcid" = c("Lys", "Arg",
                    "Asp", "Glu",
                    "Gln", "Ser",
                    "Asn", "Thr",
                    "Pro", "Tyr",
                    "His", "Cys",
                    "Gly", "Ala",
                    "Trp", "Val",
                    "Leu", "Met",
                    "Phe", "Ile"),
    "Count" = c("Lys" = sum(by_cluster$Lys), "Arg" = sum(by_cluster$Arg),
                "Asp" = sum(by_cluster$Asp), "Glu" = sum(by_cluster$Glu),
                "Gln" = sum(by_cluster$Gln), "Ser" = sum(by_cluster$Ser),
                "Asn" = sum(by_cluster$Asn), "Thr" = sum(by_cluster$Thr),
                "Pro" = sum(by_cluster$Pro), "Tyr" = sum(by_cluster$Tyr),
                "His" = sum(by_cluster$His), "Cys" = sum(by_cluster$Cys),
                "Gly" = sum(by_cluster$Gly), "Ala" = sum(by_cluster$Ala),
                "Trp" = sum(by_cluster$Trp), "Val" = sum(by_cluster$Val),
                "Leu" = sum(by_cluster$Leu), "Met" = sum(by_cluster$Met),
                "Phe" = sum(by_cluster$Phe), "Ile" = sum(by_cluster$Ile))
  )
aa.freqs.wmax

# The order of the freq data and the linear model need to be the same.
match(names(log.aaonly.summary$coefficients[,1]), aa.freqs.wmax$AminoAcid)
aa.freqs.wmax[match(names(log.aaonly.summary$coefficients[,1]), aa.freqs.wmax$AminoAcid),]
aa.freqs.wmax[match(names(log.aaonly.summary$coefficients[,1]), aa.freqs.wmax[,1]),]
aa.freqs.ordered <- aa.sort.lm(peptide.mixed.nb.log.aaonly.lm, aa.freqs.wmax)
aa.freqs.ordered
as.list(aa.freqs.ordered[,"Count"])[[1]]

marginals.clusters.log.nb.wmax <- custom.contrasts(aa.freqs.wmax, peptide.mixed.nb.log.aaonly.lm)
marginals.clusters.log.nb.wmax

# Sorting the table by ascending marginal effect.
marginals.ordered <- arrange(marginals.clusters.log.nb.wmax, Mean)

# Exporting marginal effects as a table.
write_tsv(marginals.ordered, path = "Scripts/RandomPeptides/Data/supplemental_table_2.tsv")

# Importing data for James et al.'s phylostratigraphy slopes.
# Genes first.
gene.slopes <- read_csv(file = "Data/AAcomp_Fullgene_Nontrans.csv")
gene.slopes.only <- gene.slopes[, c(2,3:4,6:7,15:16)]
gene.slopes.only
names(gene.slopes.only) <- c("AA", "AnimalYoungGeneSlope", "AnimalYoungGeneError",
                             "PlantYoungGeneSlope", "PlantYoungGeneError",
                             "AncientGeneSlope", "AncientGeneError")
gene.slopes.only

# And now Pfams.
pfam.slopes <- read_csv(file = "Data/AAcomp_Pfam_Nontrans.csv")
pfam.slopes.only <- pfam.slopes[, c(1:4, 6:7, 15:16)]
names(pfam.slopes.only)
names(pfam.slopes.only) <- c("OneLetter", "AA", "AnimalYoungPfamSlope", "AnimalYoungPfamError",
                             "PlantYoungPfamSlope", "PlantYoungPfamError",
                             "AncientPfamSlope", "AncientPfamError")
pfam.slopes.only

slopes.aa.props <- merge(gene.slopes.only, pfam.slopes.only, by = "AA")
slopes.aa.props

# Renaming the columns of the marginals data frame to match those from James et al.'s data frame.
names(marginals.clusters.log.nb.wmax)
colname.marginals <- "MarginalLogNBWMax"
colname.stderr <- "MarginalLogNBWMaxErr"
names(marginals.clusters.log.nb.wmax) <- c("AA", colname.marginals, colname.stderr)
names(marginals.clusters.log.nb.wmax)
marginals.clusters.log.nb.wmax

# Merging the James et al. data frame and the marginal effects data frames.
slopes.aa.props <- merge(slopes.aa.props, marginals.clusters.log.nb.wmax, by = "AA")
slopes.aa.props

# Looking at the errors to compare.
slopes.aa.props[, c("AnimalYoungPfamError", "PlantYoungPfamError", "AncientPfamError")]

# Adding size and amino acid cost.
source(file = "Scripts/RandomPeptides/Metrics/aa_comp_metrics.R")
properties.df <- tibble(
  "OneLetter" = slopes.aa.props$OneLetter#,
  # "Size" = c(88.6, 173.4, 114.1, 111.1, 108.5, 143.8, 138.4, 60.1, 153.2, 166.7,
  #            166.7, 168.6, 162.9, 189.9, 112.7, 89.0, 116.1, 227.8, 193.6, 140.0)
  # ,
  # "Cost-Ecoli" = c(11.7, 27.3, 14.7, 12.7, 24.7, 16.3, 15.3, 11.7, 38.3, 32.3,
  #                  27.3, 30.3, 34.3, 52.0, 20.3, 11.7, 18.7, 74.3, 50.0, 23.3),
  # "RSA" = c(0.218, 0.361, 0.342, 0.366, 0.101, 0.364, 0.411, 0.269, 0.269, 0.125,
  #           0.147, 0.446, 0.159, 0.136, 0.342, 0.278, 0.272, 0.163, 0.187, 0.143),
  # "DisorderPropensity" = c(0.450, 0.394, 0.285, 0.407, 0.000, 0.665, 0.781, 0.437, 0.259, 0.090,
  #                          0.195, 0.588, 0.291, 0.117, 1.000, 0.713, 0.401, 0.004, 0.113, 0.263)
)
properties.df
properties.df$Size <- mean.metric.calculator(properties.df$OneLetter, metric = "volume")
properties.df$CostEcoli <- mean.metric.calculator(properties.df$OneLetter, metric = "cost-ecoli")
properties.df$RSA <- mean.metric.calculator(properties.df$OneLetter, metric = "RSA-hydrophilicity")
properties.df$DisorderPropensity <- mean.metric.calculator(properties.df$OneLetter, metric = "disorder")
properties.df$stickiness <- mean.metric.calculator(properties.df$OneLetter, metric = "stickiness")
properties.df$weight <- mean.metric.calculator(properties.df$OneLetter, metric = "weight")
properties.df$pI <- mean.metric.calculator(properties.df$OneLetter, metric = "pI")
properties.df

# Also calculating the frequencies in the E.coli proteome, based off the K-12 reference.
ecoli <- readAAStringSet(filepath = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000000625/UP000000625_83333.fasta.gz")
ecoli_freqs <- alphabetFrequency(ecoli, as.prob = F)
ecoli_totals <- colSums(ecoli_freqs)
ecoli_freqs_probs <- ecoli_totals / sum(ecoli_totals)

# Making sure the amino acids match.
identical(names(ecoli_freqs_probs)[1:20], properties.df$OneLetter)
properties.df$FreqsEcoli <- ecoli_freqs_probs[1:20]

# Note that the frequencies don't quite sum to 1 because there are a few unknown amino acids.
sum(properties.df$FreqsEcoli)

# Adjusting for codon number
properties.df$Codon_freqs <-
  c(
    4/61,
    6/61,
    2/61,
    2/61,
    2/61,
    2/61,
    2/61,
    4/61,
    2/61,
    3/61,
    6/61,
    2/61,
    1/61,
    2/61,
    4/61,
    6/61,
    4/61,
    1/61,
    2/61,
    4/61
  )
# Making sure codon frequencies sum to 1
sum(properties.df$Codon_freqs)
# Double checking values to make sure they're right.
properties.df %>% select(OneLetter, Codon_freqs)
# Taking the difference in E. coli AA freqs and codon freqs to get the difference than expected.
properties.df$Diff_freqs <- properties.df$FreqsEcoli - properties.df$Codon_freqs
# Checking all three together
properties.df %>% select(OneLetter, FreqsEcoli, Codon_freqs, Diff_freqs)

cor.test(properties.df$Size, properties.df$CostEcoli, method = "spearman")
cor.test(properties.df$Size, properties.df$DisorderPropensity, method = "spearman")
cor.test(properties.df$CostEcoli, properties.df$DisorderPropensity, method = "spearman")

cor.test(properties.df$Size, properties.df$CostEcoli, method = "pearson")
cor.test(properties.df$Size, properties.df$DisorderPropensity, method = "pearson")
cor.test(properties.df$CostEcoli, properties.df$DisorderPropensity, method = "pearson")

cov(properties.df[, -1])
cor(properties.df[, -1], method = "pearson")
cor(properties.df[, -1], method = "spearman")

slopes.aa.props <- merge(slopes.aa.props, properties.df, by = "OneLetter")
slopes.aa.props

# Adding in a variable for highly hydrophobic, where I, L, M, F, and V are considered highly hydrophobic.
# Note: This is deprecated. These results are very sensitive to whether W is considered hydrophobic.
slopes.aa.props$hydrophobic <- ifelse(slopes.aa.props$OneLetter %in% c("I", "L", "M", "F", "V"), 1, 0)

marginals.hydrophobic.lm <- lm(
  data = slopes.aa.props,
  formula = MarginalLogNBWMax ~ Size + DisorderPropensity + hydrophobic,
  weights = 1 / (MarginalLogNBWMaxErr^2)
)
summary(marginals.hydrophobic.lm)
drop1(marginals.hydrophobic.lm, test = "Chisq")

# Checking other predictors.
cor(slopes.aa.props %>% select(DisorderPropensity, Size, CostEcoli,
                               stickiness, RSA, weight,
                               pI, FreqsEcoli, Codon_freqs,
                               Diff_freqs, hydrophobic))

marginals.lm <- lm(
  data = slopes.aa.props,
  formula = MarginalLogNBWMax ~ DisorderPropensity
  #+ hydrophobic
  + Size
  #+ CostEcoli
  #+ stickiness
  #+ RSA
  #+ weight
  #+ pI
  #+ FreqsEcoli
  #+ Codon_freqs
  #+ Diff_freqs
  ,
  weights = 1 / (MarginalLogNBWMaxErr^2)
)
summary(marginals.lm)
drop1(marginals.lm, test = "Chisq")

# Comparing to an RSA model.
rsa.lm <- lm(
  data = slopes.aa.props,
  formula = MarginalLogNBWMax ~ #DisorderPropensity
  #+ hydrophobic
  + Size
  #+ CostEcoli
  #+ stickiness
  + RSA
  #+ weight
  #+ pI
  #+ FreqsEcoli
  #+ Codon_freqs
  #+ Diff_freqs
  ,
  weights = 1 / (MarginalLogNBWMaxErr^2)
)
summary(rsa.lm)
drop1(rsa.lm, test = "Chisq")

intercept.lm <- lm(
  data = slopes.aa.props,
  formula = MarginalLogNBWMax ~ 1
  ,
  weights = 1 / (MarginalLogNBWMaxErr^2)
)
summary(intercept.lm)

anova(intercept.lm, rsa.lm, test = "LRT")
anova(intercept.lm, marginals.lm)

# Checking weighted correlations.
# First, making a function for a weighted Pearson's correlation where both X and Y have their own
# sources of error.
custom.weighted.pearson <- function(x, y, x.var, y.var){
  require(weights)
  x.z <- z.w.transform(x, x.var)
  y.z <- z.w.transform(y, y.var)
  weights.z <- 1/(((x.var^2)/(sd.weighted(x, x.var)^2)) +
                    ((y.var^2)/(sd.weighted(y, y.var)^2)))
  w.pearson <- wtd.cor(x.z, y.z, weight = weights.z)
  return(w.pearson)
}

# Standard normal transform using point associated standard errors.
mean.weighted <- function(x, x.err){
  w.x <- (1/(x.err^2))/sum(1/(x.err^2))
  mean.x <- sum(w.x * x)
  return(mean.x)
}
sd.weighted <- function(x, x.err){
  w.x <- (1/(x.err^2))/sum(1/(x.err^2))
  mean.x <- mean.weighted(x, x.err)
  var.x <- sum(w.x * ((x - mean.x)^2))
  sd.x <- sqrt(var.x)
  return(sd.x)
}
z.w.transform <- function(x, x.err){
  mean.x <- mean.weighted(x, x.err)
  sd.x <- sd.weighted(x, x.err)
  z <- (x - mean.x)/sd.x
  return(z)
}

# Calculating the weighted Pearson's correlation between marginal effects and phylostratigraphy slopes
# for animal, plant, and ancient Pfams.
animalpfam.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$AnimalYoungPfamSlope,
  y.var = slopes.aa.props$AnimalYoungPfamError
)
animalpfam.cor
plantpfam.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$PlantYoungPfamSlope,
  y.var = slopes.aa.props$PlantYoungPfamError
)
plantpfam.cor
ancientpfam.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$AncientPfamSlope,
  y.var = slopes.aa.props$AncientPfamError
)
ancientpfam.cor

# Adding carbons per amino acid side chain.
carbons.df <- tibble(
  "AA" = c("Ala", "Arg", "Asn", "Asp", "Cys",
           "Gln", "Glu", "Gly", "His", "Ile",
           "Leu", "Lys", "Met", "Phe", "Pro",
           "Ser", "Thr", "Trp", "Tyr", "Val"),
  "Carbons" = c(1, 4, 2, 2, 1, 3, 3, 0, 4, 4, 4, 4, 3, 7, 3, 1, 2, 9, 7, 3)
)
carbons.df

slopes.aa.props <- merge(slopes.aa.props, carbons.df, by = "AA")
slopes.aa.props

properties.carbons.df <- merge(properties.df, slopes.aa.props[, c("OneLetter", "Carbons")], by = "OneLetter")
properties.carbons.df

# Checking how carbons correlate with other properties.
cor(properties.carbons.df[, -1], method = "pearson")
cor(properties.carbons.df[, -1], method = "spearman")

carbons.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$Carbons,
  y.var = rep(1, 20)
)
carbons.cor

size.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$Size,
  y.var = rep(1, 20)
)
size.cor

cost.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$CostEcoli,
  y.var = rep(1, 20)
)
cost.cor

rsa.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$RSA,
  y.var = rep(1, 20)
)
rsa.cor

disprop.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$DisorderPropensity,
  y.var = rep(1, 20)
)
disprop.cor

weight.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$weight,
  y.var = rep(1, 20)
)
weight.cor

stickiness.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$stickiness,
  y.var = rep(1, 20)
)
stickiness.cor

pI.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$pI,
  y.var = rep(1, 20)
)
pI.cor

hydrophobicity.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$hydrophobic,
  y.var = rep(1, 20)
)
hydrophobicity.cor

ecoli_freqs.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$FreqsEcoli,
  y.var = rep(1, 20)
)
ecoli_freqs.cor

codon_freqs.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$Codon_freqs,
  y.var = rep(1, 20)
)
codon_freqs.cor

diff_freqs.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$Diff_freqs,
  y.var = rep(1, 20)
)
diff_freqs.cor

# Checking how much the model changes if tryptophan is excluded.
marginals.noW.lm <- lm(
  data = slopes.aa.props %>% filter(OneLetter != "W"),
  formula = MarginalLogNBWMax ~ Size + DisorderPropensity
  ,
  weights = 1 / (MarginalLogNBWMaxErr^2)
)
summary(marginals.noW.lm)

# Including W as a hydrophobic aminoa cid.
slopes.aa.props <- slopes.aa.props %>% mutate(
  w_hydrop = if_else(
    OneLetter == "W", 1, hydrophobic
  )
)
marginals.whydro.lm <- lm(
  data = slopes.aa.props,
  formula = MarginalLogNBWMax ~ Size + DisorderPropensity + w_hydrop
  ,
  weights = 1 / (MarginalLogNBWMaxErr^2)
)
summary(marginals.whydro.lm)

# Checking the mean and median fitness cost from Mehlhoff et al. (2020).
mehlhoff.data <- read_tsv(file = "Data/MeanMedian_Fitness.txt")
mehlhoff.data

names(mehlhoff.data)
names(mehlhoff.data)[1] <- "OneLetter"

slopes.aa.props <- merge(slopes.aa.props, mehlhoff.data, by = "OneLetter")

marginals.mehlhoff.lm <- lm(
  data = slopes.aa.props,
  formula = MarginalLogNBWMax ~ Size + DisorderPropensity + MeanFit,
  weights = 1 / (MarginalLogNBWMaxErr^2)
)
summary(marginals.mehlhoff.lm)
drop1(marginals.mehlhoff.lm, test = "Chisq")

mehlhoff.mean.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$MeanFit,
  y.var = rep(1, 20)
)
mehlhoff.mean.cor

mehlhoff.median.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$MedianFit,
  y.var = rep(1, 20)
)
mehlhoff.median.cor

properties.carbons.df <- merge(properties.carbons.df, mehlhoff.data, by = "OneLetter")
properties.carbons.df
cor(properties.carbons.df[, -1])
cor(properties.carbons.df[-c(2, 6, 13), -1], method = "spearman")

# Making the figures.
partB.name <- paste("Scripts/Figures/NonTransmembranePlantYoungPfam_MarginalLogNB_notrel_", today.date, ".png", sep = "")
plant.cor.label <- paste("R =", signif(plantpfam.cor[1], 2), "p =", signif(plantpfam.cor[4], 1))
png(filename = partB.name, width = 600, height = 600)
ggplot(data = slopes.aa.props, aes(x = MarginalLogNBWMax, y = PlantYoungPfamSlope)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin=PlantYoungPfamSlope - PlantYoungPfamError,
                    ymax=PlantYoungPfamSlope + PlantYoungPfamError),
                size = 1) +
  geom_errorbarh(aes(xmin=MarginalLogNBWMax - MarginalLogNBWMaxErr,
                     xmax=MarginalLogNBWMax + MarginalLogNBWMaxErr), size = 1) +
  ylab("\U0394% point change per BY") +
  xlab("Effect on genotype freq / cycle") +
  #scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1),
  #                   labels = round(2^c(-0.2, -0.1, 0, 0.1), digits = 2)) +
  scale_y_continuous(breaks = c(1e-5, 0, -1e-5, -2e-5, -3e-5),
                     limits = c(-3e-5, 1.8e-5),
                     labels = c("1%", "0%", "-1%", "-2%", "-3%")) +
  geom_text(aes(label = OneLetter), hjust = -0.2, vjust = -0.4, size = 9, color = "grey30") +
  annotate("text", label = plant.cor.label, x = -0.05, y = -1.5e-5, size = 10,
           parse = F) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

partA.name <- paste("Scripts/Figures/NonTransmembraneAnimalYoungPfam_MarginalLogNB_notrel_", today.date, ".png", sep = "")
animal.cor.label <- paste("R =", signif(animalpfam.cor[1], 2), "p =", signif(animalpfam.cor[4], 1))
png(filename = partA.name, width = 600, height = 600)
ggplot(data = slopes.aa.props, aes(x = MarginalLogNBWMax, y = AnimalYoungPfamSlope)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin=AnimalYoungPfamSlope - AnimalYoungPfamError,
                    ymax=AnimalYoungPfamSlope + AnimalYoungPfamError),
                size = 1) +
  geom_errorbarh(aes(xmin=MarginalLogNBWMax - MarginalLogNBWMaxErr,
                     xmax=MarginalLogNBWMax + MarginalLogNBWMaxErr), size = 1) +
  ylab("\U0394% point change per BY") +
  xlab("Effect on genotype freq / cycle") +
  #scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1),
  #                   labels = round(2^c(-0.2, -0.1, 0, 0.1), digits = 2)) +
  scale_y_continuous(breaks = c(1e-5, 0, -1e-5, -2e-5, -3e-5),
                     limits = c(-3e-5, 1.8e-5),
                     labels = c("1%", "0%", "-1%", "-2%", "-3%")) +
  geom_text(aes(label = OneLetter), hjust = -0.2, vjust = -0.4, size = 9, color = "grey30") +
  annotate("text", label = animal.cor.label, x = -0.05, y = -1.5e-5, size = 10,
           parse = F) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

partC.name <- paste("Scripts/Figures/NonTransmembraneAncientPfam_MarginalLogNB_notrel_", today.date, ".png", sep = "")
ancient.cor.label <- paste("R =", signif(ancientpfam.cor[1], 2), "p =", signif(ancientpfam.cor[4], 1))
png(filename = partC.name, width = 600, height = 600)
ggplot(data = slopes.aa.props, aes(x = MarginalLogNBWMax, y = AncientPfamSlope)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin=AncientPfamSlope - AncientPfamError,
                    ymax=AncientPfamSlope + AncientPfamError),
                size = 1) +
  geom_errorbarh(aes(xmin=MarginalLogNBWMax - MarginalLogNBWMaxErr,
                     xmax=MarginalLogNBWMax + MarginalLogNBWMaxErr), size = 1) +
  ylab("\U0394% point change per BY") +
  xlab("Effect on genotype freq / cycle") +
  #scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1),
  #                   labels = round(2^c(-0.2, -0.1, 0, 0.1), digits = 2)) +
  scale_y_continuous(breaks = c(1e-5, 0, -1e-5, -2e-5, -3e-5),
                     limits = c(-3e-5, 1.8e-5),
                     labels = c("1%", "0%", "-1%", "-2%", "-3%")) +
  geom_text(aes(label = OneLetter), hjust = -0.2, vjust = -0.4, size = 9, color = "grey30") +
  annotate("text", label = ancient.cor.label, x = -0.05, y = -1.5e-5, size = 10,
           parse = F) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

# Also checking a figure for marginals vs number of side chain carbons.
carbons.name <- paste("Scripts/Figures/carbons_vs_marginals_", today.date, ".png", sep = "")
carbons.cor.label <- paste("R =", signif(carbons.cor[1], 2), "p =", signif(carbons.cor[4], 1))
png(filename = carbons.name, width = 600, height = 600)
ggplot(data = slopes.aa.props, aes(x = MarginalLogNBWMax, y = Carbons)) +
  geom_point(size = 4) +
  #geom_errorbar(aes(ymin=AncientPfamSlope - AncientPfamError,
  #                  ymax=AncientPfamSlope + AncientPfamError),
  #              size = 1) +
  geom_errorbarh(aes(xmin=MarginalLogNBWMax - MarginalLogNBWMaxErr,
                     xmax=MarginalLogNBWMax + MarginalLogNBWMaxErr), size = 1) +
  ylab("Side chain carbons") +
  xlab("Effect on genotype freq / cycle") +
  #scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1),
  #                   labels = round(2^c(-0.2, -0.1, 0, 0.1), digits = 2)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8)) +
  geom_text(aes(label = OneLetter), hjust = -0.2, vjust = -0.4, size = 9, color = "grey30") +
  annotate("text", label = carbons.cor.label, x = -0.05, y = 8.5, size = 10,
           parse = F) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

Size.name <- paste("Scripts/Figures/Size_vs_marginals_", today.date, ".png", sep = "")
Size.cor.label <- paste("R =", signif(size.cor[1], 2), "p =", signif(size.cor[4], 1))
png(filename = Size.name, width = 600, height = 600)
ggplot(data = slopes.aa.props, aes(y = MarginalLogNBWMax, x = Size)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin=MarginalLogNBWMax - MarginalLogNBWMaxErr,
                     ymax=MarginalLogNBWMax + MarginalLogNBWMaxErr), size = 1) +
  xlab(paste("Volume ", "(", "\uc5\u00b3", ")", sep = "")) +
  ylab("Effect on genotype freq / cycle") +
  #scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1),
  #                   labels = round(2^c(-0.2, -0.1, 0, 0.1), digits = 2)) +
  #scale_y_continuous(breaks = c(0, 2, 4, 6, 8)) +
  geom_text(aes(label = OneLetter), hjust = -0.2, vjust = -0.4, size = 9, color = "grey30") +
  annotate("text", label = Size.cor.label, y = -0.05, x = 125, size = 10,
           parse = F) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

cost.name <- paste("Scripts/Figures/cost_vs_marginals_", today.date, ".png", sep = "")
cost.cor.label <- paste("R =", signif(cost.cor[1], 2), "p =", signif(cost.cor[4], 1))
png(filename = cost.name, width = 600, height = 600)
ggplot(data = slopes.aa.props, aes(y = MarginalLogNBWMax, x = CostEcoli)) +
  geom_point(size = 4) +
  #geom_errorbar(aes(ymin=AncientPfamSlope - AncientPfamError,
  #                  ymax=AncientPfamSlope + AncientPfamError),
  #              size = 1) +
  geom_errorbar(aes(ymin=MarginalLogNBWMax - MarginalLogNBWMaxErr,
                     ymax=MarginalLogNBWMax + MarginalLogNBWMaxErr), size = 1) +
  xlab("Amino acid cost in E. coli") +
  ylab("Effect on genotype freq / cycle") +
  #scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1),
  #                   labels = round(2^c(-0.2, -0.1, 0, 0.1), digits = 2)) +
  #scale_y_continuous(breaks = c(0, 2, 4, 6, 8)) +
  geom_text(aes(label = OneLetter), hjust = -0.2, vjust = -0.4, size = 9, color = "grey30") +
  annotate("text", label = cost.cor.label, y = -0.07, x = 55, size = 10,
           parse = F) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

disprop.name <- paste("Scripts/Figures/disprop_vs_marginals_", today.date, ".png", sep = "")
disprop.cor.label <- paste("R =", signif(disprop.cor[1], 2), "p =", signif(disprop.cor[4], 1))
png(filename = disprop.name, width = 600, height = 600)
ggplot(data = slopes.aa.props, aes(y = MarginalLogNBWMax, x = DisorderPropensity)) +
  geom_point(size = 4) +
  #geom_errorbar(aes(ymin=AncientPfamSlope - AncientPfamError,
  #                  ymax=AncientPfamSlope + AncientPfamError),
  #              size = 1) +
  geom_errorbar(aes(ymin=MarginalLogNBWMax - MarginalLogNBWMaxErr,
                     ymax=MarginalLogNBWMax + MarginalLogNBWMaxErr), size = 1) +
  xlab("Disorder propensity") +
  ylab("Effect on genotype freq / cycle") +
  #scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1),
  #                   labels = round(2^c(-0.2, -0.1, 0, 0.1), digits = 2)) +
  #scale_y_continuous(breaks = c(0, 2, 4, 6, 8)) +
  geom_text(aes(label = OneLetter), hjust = -0.2, vjust = -0.4, size = 9, color = "grey30") +
  annotate("text", label = disprop.cor.label, y = -0.07, x = 0.5, size = 10,
           parse = F) +
  theme_bw(base_size = 28) +
  theme(legend.position = c(0.905, 0.155), legend.background = element_rect(color = "black"))
dev.off()

# Checking marginals vs Mehlhoff et al.'s mean and median fitness.
mehlhoff.mean.cor.label <- paste("R =", signif(mehlhoff.mean.cor[1], 2), "p =", signif(mehlhoff.mean.cor[4], 1))
ggplot(data = slopes.aa.props, aes(y = MarginalLogNBWMax, x = (1 - MeanFit))) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin=MarginalLogNBWMax - MarginalLogNBWMaxErr,
                    ymax=MarginalLogNBWMax + MarginalLogNBWMaxErr), size = 1) +
  xlab("Fitness change") +
  ylab("Effect on genotype freq / cycle") +
  geom_text(aes(label = OneLetter), hjust = -0.2, vjust = -0.4, size = 9, color = "grey30") +
  annotate("text", label = mehlhoff.mean.cor.label, y = -0.07, x = 0.055, size = 10,
          parse = F) +
  theme_bw(base_size = 28)

