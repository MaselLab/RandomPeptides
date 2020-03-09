# Figure 3 script.

# Load packages
library(tidyverse)
library(lme4)
library(Hmisc)

# Load data.
peptide.data <- read.table(file = "Data/peptide_data_clusters_2-14-20.tsv", header = T, stringsAsFactors = F)

# Date of file generation.
today.date <- "3-9-2020"

# Calculating the marginal effects.
# Amino acid frequency only model.
peptide.mixed.nb.log.aaonly.lm <-
  lmer(
    data = peptide.data,
    formula = log2(Fitness.nb) ~
      Leu + Pro + Met + Trp + Ala +
      Val + Phe + Ile + Gly + Ser +
      Thr + Cys + Asn + Gln + Tyr +
      His + Asp + Glu + Lys + Arg +
      (1|Cluster) +
      0,
    weights = Weight.nb
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
            Leu.freq = wtd.mean(Leu.freq, weights = Weight.nb), Phe.freq = wtd.mean(Phe.freq, weights = Weight.nb),
            Met.freq = wtd.mean(Met.freq, weights = Weight.nb), Val.freq = wtd.mean(Val.freq, weights = Weight.nb),
            Ile.freq = wtd.mean(Ile.freq, weights = Weight.nb), Lys.freq = wtd.mean(Lys.freq, weights = Weight.nb),
            His.freq = wtd.mean(His.freq, weights = Weight.nb), Arg.freq = wtd.mean(Arg.freq, weights = Weight.nb),
            Glu.freq = wtd.mean(Glu.freq, weights = Weight.nb), Asp.freq = wtd.mean(Asp.freq, weights = Weight.nb),
            Gln.freq = wtd.mean(Gln.freq, weights = Weight.nb), Asn.freq = wtd.mean(Asn.freq, weights = Weight.nb),
            Gly.freq = wtd.mean(Gly.freq, weights = Weight.nb), Ala.freq = wtd.mean(Ala.freq, weights = Weight.nb),
            Pro.freq = wtd.mean(Pro.freq, weights = Weight.nb), Ser.freq = wtd.mean(Ser.freq, weights = Weight.nb),
            Trp.freq = wtd.mean(Trp.freq, weights = Weight.nb), Tyr.freq = wtd.mean(Tyr.freq, weights = Weight.nb),
            Thr.freq = wtd.mean(Thr.freq, weights = Weight.nb), Cys.freq = wtd.mean(Cys.freq, weights = Weight.nb),
            Clustering.Six = wtd.mean(Clustering.Six, weights = Weight.nb),
            WaltzBinary = wtd.mean(WaltzBinary, weights = Weight.nb),
            WaltzAAsInAPRs = wtd.mean(WaltzAAsInAPRs, weights = Weight.nb),
            Waltz.delta = wtd.mean(Waltz.delta, weights = Weight.nb),
            TangoBinary = wtd.mean(TangoBinary, weights = Weight.nb),
            TangoAAsInAPRs = wtd.mean(TangoAAsInAPRs, weights = Weight.nb),
            Tango.delta = wtd.mean(Tango.delta, weights = Weight.nb),
            Cys.squared = wtd.mean(Cys.squared, weights = Weight.nb),
            CamSol.avg = wtd.mean(CamSol.avg, weights = Weight.nb),
            AnchorAvg = wtd.mean(AnchorAvg, weights = Weight.nb),
            ISD.delta = wtd.mean(ISD.delta, weights = Weight.nb),
            mean.run.norm = wtd.mean(mean.run.norm, weights = Weight.nb),
            max.run.length = wtd.mean(max.run.length, weights = Weight.nb),
            net.charge = wtd.mean(net.charge, weights = Weight.nb),
            pI = wtd.mean(pI, weights = Weight.nb), pI.ecoli = wtd.mean(pI.ecoli, weights = Weight.nb),
            charge.pos = wtd.mean(charge.pos, weights = Weight.nb), charge.neg = wtd.mean(charge.neg)
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
aa.freqs.wmax[match(names(log.aaonly.summary$coefficients[,1]), aa.freqs.wmax[,"AminoAcid"]),]
aa.freqs.ordered <- aa.sort.lm(peptide.mixed.nb.log.aaonly.lm, aa.freqs.wmax)
aa.freqs.ordered
as.list(aa.freqs.ordered[,"Count"])[[1]]

marginals.clusters.log.nb.wmax <- custom.contrasts(aa.freqs.wmax, peptide.mixed.nb.log.aaonly.lm)
marginals.clusters.log.nb.wmax

# Exporting marginal effects as a table.
write_tsv(marginals.clusters.log.nb.wmax, path = "Data/supplemental_table_2.tsv")

# Importing data for James et al.'s phylostratigraphy slopes.
slopes.aa.props <- read_tsv(file = "Data/phylostratigraphy_aasummary_3-4-2020.tsv")
slopes.aa.props

# Renaming the columns of the marginals data frame to match those from James et al.'s data frame.
names(marginals.clusters.log.nb.wmax)
colname.marginals <- "MarginalLogNBWMax"
colname.stderr <- "MarginalLogNBWMaxErr"
names(marginals.clusters.log.nb.wmax) <- c("AA", colname.marginals, colname.stderr)
names(marginals.clusters.log.nb.wmax)
marginals.clusters.log.nb.wmax

# Deleting any columns in the James et al. data frame with these names.
slopes.aa.props$MarginalLogNBWMax <- NULL
slopes.aa.props$MarginalLogNBWMaxErr <- NULL

# Merging the James et al. data frame and the marginal effects data frames.
slopes.aa.props <- merge(slopes.aa.props, marginals.clusters.log.nb.wmax, by = "AA")

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
  y = slopes.aa.props$AnimalNonTransmembraneYoungPfamSlopes,
  y.var = slopes.aa.props$AnimalNonTransmembraneYoungPfamError
)
animalpfam.cor
plantpfam.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$PlantNonTransmembraneYoungPfamSlopes,
  y.var = slopes.aa.props$PlantNonTransmembraneYoungPfamError
)
plantpfam.cor
ancientpfam.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$AllNonTransmembraneOldPfamSlopes,
  y.var = slopes.aa.props$AllNonTransmembraneOldPfamError
)
ancientpfam.cor

# And now checking phylostratigraphy slopes for full genes.
##################################################################################
# NOTE: THESE HAVE ERRORS! FOR SOME REASON THE PLANT SLOPES ARE THE SAME AS ANIMAL
# NEED TO GET THE UPDATED SLOPES FROM JENNY
##################################################################################
animalgene.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$`Animal_<1496_Slope`,
  y.var = slopes.aa.props$`Animal_<1496_Error`
)
animalgene.cor
plantgene.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$`Plant_<1496_Slope`,
  y.var = slopes.aa.props$`Plant_<1496_Error`
)
plantgene.cor
ancientgene.cor <- custom.weighted.pearson(
  x = slopes.aa.props$MarginalLogNBWMax,
  x.var = slopes.aa.props$MarginalLogNBWMaxErr,
  y = slopes.aa.props$`All_>2101_Slope`,
  y.var = slopes.aa.props$`All_>2101_Error`
)
ancientgene.cor

# Making the figures.
partB.name <- paste("Scripts/Figures/NonTransmembranePlantYoungPfam_MarginalLogNB_notrel_", today.date, ".png", sep = "")
plant.cor.label <- paste("R =", signif(plantpfam.cor[1], 2), "p =", signif(plantpfam.cor[4], 1))
png(filename = partB.name, width = 600, height = 600)
ggplot(data = slopes.aa.props, aes(x = MarginalLogNBWMax, y = PlantNonTransmembraneYoungPfamSlopes)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin=PlantNonTransmembraneYoungPfamSlopes - PlantNonTransmembraneYoungPfamError,
                    ymax=PlantNonTransmembraneYoungPfamSlopes + PlantNonTransmembraneYoungPfamError),
                size = 1) +
  geom_errorbarh(aes(xmin=MarginalLogNBWMax - MarginalLogNBWMaxErr,
                     xmax=MarginalLogNBWMax + MarginalLogNBWMaxErr), size = 1) +
  ylab("\U0394% points change per BY") +
  xlab("Multiplier effect on genotype freq / cycle\n(fold change units)") +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1),
                     labels = round(2^c(-0.2, -0.1, 0, 0.1), digits = 2)) +
  scale_y_continuous(breaks = c(1e-5, 0, -1e-5, -2e-5, -3e-5),
                     limits = c(-3e-5, 1.8e-5),
                     labels = c("1%", "0%", "-1%", "-2%", "-3%")) +
  geom_text(aes(label = OneLetter), hjust = -0.2, vjust = -0.4, size = 9, color = "grey30") +
  annotate("text", label = plant.cor.label, x = -0.12, y = -1.5e-5, size = 10,
           parse = F) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

partA.name <- paste("Scripts/Figures/NonTransmembraneAnimalYoungPfam_MarginalLogNB_notrel_", today.date, ".png", sep = "")
animal.cor.label <- paste("R =", signif(animalpfam.cor[1], 2), "p =", signif(animalpfam.cor[4], 1))
png(filename = partA.name, width = 600, height = 600)
ggplot(data = slopes.aa.props, aes(x = MarginalLogNBWMax, y = AnimalNonTransmembraneYoungPfamSlopes)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin=AnimalNonTransmembraneYoungPfamSlopes - AnimalNonTransmembraneYoungPfamError,
                    ymax=AnimalNonTransmembraneYoungPfamSlopes + AnimalNonTransmembraneYoungPfamError),
                size = 1) +
  geom_errorbarh(aes(xmin=MarginalLogNBWMax - MarginalLogNBWMaxErr,
                     xmax=MarginalLogNBWMax + MarginalLogNBWMaxErr), size = 1) +
  ylab("\U0394% points change per BY") +
  xlab("Multiplier effect on genotype freq / cycle\n(fold change units)") +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1),
                     labels = round(2^c(-0.2, -0.1, 0, 0.1), digits = 2)) +
  scale_y_continuous(breaks = c(1e-5, 0, -1e-5, -2e-5, -3e-5),
                     limits = c(-3e-5, 1.8e-5),
                     labels = c("1%", "0%", "-1%", "-2%", "-3%")) +
  geom_text(aes(label = OneLetter), hjust = -0.2, vjust = -0.4, size = 9, color = "grey30") +
  annotate("text", label = animal.cor.label, x = -0.1, y = -1.5e-5, size = 10,
           parse = F) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

partC.name <- paste("Scripts/Figures/NonTransmembraneAncientPfam_MarginalLogNB_notrel_", today.date, ".png", sep = "")
ancient.cor.label <- paste("R =", signif(ancientpfam.cor[1], 2), "p =", signif(ancientpfam.cor[4], 1))
png(filename = partC.name, width = 600, height = 600)
ggplot(data = slopes.aa.props, aes(x = MarginalLogNBWMax, y = AllNonTransmembraneOldPfamSlopes)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin=AllNonTransmembraneOldPfamSlopes - AllNonTransmembraneOldPfamError,
                    ymax=AllNonTransmembraneOldPfamSlopes + AllNonTransmembraneOldPfamError),
                size = 1) +
  geom_errorbarh(aes(xmin=MarginalLogNBWMax - MarginalLogNBWMaxErr,
                     xmax=MarginalLogNBWMax + MarginalLogNBWMaxErr), size = 1) +
  ylab("\U0394% points change per BY") +
  xlab("Multiplier effect on genotype freq / cycle\n(fold change units)") +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1),
                     labels = round(2^c(-0.2, -0.1, 0, 0.1), digits = 2)) +
  scale_y_continuous(breaks = c(1e-5, 0, -1e-5, -2e-5, -3e-5),
                     limits = c(-3e-5, 1.8e-5),
                     labels = c("1%", "0%", "-1%", "-2%", "-3%")) +
  geom_text(aes(label = OneLetter), hjust = -0.2, vjust = -0.4, size = 9, color = "grey30") +
  annotate("text", label = ancient.cor.label, x = -0.1, y = -1.5e-5, size = 10,
           parse = F) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()
