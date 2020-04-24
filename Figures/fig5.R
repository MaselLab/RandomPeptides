# Looking at the fitness of random linker regions from Frumkin et al. 2017.

# Load packages.
library(tidyverse)
library(MASS)
library(ggpubr)

# Box-Cox transformation.
bc.transform <- function(data.vector, lambda){
  data.transformed <- ((data.vector ^ lambda) - 1) / lambda
  return(data.transformed)
}

bc.back <- function(data.transformed, lambda){
  data.orginal <- ((data.transformed * lambda) + 1) ^ (1/lambda)
  return(data.orginal)
}

# Global variables for quick editing.
todays.date <- "4-6-2020"
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Load Frumkin et al. data, kindly sent to us by Dvir Schirman.
frumkin.data <- read_csv(file = "Data/peptide_fitness.csv")
frumkin.data

# Load fitness predictor.
source(file = "Scripts/predict_fitness.R")

# Predicting fitness from AA sequence.
frumkin.data$Predicted.fitness <- fitness.calculator(aa.sequence = frumkin.data$AA_seq)

# Looking at the relationship between fitness and predicted fitness.
ggplot(
  data = frumkin.data,
  aes(
    x = Predicted.fitness,
    y = fitness_residuals_mean
  )
) +
  geom_point() +
  geom_smooth()

# Looking at the distributions.
quantile(frumkin.data$fitness_residuals_mean)
hist(frumkin.data$fitness_residuals_mean)
fit.resid.min.1 <- 1 - min(frumkin.data$fitness_residuals_mean)
hist(frumkin.data$fitness_residuals_mean + fit.resid.min.1)
boxcox(frumkin.data$fitness_residuals_mean + fit.resid.min.1 ~ 1, lambda = seq(100, 400, by = 0.1))
fit.resid.lambda = 250
hist(bc.transform(frumkin.data$fitness_residuals_mean + fit.resid.min.1, lambda = fit.resid.lambda))

hist(frumkin.data$Predicted.fitness)
boxcox(frumkin.data$Predicted.fitness ~ 1)
# Lambda = 0
hist(log(frumkin.data$Predicted.fitness))

ggplot(
  data = frumkin.data,
  aes(
    x = Predicted.fitness,
    y = fitness_residuals_mean
  )
) +
  geom_point() +
  geom_smooth(method = "loess", color = "red") +
  geom_smooth(method = "lm") +
  ylab("Fitness residuals") +
  xlab("Predicted fitness") +
  theme_bw(base_size = 28)

# Making a better histogram of transformed mean fitness residuals.
ggplot(
  data = frumkin.data,
  aes(
    x = bc.transform(frumkin.data$fitness_residuals_mean + fit.resid.min.1, lambda = fit.resid.lambda)
  )
) +
  geom_histogram(bins = 11) +
  xlab("Fitness residuals") +
  scale_x_continuous(breaks = bc.transform(c(-0.01, 0, 0.002) + fit.resid.min.1, lambda = 250),
                     labels = c(-0.01, 0, 0.002)) +
  theme_bw(base_size = 28)

# Spearman's correlation.
linker.spearman <- cor.test(frumkin.data$fitness_residuals_mean, frumkin.data$Predicted.fitness, method = "spearman")
linker.spearman
linker.cor.label <- paste("Spearman's rho = ", round(linker.spearman$estimate, 2), ", ", 
                          "p = ", round(linker.spearman$p.value, 3), sep = "")
linker.cor.label

# Checking the relationship with predicted fitness.
path.title.plot <- paste("Scripts/Figures/linker_fitness_predicted_", todays.date, ".png", sep = "")
png(filename = path.title.plot, width = 500, height = 500)
ggplot(
  data = frumkin.data,
  aes(
    x = Predicted.fitness,
    y = fitness_residuals_mean
  )
) +
  geom_point() +
  #geom_smooth(method = "loess", color = cbbPalette[2]) +
  #geom_smooth(method = "lm", color = cbbPalette[6]) +
  xlab("Predicted fitness") +
  ylab("Fitness residual") +
  #stat_cor(method = "spearman", label.x = 0.1, label.y = -0.007, size = 5) +
  annotate("text", label = linker.cor.label, x = 0.35, y = -0.007, size = 8,
           parse = F) +
  # scale_y_continuous(breaks = bc.transform(c(-0.01, 0, 0.002) + fit.resid.min.1, lambda = 250),
  #                    labels = c(-0.01, 0, 0.002)) +
  # scale_x_continuous(breaks = log(c(0.05, 0.15, 0.4)),
  #                    labels = c(0.05, 0.15, 0.4)) +
  theme_bw(base_size = 28)
dev.off()

# Building a simple linear regression model.
linker.lm <- lm(
  data = frumkin.data,
  formula = frumkin.data$fitness_residuals_mean ~ Predicted.fitness
)
summary(linker.lm)
#plot(linker.lm)
drop1(linker.lm, test = "Chisq")
