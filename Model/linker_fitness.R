# Looking at the fitness of random linker regions from Frumkin et al. 2017.

# Load packages.
library(tidyverse)
library(MASS)

# Global variables for quick editing.
todays.date <- "3-16-2020"
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
hist(frumkin.data$fitness_residuals_mean)
fit.resid.min.1 <- 1 - min(frumkin.data$fitness_residuals_mean)
hist(frumkin.data$fitness_residuals_mean + fit.resid.min.1)
boxcox(frumkin.data$fitness_residuals_mean + fit.resid.min.1 ~ 1, lambda = seq(100, 400, by = 0.1))
hist((frumkin.data$fitness_residuals_mean + fit.resid.min.1)^250)

hist(frumkin.data$Predicted.fitness)
boxcox(frumkin.data$Predicted.fitness ~ 1)
hist(log(frumkin.data$Predicted.fitness))

# Checking the relationship with predicted fitness.
path.title.plot <- paste("Scripts/Figures/linker_fitness_predicted_", todays.date, ".png", sep = "")
png(filename = path.title.plot, width = 500, height = 500)
ggplot(
  data = frumkin.data,
  aes(
    x = log(Predicted.fitness),
    y = (frumkin.data$fitness_residuals_mean + fit.resid.min.1)^250
  )
) +
  geom_point() +
  geom_smooth(method = "loess", color = cbbPalette[2]) +
  geom_smooth(method = "lm", color = cbbPalette[6]) +
  xlab("Predicted fitness") +
  ylab("Fitness residual") +
  scale_y_continuous(breaks = (c(-0.005, 0, 0.002) + fit.resid.min.1)^250,
                     labels = c(-0.005, 0, 0.002)) +
  scale_x_continuous(breaks = log(c(0.05, 0.15, 0.4)),
                     labels = c(0.05, 0.15, 0.4)) +
  theme_bw(base_size = 28)
dev.off()

# Building a simple linear regression model.
linker.lm <- lm(
  data = frumkin.data,
  formula = I((frumkin.data$fitness_residuals_mean + fit.resid.min.1)^250) ~ log(Predicted.fitness)
)
summary(linker.lm)
plot(linker.lm)
drop1(linker.lm, test = "Chisq")
