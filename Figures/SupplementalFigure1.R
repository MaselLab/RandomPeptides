library(seqinr)
library(tidyverse)
library(readxl)
library(readr)
library(dplyr)
library(plotrix)

df.negbin_rev <- NegBinOut_rev

ggplot(data = df.negbin_rev) +
  geom_point(mapping = aes(y = 1/(WEIGHT), x = `FIT EST`)) +
  theme_bw(base_size = 28) +
  labs(x = 'Fitness', y = "Variance estimate")