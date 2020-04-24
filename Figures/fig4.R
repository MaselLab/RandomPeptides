# Weighted T-test of GC content 1st or 2nd position of amino acid vs marginal effects.

# Load libraries.
library(tidyverse)
library(weights)

# Global variables.
todays.date <- "4-23-2020"
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Load marginal effects.
marginals <- read_tsv("Scripts/RandomPeptides/Data/supplemental_table_2.tsv")
marginals

# Adding in weights.
marginals$Weight <- 1/(marginals$Std.Err ^ 2)
marginals

# Adding first and second position GC or AT preference.
gc.first.second <- tibble(
  "AminoAcid" = marginals$AminoAcid,
  "GC.first" = c(NA, "GC", "AT", "AT", "GC",
                 "GC", "AT", "AT", "GC", "AT",
                 "AT", "AT", "AT", "GC", "AT",
                 "GC", "GC", "GC", "AT", NA),
  "GC.second" = c("AT", "GC", "AT", "GC", "GC",
                  "AT", "AT", "AT", "GC", "GC",
                  "GC", "GC", "AT", "AT", "AT",
                  "AT", "AT", "AT", "AT", "GC")
)
gc.first.second

marginals.gc <- merge(marginals, gc.first.second, by = "AminoAcid")
marginals.gc

# Looking at the standard errors between the groups.
summary(marginals.gc[marginals.gc$GC.first == "AT",]$Std.Err)
summary(marginals.gc[marginals.gc$GC.first == "GC",]$Std.Err)
summary(marginals.gc[marginals.gc$GC.second == "AT",]$Std.Err)
summary(marginals.gc[marginals.gc$GC.second == "GC",]$Std.Err)

# T-tests.
first.pos.ttest <- wtd.t.test(y = marginals.gc[marginals.gc$GC.first == "AT",]$MarginalEffect,
                              x = marginals.gc[marginals.gc$GC.first == "GC",]$MarginalEffect,
                              weighty = marginals.gc[marginals.gc$GC.first == "AT",]$Weight,
                              weight = marginals.gc[marginals.gc$GC.first == "GC",]$Weight)
first.pos.ttest
2 ^ first.pos.ttest$additional[[1]]

second.pos.ttest <- wtd.t.test(y = marginals.gc[marginals.gc$GC.second == "AT",]$MarginalEffect,
                               x = marginals.gc[marginals.gc$GC.second == "GC",]$MarginalEffect,
                               weighty = marginals.gc[marginals.gc$GC.second == "AT",]$Weight,
                               weight = marginals.gc[marginals.gc$GC.second == "GC",]$Weight)
second.pos.ttest
2 ^ second.pos.ttest$additional[[1]]

# Checking confidence intervals.
t.stat.first <- qt(0.975, df = first.pos.ttest$coefficients[[2]])
t.stat.first
2 ^ (first.pos.ttest$additional[[1]] + t.stat.first * first.pos.ttest$additional[[4]]) # 1.11
2 ^ (first.pos.ttest$additional[[1]] - t.stat.first * first.pos.ttest$additional[[4]]) # 0.98

t.stat.second <- qt(0.975, df = second.pos.ttest$coefficients[[2]])
t.stat.second
2 ^ (second.pos.ttest$additional[[1]] + t.stat.second * second.pos.ttest$additional[[4]]) # 1.13
2 ^ (second.pos.ttest$additional[[1]] - t.stat.second * second.pos.ttest$additional[[4]]) # 1.03

# Adding one letter amino acid abbreviations for graphing.
marginals.gc$OneLetter <-
  c("A", "R", "N", "D", "C",
    "Q", "E", "G", "H", "I",
    "L", "K", "M", "F", "P",
    "S", "T", "W", "Y", "V")
marginals.gc

# Plotting these.
codon.first.name <- paste("Scripts/Figures/GCvsAT_first_marginals_", todays.date, ".png", sep = "")
png(filename = codon.first.name, width = 600, height = 600)
ggplot(
  data = marginals.gc,
  aes(
    x = GC.first,
    y = MarginalEffect,
    size = Weight
  )
) +
  geom_point(alpha = 0.4) +
  geom_text(aes(label = OneLetter), hjust = -0.4, vjust = -0.01, size = 9, color = "grey30") +
  theme_bw(base_size = 28) +
  xlab("GC vs AT,\n1st nucleotide in codon") +
  ylab("Effect on genotype freq / cycle\n(fold change)") +
  scale_y_continuous(breaks = c(-0.2, -0.1, 0, 0.1),
                     labels = round(2^c(-0.2, -0.1, 0, 0.1), digits = 2)) +
  theme(legend.position = "none")
dev.off()

codon.second.name <- paste("Scripts/Figures/GCvsAT_second_marginals_", todays.date, ".png", sep = "")
png(filename = codon.second.name, width = 600, height = 600)
ggplot(
  data = marginals.gc,
  aes(
    x = GC.second,
    y = MarginalEffect,
    size = Weight
  )
) +
  geom_point(alpha = 0.4) +
  geom_text(aes(label = OneLetter), hjust = -0.4, vjust = -0.01, size = 9, color = "grey30") +
  theme_bw(base_size = 28) +
  xlab("GC vs AT,\n2nd nucleotide in codon") +
  ylab("Effect on genotype freq / cycle\n(fold change)") +
  scale_y_continuous(breaks = c(-0.2, -0.1, 0, 0.1),
                     labels = round(2^c(-0.2, -0.1, 0, 0.1), digits = 2)) +
  theme(legend.position = "none")
dev.off()

# Plotting the confidence intervals.
first.second.ttest.df <- tibble(
  "pos" = factor(c("first", "second")),
  "diff" = c(first.pos.ttest$additional[[1]], second.pos.ttest$additional[[1]]),
  "se" = c(first.pos.ttest$additional[[4]], second.pos.ttest$additional[[4]])
)
first.second.ttest.df

first.second.comparison <- paste("Scripts/Figures/GCvsAT_mean_comparison_", todays.date, ".png", sep = "")
png(filename = first.second.comparison, width = 600, height = 600)
ggplot(
  data = first.second.ttest.df,
  aes(
    x = pos,
    y = diff
  )
) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = diff - se, ymax = diff + se), width = 0.4, size = 2) +
  ylab("Effect of GC - effect of AT\non genotype freq / cycle (fold change)") +
  scale_y_continuous(breaks = c(0.04, 0.08, 0.12, 0.16),
                     labels = round(2^c(0.04, 0.08, 0.12, 0.16), digits = 2)) +
  xlab("Codon position") +
  scale_x_discrete(labels = c("1st", "2nd")) +
  theme_bw(base_size = 28)
dev.off()

