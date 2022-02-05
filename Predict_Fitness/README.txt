All scripts for generating our fitness estimates go here.

Mathematica - Code
	-NegBinCorrection.nb: takes in data the Python Sorted data of the subpopulations growth over time "Growth_Vector _Data.m", as well as the peptide IDs "PEPIDS.tsv" and the initial estimates "NegBinOut_w _InitFreq.tsv" to generate corrected estimates for the fitness, weight, and inflation factor for the data PCR reads data (generates NegBinOut-.tsv files).
	-NegBinWeights.nb: this code is a cut down version of the previous made to only generate weights and therefore takes in many of the outputs of the previous file.
	-PoissonCorrection.nb: this is the naive model where we do not take into account PCR inflation and therefore assume a regular Poisson model of the PCR reads, not "overdispersed" or approximately Negative Binomial as above.

NegBinomial - Output
	-NegBinOutFitCorrect.tsv - this file contains the mean fitness correction (i.e. over time the average fitness is "moving" as some lineages die off and others become much more abundant)
	-NegBinOutInflation.tsv - estimated PCR inflation factor
	-NegBinOutWeights_Full.tsv - fitness estimates and weights (observed fisher information)
	-NegBinOut_w_InitFreq.tsv - the regular fitness estimates output of the NegBinCorrection.nb, but including initial frequencies of the lineages. 
	-NegBin_Weights.tsv - estimated weights for the regression model (observed fisher information)
	-Neme_Data.csv - original data from Neme et al.

Python Sorting - Code & Output
	-Growth_Vectors_Parse.py: processes the Neme_Data.csv so that it is useable for the scripts above.
R Code - Graphing
	-Neme_vs_MLE_Histogram.R: produces a histogram identifying possibly "beneficial" lineages using MLE.
