# RandomPeptides
This repository contains scripts for the paper, "Random peptides rich in small and disorder promoting amino acids are less likely to be harmful", by Luke Kosinski, Nathan Aviles, Kevin Gomez, and Joanna Masel.

The repository is organized into the several folders. For those interested in our fitness estimation technique, from which fitness can be estimated from counts while accounting for changing mean fitness, please see the "Predict fitness" folder and its sub-readme. For those interested in estimating fitness using amino acid frequencies, or for a script that can be used to summarize amino acid sequencies using various metrics, please see the "Metrics" folder.

A description of each folder is given below:

	Data: Contains the data generated from the present study. The data from Neme et al. (2017) can be found at Dryad http://dx.doi.org/10.5061/dryad.6f356, and the original sequences are available at the European Nucleotide Archive (ENA) under the project number PRJEB19640.
	Figures: All the scripts used to generate the figures found in the paper can be found here.
	Metrics: Contains two scripts, one for predicting fitness using amino acid frequencies, and another for summarizing amino acid sequence properties using a variety of metrics.
	Model: Contains scripts related to generating and testing the regression models found in the paper.
	Predict fitness: Contains all the code used for predicting fitness from sequencing reads. This folder has its own readme; please refer to that readme for more information.
