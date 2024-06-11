Censored Panel Quantile Regression with Fixed Effects via an Asymmetric Link Function 

Fulden Komuryakan & Selahattin Guris

Overview

This repository contains the R codes used for the simulation analyses in the manuscript titled "Censored Panel Quantile Regression with Fixed Effects via an Asymmetric Link Function" by Fulden Komuryakan & Selahattin Guris. The analysis includes 2- and 3-step estimators for censored quantile regression models with fixed effects. In the first step of these estimators, the informative subset is determined by estimating the censoring probability. Previous studies commonly assume symmetric link functions, which may fail to capture the best subset when the censoring probability distribution is skewed. Given that asymmetrically censored data are prevalent in empirical analyses, we propose an asymmetric link function to more accurately determine the subset.

Contents

1. Overview
2. Libraries
3. Data Generating Process
4. Simulation
4.1. Two-step estimator
4.1.1. Generalized linear models
4.1.1.1. Logit
4.1.1.2. Clog-log
4.1.2. Generalized additive models
4.1.2.1. Logit
4.1.2.2. Clog-log
4.2. Three-step estimator
4.2.1. Logit
4.2.2. Clog-log
4.3. Omniscient and Naive

Installation

To run the R scripts, you need to have R installed on your system. You can download and install R from CRAN.

Additionally, you will need the following R packages:

quantreg
MASS
splines
gam
mfx

You can install them using the following commands:

install.packages("quantreg")
install.packages("MASS")
install.packages("splines")
install.packages("gam")
install.packages("mfx")

Contact

For any questions or comments, please contact Fulden Komuryakan at fkomuryakan@bandirma.edu.tr.

