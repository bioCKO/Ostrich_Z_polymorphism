#**********************************************************************************************#
#**********************************************************************************************#
# Statistical analysis using the black subspecies**********************************************#
# Homa Papoli**********************************************************************************#
#**********************************************************************************************#
rm(list=ls())
library(psych)
#***********************************************************************************************
# Read data into R
#***********************************************************************************************
black1Mb_stat <- read.table("../data/black.1000Kb.INTERGENIC.rec.resampling.Z.forR.txt", header = T)
head(black1Mb_stat)
#blue1Mb_stat <- read.table("../data/blue.1000Kb.INTERGENIC.rec.resampling.Z.forR.txt", header = T)
#red1Mb_stat <- read.table("../data/red.1000Kb.INTERGENIC.rec.resampling.Z.forR.txt", header = T)
#***********************************************************************************************
# Separating dataset to PAR and non-PAR
blackPAR1Mb <- black1Mb_stat[black1Mb_stat$ChrZ_start < 53065700,]
blacknonPAR1Mb <- black1Mb_stat[black1Mb_stat$ChrZ_start > 53065700,]

#bluePAR1Mb <- blue1Mb_stat[blue1Mb_stat$ChrZ_start < 53065700,]
#bluenonPAR1Mb <- blue1Mb_stat[blue1Mb_stat$ChrZ_start > 53065700,]

#redPAR1Mb <- red1Mb_stat[red1Mb_stat$ChrZ_start < 53065700,]
#rednonPAR1Mb <- red1Mb_stat[red1Mb_stat$ChrZ_start > 53065700,]
#***********************************************************************************************
# Exploratory data analysis
#***********************************************************************************************
# Plot pairwise correlations between all variables (response and explanatory)
#***********************************************************************************************
library("PerformanceAnalytics")
# Z
black1000Kb_matrix <- black1Mb_stat[c(7,11,12,13,16,17,18,19)]
chart.Correlation(black1000Kb_matrix, histogram=TRUE, pch=19)
# PAR
black1000Kb_PAR <- blackPAR1Mb[c(7,11,12,13,16,17,18,19)]
chart.Correlation(black1000Kb_PAR, histogram=TRUE, pch=19)
# non-PAR
black1000Kb_nonPAR <- blacknonPAR1Mb[c(7,11,12,13,16,17,18,19)]
chart.Correlation(black1000Kb_nonPAR, histogram=TRUE, pch=19)
#***********************************************************************************************
# Check for outlier
#***********************************************************************************************
par(mfrow=c(2,4))  # divide graph area in 2 columns
boxplot(black1Mb_stat$pi_per_window, ylab = "Diversity (pi)")
boxplot(black1Mb_stat$CDS_density, ylab = "CDS density")
boxplot(black1Mb_stat$GC_density, ylab = "GC density")
boxplot(black1Mb_stat$repeat_density, ylab = "Repeat density")
boxplot(black1Mb_stat$CM_per_bp, ylab = "CM per bp")
boxplot(black1Mb_stat$male_CM_per_bp, ylab = "male CM per bp")
boxplot(black1Mb_stat$female_CM_per_bp, ylab = "female CM per bp")
boxplot(black1Mb_stat$CM_from_PAR_boundary, ylab = "CM from PAR boundary")
#***********************************************************************************************
# Conclusion of exploratory data analysis:
# Z chromosomes
# It seems that one explanatory variable, repeat density, can be safely removed. It shows no 
# correlation with diversity or any other of the explanatory variables.
# There is signiciant correlation between diversity (pi) and 2 explanatory variables:
# 1. CM_from_PAR_boundary
# 2. GC_density
#***********************************************************************************************

#***********************************************************************************************
# 1. CM_from_PAR_boundary

# CM_from_PAR_boundary is correlated with GC_density. Biological explanation for this correlation
# is difficult unless areas furthest away from the PAR boundary also experience highest rate of
# recombiantion. This is plausible given that areas furthest away from the PAR boundary are closer
# to the chromosomes' end which experience higher rates of recombination. We can check this using
# coplot of GC_density as a function of CM_from_PAR_boundary in different levels of recombination
# rate (CM_per_bp).

# For PAR

given.CM_per_bp <- co.intervals(black1000Kb_PAR$CM_per_bp, number = 3)
coplot(GC_density ~ CM_from_PAR_boundary | CM_per_bp, data = black1000Kb_PAR, 
       given.v = given.CM_per_bp, rows = 1)
# For PAR it is not quite clear that higher GC values further away from the PAR boundary coincide
# with higher recombination rate.

# For nonPAR

given.CM_per_bp <- co.intervals(black1000Kb_nonPAR$CM_per_bp, number = 3)
coplot(GC_density ~ CM_from_PAR_boundary | CM_per_bp, data = black1000Kb_nonPAR, 
       given.v = given.CM_per_bp, rows = 1)
# For nonPAR, we see higher values of GC further away from PAR boundary coinciding with higher
# rate of recombination. The positive slope between GC density and cM from PAR boundary appears
# only in regions with higher recombination.

#***********************************************************************************************
# 2. GC density

# The strong correlation between diversity (pi) and GC density can be explained biologically 
# through an increased mutation rate in GC sites. As neutral genetic diversity is a function of
# mutation rate and Ne, regions with higher mutation rate such as GC dense regions may show an 
# increase in genetic diversity. 
# GC density is signifcantly correlated with recombination rate and genetic distance from the
# PAR boundary. GC density is also correlated with CDS density. The level of 
# GC content in avian coding sequence, as shown in various neognathae species, is higher than in
# nocoding sequences so in a given window, it is expected to see higher density of GC with 
# higher density of coding sequence.

# Non-significant but important explanatory variables
#***********************************************************************************************
# Recombination rate
#***********************************************************************************************
# The correlation matrix shows a significant correlation between CM_per_bp, male_CM_per_bp and 
# female_CM_per_bp. This is expected given that CM_per_bp is just the sex average recombination
# rate. For the model, we use the sex averaged recombination rate (CM_per_bp) given that genetic 
# diversity on the Z is influenced by the average recombination rate in both males and females 
# (non-PAR Z spends 1/3 of the time in females and 2/3 of the time in male; PAR Z spends equal
# amount of time in males and females.).

# CM_per_bp is correlated with two other explanatory variables: GC_density and CDS_density. Both
# of these correlations are expected biologically. The correlation between CM_per_pb and
# GC_density is weak (r = 0.2). It has been shown that DNA repair at recombination double-strand
# breaks is biased towards GC. That is if we have a double strand break where we have heterozygous
# in the form A/G or T/G, A/C or T/C, the resulting repaired sequence will likely be GG, CC, GC
# or CG. 

# The correlation between CM_per_bp and CDS_density is significant. This is also biologically
# meaningful. In birds, it has been shown that recombination hotspots occur close to genes. It 
# is therefore expected to see higher recombination rates in places with higher CDS density.

#***********************************************************************************************
# Considering the above, we use the following explanatory variables for Z, PAR and non-PAR:
# Response variable: Diversity (pi)

# Explanatory variables: 
# 1. cM from PAR boundary
# Correlates:
# GC density 

# 2. GC density
# Correlates:
# Recombination rate
# cM from PAR boundary
# CDS density

# 3. Recombination rate
# Correlates:
# GC density
# CDS density

# 4. CDS density
# Correlates:
# GC density
# Recombination rate
#***********************************************************************************************
# Linear regression
#***********************************************************************************************
library(nlme)
# General linear model with all explanatory variables
scaled_pi_Z = scale(black1Mb_stat$Resampling_mean)
################################################################################################
full_Z <- lm(scaled_pi_Z ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) + 
               scale(CDS_density), na.action="na.exclude", data = black1Mb_stat)
summary(full_Z)
################################################################################################
scaled_pi_PAR = scale(black1000Kb_PAR$Resampling_mean)
PAR_Z <- lm(scaled_pi_PAR ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) + 
              scale(CDS_density), na.action="na.exclude", data = black1000Kb_PAR)
summary(PAR_Z)

# Add CM_from_PAR_boundary as quadratic term
CMPAR2 <- black1000Kb_PAR$CM_from_PAR_boundary^2
PAR_Z_quadratic <- lm(scaled_pi_PAR ~ scale(CM_from_PAR_boundary) + CMPAR2 + scale(GC_density) + 
                      scale(CM_per_bp) + scale(CDS_density), na.action="na.exclude", 
                      data = black1000Kb_PAR)
summary(PAR_Z_quadratic)
# Model with only CM_from_PAR_boundary
PAR_Z_quadratic_1 <- lm(scaled_pi_PAR ~ scale(CM_from_PAR_boundary) + CMPAR2, na.action="na.exclude", 
                        data = black1000Kb_PAR)
summary(PAR_Z_quadratic_1)
################################################################################################
scaled_pi_nonPAR = scale(black1000Kb_nonPAR$Resampling_mean)
nonPAR_Z <- lm(scaled_pi_nonPAR ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) +
                 scale(CDS_density), na.action="na.exclude", data = black1000Kb_nonPAR)
summary(nonPAR_Z)
#***********************************************************************************************
# Check assumptions of GLM for the Z chromosome
#***********************************************************************************************
# 1. Observations are independent of another
#***********************************************************************************************
# We run the model with na.action=na.omit again because acf() does not work with 
# na.action="na.exclude".
par(mfrow=c(1,1))
full_Z <- lm(scaled_pi_Z ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) + 
            scale(CDS_density), na.action=na.omit, data = black1Mb_stat)
E <- residuals(full_Z)
acf(E) 
# Residuals are not independent but there is auto-coorelation.
# Heterozygosity values along the chromosome are not independent of each other due to a feature 
# called Identity Disequilibrium (ID), that is the correlation in heterozygosity/homozygosity 
# among loci. This correlation comes about due to some form of consanguineous mating. 

#***********************************************************************************************
# 2. The response variable is continuous
#***********************************************************************************************
# The response variable is heterozygosity. Heterozygosity should have a binomial distribution, 
# that is the data are the number of successes and the number of failures. However, we can 
# use (average number of pairwise differences)/(average number of pairwise comparisons) as a 
# measure of heterozygosity.
#***********************************************************************************************
# 3. Error has a normal distribution
#***********************************************************************************************
ks.test(residuals(full_Z), "pnorm") 
#***********************************************************************************************
# 4. Homogeniety of variance
#***********************************************************************************************
par(mfrow=c(2,2))
plot(full_Z)

# Residuals vs fitted AND Scale-Location

# Residuals and standardized residuals as a function of fitted values show more or 
# less a straight line which may argue for linearity and homogeneity of variance.
# However, for the Scaled Location graph, there is a positive slope in residuals between -0.5
# and 0 of the Fitted values!

# Normal Q-Q
# All the points fall approximately along the reference line so we can accept normality.

# Residuals vs Leverage
# Standardized residuals can be interpreted as the number of standard errors away from 
# the regression line. Two points, 48 and 51 are reported as outliers, however, no 
# outlier exceeds 3 standard deviations.

#***********************************************************************************************
# GLS: Generalized least squares
#***********************************************************************************************
# Due to lack of independence among residuals, we use GLS model.
#***********************************************************************************************
library(nlme)
# I used number of rows (obs) in the corAR1 formula. 
black1Mb_stat$obs = as.numeric(rownames(black1Mb_stat))
head(black1Mb_stat)
# This is because when I used the physical position, Phi1 taken from summary(model) always returned 0. 

Z0 <- gls(scaled_pi_Z ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) +
            scale(CDS_density), na.action = na.omit, data = black1Mb_stat)

Z1 <- gls(scaled_pi_Z ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) +
            scale(CDS_density), na.action = na.omit, data = black1Mb_stat, 
          correlation = corAR1(form =~ obs))
Z2 <- gls(scaled_pi_Z ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) +
            scale(CDS_density), na.action = na.omit, data = black1Mb_stat, 
          correlation = corSpher(form =~ obs))
Z3 <- gls(scaled_pi_Z ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) +
            scale(CDS_density), na.action = na.omit, data = black1Mb_stat, 
          correlation = corLin(form =~ obs))
Z4 <- gls(scaled_pi_Z ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) +
            scale(CDS_density), na.action = na.omit, data = black1Mb_stat, 
          correlation = corRatio(form =~ obs))
Z5 <- gls(scaled_pi_Z ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) +
            scale(CDS_density), na.action = na.omit, data = black1Mb_stat, 
          correlation = corGaus(form =~ obs))
Z6 <- gls(scaled_pi_Z ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) +
            scale(CDS_density), na.action = na.omit, data = black1Mb_stat, 
          correlation = corExp(form =~ obs))

AIC(Z0, Z1, Z2, Z3, Z4, Z5, Z6)

# All models with spatial structure have a lower AIC, showing that adding a spatial correlation
# structure improves the model. This is particularly the case for models with corAR1, corExp and
# corSpher and corRatio.

anova(Z0, Z1)
# The value of adding a spatial structure is also confirmed by anova.

# Plot residuals against fitted values
par(mfrow=c(1,3))
plot(full_Z$fitted,full_Z$residuals)
abline(h=0)
title(main="GLM")
plot(Z0$fitted,Z0$residuals)
abline(h=0)
title(main="GLS no correlation")
plot(Z1$fitted,Z1$residuals)
abline(h=0)
title(main="GLS correlation")

# Plot ACF graph for residuals
E_lm <- residuals(full_Z)
acf(E_lm, ylab="ACF glm")
E_gls <- residuals(Z0)
acf(E_gls, ylab = "gls no corr.") 
E_gls_cor <- residuals(Z1)
acf(E_gls_cor, ylab = "gls corr.") 

# Dropping less significant values one by one for all Z, PAR and non-PAR:

Z1 <- gls(scaled_pi_Z ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) +
            scale(CDS_density), na.action = "na.exclude", data = black1Mb_stat, 
          correlation = corAR1(form =~ obs))
summary(Z1)
Z1A <- gls(scaled_pi_Z ~ scale(GC_density) + scale(CM_per_bp) + scale(CDS_density), 
           na.action = "na.exclude", data = black1Mb_stat, 
           correlation = corAR1(form =~ obs))
summary(Z1A)
Z2A <- gls(scaled_pi_Z ~ scale(GC_density) + scale(CDS_density) , 
           na.action = "na.exclude", data = black1Mb_stat, 
           correlation = corAR1(form =~ obs))
summary(Z2A)
Z3A <- gls(scaled_pi_Z ~ scale(GC_density), na.action = "na.exclude", data = black1Mb_stat, 
           correlation = corAR1(form =~ obs))
summary(Z3A)

AIC(Z1, Z1A, Z2A, Z3A)
#***********************************************************************************************
# PAR
#***********************************************************************************************
black1000Kb_PAR$obs = as.numeric(rownames(black1000Kb_PAR))
scaled_pi = scale(black1000Kb_PAR$Resampling_mean)
PAR1 <- gls(scaled_pi_PAR ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) + 
            scale(CDS_density), na.action = "na.exclude", data = black1000Kb_PAR, 
            correlation = corAR1(form =~ obs))
summary(PAR1)
PAR1A <- gls(scaled_pi_PAR ~ scale(CM_from_PAR_boundary) + scale(CM_per_bp) + 
               scale(CDS_density), na.action = "na.exclude", data = black1000Kb_PAR, 
             correlation = corAR1(form =~ obs))
summary(PAR1A)
PAR2A <- gls(scaled_pi_PAR ~ scale(CM_from_PAR_boundary) + 
               scale(CDS_density), na.action = "na.exclude", data = black1000Kb_PAR, 
             correlation = corAR1(form =~ obs))
summary(PAR2A)
PAR3A <- gls(scaled_pi_PAR ~ scale(CDS_density), na.action = "na.exclude", data = black1000Kb_PAR, 
             correlation = corAR1(form =~ obs))
summary(PAR3A)

AIC(PAR1, PAR1A, PAR2A, PAR3A)
#***********************************************************************************************
# nonPAR
#***********************************************************************************************
black1000Kb_nonPAR$obs = as.numeric(rownames(black1000Kb_nonPAR))
nonPAR1 <- gls(scaled_pi_nonPAR ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) + 
                scale(CDS_density), na.action = "na.exclude", data = black1000Kb_nonPAR, 
               correlation = corAR1(form =~ obs))
summary(nonPAR1)

nonPAR1A <- gls(scaled_pi_nonPAR ~ scale(CM_from_PAR_boundary) + scale(CM_per_bp) + 
                scale(GC_density), na.action = "na.exclude", data = black1000Kb_nonPAR, 
                correlation = corAR1(form =~ obs))
summary(nonPAR1A)

nonPAR2A <- gls(scaled_pi_nonPAR ~ scale(CM_from_PAR_boundary) +
                scale(GC_density), na.action = "na.exclude", data = black1000Kb_nonPAR, 
                correlation = corAR1(form =~ obs))
summary(nonPAR2A)

nonPAR3A <- gls(scaled_pi_nonPAR ~ scale(GC_density), na.action = "na.exclude", 
                data = black1000Kb_nonPAR, correlation = corAR1(form =~ obs))
summary(nonPAR3A)

AIC(nonPAR1, nonPAR1A, nonPAR2A, nonPAR3A)
#***********************************************************************************************
