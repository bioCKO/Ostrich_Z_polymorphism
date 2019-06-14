#*************************************************************************************
# R script for project ostrich Z polymorphism
# Written by Homa Papoli Yazdi
# January 2019
#*************************************************************************************
#*************************************************************************************
#*************************************************************************************
#*************************************************************************************
# Figure 1: Genetic diversity across Z chromosome
#*************************************************************************************
# Plot black, blue and red intergenic diversity with 
# the confidence interval of the resampling mean and
# the loess regression line and the confidence interval
# of the predicted value.
# About loess:
# 1. Use a type of sliding window to divide the data into smaller blobs.
# 2. At each data point, use a type of least squares to fit a line.
#*************************************************************************************
# Clean the work environment
rm(list=ls())
# Set working environment
#setwd("/proj/uppstore2017180/private/homap/ostrich_Z_diversity/scripts")
#*************************************************************************************
# Read the data
#*************************************************************************************
# Male data
#*************************************************************************************
# Black subspecies
black_intergenic <- read.table("../data/black.500Kb.INTERGENIC.mean.CI.ChrZ.txt")
colnames(black_intergenic) = c("Chr", "Scaffold", "S_Start", "S_End", "Start", "End", 
                               "Resampling_mean", "Low_CI", "Up_CI", "SD")
head(black_intergenic)
# Blue subspecies
blue_intergenic <- read.table("../data/blue.500Kb.INTERGENIC.mean.CI.ChrZ.txt")
colnames(blue_intergenic) = c("Chr", "Scaffold", "S_Start", "S_End", "Start", "End", 
                              "Resampling_mean", "Low_CI", "Up_CI", "SD")
head(blue_intergenic)
# Red subspecies
red_intergenic <- read.table("../data/red.500Kb.INTERGENIC.mean.CI.ChrZ.txt")
colnames(red_intergenic) = c("Chr", "Scaffold", "S_Start", "S_End", "Start", "End", 
                             "Resampling_mean", "Low_CI", "Up_CI", "SD")
head(red_intergenic)

#*************************************************************************************
# Female data
#*************************************************************************************
black_female <- read.table("../data/female_data/black.500Kb.female.INTERGENIC.mean.CI.ChrZ.txt")
head(black_female)
blackPAR1Mb <- black_female[black_female$V5 < 49724557,]
blacknonPAR1Mb <- black_female[black_female$V5 > 53065700,]
mean(blackPAR1Mb$V7)

#*************************************************************************************
# Set graphic parameters
#*************************************************************************************
#pdf("../result/graphs/pi_Z_3_subspecies.pdf", width = 7, height = 8.7)
par(mar=c(6,6,4,4))
layout(matrix(1:3, ncol = 1), widths = 1, heights = c(2.3,2,2.3), respect = FALSE)
par(bty = 'n')
# Plot x axis in Mb
Mb = 1000000
# Find the transparent color for the loess CI
#https://www.rapidtables.com/web/color/RGB_Color.html
Tblue <- rgb(102, 178, 255, max = 255, alpha = 125)
Tred <- rgb(255, 153, 153, max = 255, alpha = 125)
#*************************************************************************************
# Plot each subspecies, CI, loess line and its CI
#*************************************************************************************
# Black subspecies
par(mar = c(0, 4.1, 4.1, 2.1))
plot(black_intergenic$Start/Mb, black_intergenic$Resampling_mean, pch = 20, ylab = "", 
     xlab = "", tck = 0, axes = FALSE)
box(which = "plot", bty = "l")
axis(2, c(0.0004, 0.0008, 0.0012), col = NA, col.ticks = 1)
arrows(black_intergenic$Start/Mb, black_intergenic$Low_CI, black_intergenic$Start/Mb, 
       black_intergenic$Up_CI, length=0.02, angle=90, code=3)
loesspi <- predict(loess(black_intergenic$Resampling_mean~black_intergenic$Start, 
                         span=0.5), se=T)
lines(black_intergenic$Start/Mb,loesspi$fit)
#lines(black_intergenic$V2/Mb, loesspi$fit - qt(0.975,loesspi$df)*loesspi$se, lty=2)
#lines(black_intergenic$V2/Mb, loesspi$fit + qt(0.975,loesspi$df)*loesspi$se, lty=2)
polygon(c(black_intergenic$Start/Mb, rev(black_intergenic$Start/Mb)), 
        c(loesspi$fit + qt(0.975,loesspi$df)*loesspi$se, 
        rev(loesspi$fit - qt(0.975,loesspi$df)*loesspi$se)), col="#99999977", 
        border="#99999977")
# Add the PAR boundary
abline(v = 53.0657, col="dimgrey", lty = 2)

# Blue subspecies
par(mar = c(0, 4.1, 0, 2.1))
plot(blue_intergenic$Start/Mb, blue_intergenic$Resampling_mean, pch = 20, col = "blue", 
     ylim = c(0.0003, 0.0014), ylab = "Diversity (pi)", xlab = "", tck = 0, axes = FALSE)
box(which = "plot", bty = "l")
axis(2, c(0.0004, 0.0008, 0.0012), col = NA, col.ticks = 1)
arrows(blue_intergenic$Start/Mb, blue_intergenic$Low_CI, blue_intergenic$Start/Mb, 
       blue_intergenic$Up_CI, length=0.02, angle=90, code=3, col = "blue")
loesspi <- predict(loess(blue_intergenic$Resampling_mean~blue_intergenic$Start, 
                  span=0.5), se=T)
lines(blue_intergenic$Start/Mb,loesspi$fit, col = "blue")
#lines(black_intergenic$V2/Mb, loesspi$fit - qt(0.975,loesspi$df)*loesspi$se, lty=2)
#lines(black_intergenic$V2/Mb, loesspi$fit + qt(0.975,loesspi$df)*loesspi$se, lty=2)
polygon(c(blue_intergenic$Start/Mb, rev(blue_intergenic$Start/Mb)), 
        c(loesspi$fit + qt(0.975,loesspi$df)*loesspi$se, 
        rev(loesspi$fit - qt(0.975,loesspi$df)*loesspi$se)), col="#66B2FF7D", 
        border="#66B2FF7D")
# Add the PAR boundary
abline(v = 53.0657, col="dimgrey", lty = 2)

par(mar = c(4.1, 4.1, 0, 2.1))
plot(red_intergenic$Start/Mb, red_intergenic$Resampling_mean, pch = 20, col = "red", 
     ylab = "", xlab = "Position (Mb)", tck = 0, axes = FALSE)
box(which = "plot", bty = "l")
axis(1, at = seq(0, 80, 10), col = NA, col.ticks = 1)
axis(2, at = c(0.0002, 0.0006, 0.001), labels = formatC(c(0.0002, 0.0006, 0.001)), 
     col = NA, col.ticks = 1)
arrows(red_intergenic$Start/Mb, red_intergenic$Low_CI, red_intergenic$Start/Mb, 
       red_intergenic$Up_CI, length=0.02, angle=90, code=3, col = "red")
loesspi <- predict(loess(red_intergenic$Resampling_mean~red_intergenic$Start, 
                         span=0.5), se=T)
lines(red_intergenic$Start/Mb,loesspi$fit, col = "red")
#lines(black_intergenic$V2/Mb, loesspi$fit - qt(0.975,loesspi$df)*loesspi$se, lty=2)
#lines(black_intergenic$V2/Mb, loesspi$fit + qt(0.975,loesspi$df)*loesspi$se, lty=2)
polygon(c(red_intergenic$Start/Mb, rev(red_intergenic$Start/Mb)), 
        c(loesspi$fit + qt(0.975,loesspi$df)*loesspi$se, 
        rev(loesspi$fit - qt(0.975,loesspi$df)*loesspi$se)), col="#FF99997D", 
        border = "#FF99997D")
# Add the PAR boundary
abline(v = 53.0657, col="dimgrey", lty = 2)
#dev.off()
#*************************************************************************************
#*************************************************************************************
#*************************************************************************************
# Figure 2: Genetic diversity at the PAR boundary using black subspecies
#*************************************************************************************
#*************************************************************************************
# Read data into R
#*************************************************************************************
# Black subspecies
black200Kb <- read.table("../data/black.200Kb.INTERGENIC.mean.CI.ChrZ.txt")
colnames(black200Kb) = c("Chr", "Scaffold", "S_Start", "S_End", "Start", "End", 
                         "Resampling_mean", "Low_CI", "Up_CI", "SD")
head(black200Kb)
#*************************************************************************************
# Set graphic parameters
#*************************************************************************************
#pdf("../result/graphs/pi_black_PAR_boundary.pdf", width = 6.51, height = 6.04)
#*************************************************************************************
# Black subspecies
PAR_b_black = black200Kb[black200Kb$Scaffold=="superscaffold36",]
plot(PAR_b_black$Start/Mb, PAR_b_black$Resampling_mean, pch = 20, xlab = "Position (Mb)", 
     ylab = "Diversity (Pi)", tck = 0, axes = FALSE)
box(which = "plot", bty = "l")
axis(1, at = seq(44, 56, 2), col = NA, col.ticks = 1, cex.axis = 0.9)
axis(2, at = c(0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.0012), cex.axis = 0.9, 
     labels = formatC(c(0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.0012)), col = NA, col.ticks = 1)
arrows(PAR_b_black$Start/Mb, PAR_b_black$Low_CI, PAR_b_black$Start/Mb, 
       PAR_b_black$Up_CI, length=0.02, angle=90, code=3)
pos_Mb = PAR_b_black$Start/Mb
lm_fit <- predict(lm(PAR_b_black$Resampling_mean~pos_Mb), se = T)
lines(pos_Mb, lm_fit$fit)
lines(pos_Mb, lm_fit$fit - qt(0.975,lm_fit$df)*lm_fit$se, lty=2)
lines(pos_Mb, lm_fit$fit + qt(0.975,lm_fit$df)*lm_fit$se, lty=2)
abline(v = 53.0657, col="dimgrey", lty = 2)

polygon(c(PAR_b_black$Start/Mb, rev(PAR_b_black$Start/Mb)), 
        c(lm_fit$fit + qt(0.975,lm_fit$df)*lm_fit$se, 
       rev(lm_fit$fit - qt(0.975,lm_fit$df)*lm_fit$se)), col="#99999977", 
       border="#99999977")
#dev.off()
#**********************************************************************************************#
#**********************************************************************************************#
#*********************** Statistical analysis using the black subspecies***********************#
#**********************************************************************************************#
#**********************************************************************************************#
rm(list=ls())
library(psych)
#***********************************************************************************************
# Read data into R
#***********************************************************************************************
black1Mb_stat <- read.table("../data/black.1000Kb.INTERGENIC.rec.resampling.Z.forR.txt", header = T)
head(black1Mb_stat)
blue1Mb_stat <- read.table("../data/blue.1000Kb.INTERGENIC.rec.resampling.Z.forR.txt", header = T)
red1Mb_stat <- read.table("../data/red.1000Kb.INTERGENIC.rec.resampling.Z.forR.txt", header = T)
#***********************************************************************************************
# Separating dataset to PAR and non-PAR
blackPAR1Mb <- black1Mb_stat[black1Mb_stat$ChrZ_start < 53065700,]
blacknonPAR1Mb <- black1Mb_stat[black1Mb_stat$ChrZ_start > 53065700,]

bluePAR1Mb <- blue1Mb_stat[blue1Mb_stat$ChrZ_start < 53065700,]
bluenonPAR1Mb <- blue1Mb_stat[blue1Mb_stat$ChrZ_start > 53065700,]

redPAR1Mb <- red1Mb_stat[red1Mb_stat$ChrZ_start < 53065700,]
rednonPAR1Mb <- red1Mb_stat[red1Mb_stat$ChrZ_start > 53065700,]
#***********************************************************************************************
# Exploratory data analysis
#***********************************************************************************************
# Genetic diversity of PAR and non-PAR in males
#***********************************************************************************************
# Black
# Z
mean(black1Mb_stat$pi_per_window)
sd(black1Mb_stat$pi_per_window)
# PAR
mean(blackPAR1Mb$pi_per_window)
sd(blackPAR1Mb$pi_per_window)
# non-PAR
mean(blacknonPAR1Mb$pi_per_window)
sd(blacknonPAR1Mb$pi_per_window)
# PAR, non-PAR t-test
t.test(blackPAR1Mb$pi_per_window, blacknonPAR1Mb$pi_per_window)
# PAR 1-24 Mb
# PAR
mean(blackPAR1Mb$pi_per_window[1:2])


sd(blackPAR1Mb$pi_per_window)


# Blue
# Z
mean(blue1Mb_stat$pi_per_window)
sd(blue1Mb_stat$pi_per_window)
# PAR
mean(bluePAR1Mb$pi_per_window)
sd(bluePAR1Mb$pi_per_window)
# non-PAR
mean(bluenonPAR1Mb$pi_per_window)
sd(bluenonPAR1Mb$pi_per_window)
# PAR, non-PAR t-test
t.test(bluePAR1Mb$pi_per_window, bluenonPAR1Mb$pi_per_window)

# Red
# Z
mean(red1Mb_stat$pi_per_window)
sd(red1Mb_stat$pi_per_window)
# PAR
mean(redPAR1Mb$pi_per_window)
sd(redPAR1Mb$pi_per_window)
# non-PAR
mean(rednonPAR1Mb$pi_per_window)
sd(rednonPAR1Mb$pi_per_window)
# PAR, non-PAR t-test
t.test(redPAR1Mb$pi_per_window, rednonPAR1Mb$pi_per_window)

# Female cM/Mb
mean(black1Mb_stat$female_CM_per_bp*10^6, na.rm = T)
sd(black1Mb_stat$female_CM_per_bp*10^6, na.rm = T)
# PAR
mean(blackPAR1Mb$female_CM_per_bp*10^6, na.rm = T)
sd(blackPAR1Mb$female_CM_per_bp*10^6, na.rm = T)
# non-PAR
mean(blacknonPAR1Mb$female_CM_per_bp*10^6, na.rm = T)
sd(blacknonPAR1Mb$female_CM_per_bp*10^6, na.rm = T)
# 10 Mb before the PAR boundary
mean(blackPAR1Mb$female_CM_per_bp[46:56]*10^6)
mean(blackPAR1Mb$male_CM_per_bp[46:56]*10^6)
t.test(blackPAR1Mb$female_CM_per_bp[46:56]*10^6, blackPAR1Mb$male_CM_per_bp[46:56]*10^6)


# Male cM/Mb
mean(black1Mb_stat$male_CM_per_bp*10^6, na.rm = T)
sd(black1Mb_stat$male_CM_per_bp*10^6, na.rm = T)
# PAR
mean(blackPAR1Mb$male_CM_per_bp*10^6, na.rm = T)
sd(blackPAR1Mb$male_CM_per_bp*10^6, na.rm = T)
# non-PAR
mean(blacknonPAR1Mb$male_CM_per_bp*10^6, na.rm = T)
sd(blacknonPAR1Mb$male_CM_per_bp*10^6, na.rm = T)

# Sex averaged cM/Mb
mean(black1Mb_stat$CM_per_bp*10^6, na.rm = T)
sd(black1Mb_stat$CM_per_bp*10^6, na.rm = T)
# PAR
mean(blackPAR1Mb$CM_per_bp*10^6, na.rm = T)
sd(blackPAR1Mb$CM_per_bp*10^6, na.rm = T)
# non-PAR
mean(blacknonPAR1Mb$CM_per_bp*10^6, na.rm = T)
sd(blacknonPAR1Mb$CM_per_bp*10^6, na.rm = T)

# GC density
mean(black1Mb_stat$GC_density*100, na.rm = T)
sd(black1Mb_stat$GC_density*100, na.rm = T)
# PAR
mean(blackPAR1Mb$GC_density*100, na.rm = T)
sd(blackPAR1Mb$GC_density*100, na.rm = T)
# non-PAR
mean(blacknonPAR1Mb$GC_density*100, na.rm = T)
sd(blacknonPAR1Mb$GC_density*100, na.rm = T)
t.test(blackPAR1Mb$GC_density*10, blacknonPAR1Mb$GC_density*10)

# Gene density
mean(black1Mb_stat$CDS_density*100, na.rm = T)
sd(black1Mb_stat$CDS_density*100, na.rm = T)
# PAR
mean(blackPAR1Mb$CDS_density*100, na.rm = T)
sd(blackPAR1Mb$CDS_density*100, na.rm = T)
# non-PAR
mean(blacknonPAR1Mb$CDS_density*100, na.rm = T)
sd(blacknonPAR1Mb$CDS_density*100, na.rm = T)
t.test(blackPAR1Mb$CDS_density*100, na.rm = T, blacknonPAR1Mb$CDS_density*100, na.rm = T)

# Repeat density
mean(black1Mb_stat$repeat_density*100, na.rm = T)
sd(black1Mb_stat$repeat_density*100, na.rm = T)
# PAR
mean(blackPAR1Mb$repeat_density*100, na.rm = T)
sd(blackPAR1Mb$repeat_density*100, na.rm = T)
# non-PAR
mean(blacknonPAR1Mb$repeat_density*100, na.rm = T)
sd(blacknonPAR1Mb$repeat_density*100, na.rm = T)
t.test(blackPAR1Mb$repeat_density*100, blacknonPAR1Mb$repeat_density*100)

# ANOVA between black, blue and red subspecies
het_data <- rbind(as.matrix(black1Mb_stat$pi_per_window), as.matrix(blue1Mb_stat$pi_per_window), 
      as.matrix(red1Mb_stat$pi_per_window))
group_data <- rbind(as.matrix(rep("black", 85)), as.matrix(rep("blue", 85)), as.matrix(rep("red", 85)))
anova_df <- as.data.frame(cbind(het_data, group_data))
colnames(anova_df) <- c("pi", "group")
res.aov <- aov(as.numeric(pi) ~ group, data = anova_df)
summary(res.aov)
TukeyHSD(res.aov)
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

# Pairwise correlation significance test
# Z-GC
cor.test(black1Mb_stat$pi_per_window, black1Mb_stat$GC_density)
# PAR-GC
cor.test(blackPAR1Mb$pi_per_window, blackPAR1Mb$GC_density)
# nonPAR-GC
cor.test(blacknonPAR1Mb$pi_per_window, blacknonPAR1Mb$GC_density)

# Z-cM/Mb
cor.test(black1Mb_stat$pi_per_window, log2(black1Mb_stat$CM_per_bp+1))
# PAR-cM/Mb
cor.test(blackPAR1Mb$pi_per_window, log2(blackPAR1Mb$CM_per_bp+1))
# nonPAR-cM/Mb
cor.test(blacknonPAR1Mb$pi_per_window, log2(blacknonPAR1Mb$CM_per_bp+1))

# Z-genetic distance
cor.test(black1Mb_stat$pi_per_window, black1Mb_stat$CM_from_PAR_boundary)
# PAR-genetic distance
cor.test(blackPAR1Mb$pi_per_window, blackPAR1Mb$CM_from_PAR_boundary)
# nonPAR-genetic distance
cor.test(blacknonPAR1Mb$pi_per_window, blacknonPAR1Mb$CM_from_PAR_boundary)

# Gene density and recombination rate
cor.test(black1Mb_stat$CDS_density, black1Mb_stat$CM_per_bp)
cor.test(black1Mb_stat$CDS_density, black1Mb_stat$GC_density)

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
# Correlation among significantly correlated explanatory variables with diversity
#***********************************************************************************************
#***********************************************************************************************
# 1. CM_from_PAR_boundary

# CM_from_PAR_boundary is correlated with GC_density. Biological explanation for this correlation
# is difficult unless areas furthest away from the PAR boundary also experience highest rate of
# recombiantion. This is plausible given that areas furthest away from the PAR boundary are closer
# to the chromosomes' end which experience higher rates of recombination. We can check this using
# coplot of GC_density as a function of CM_from_PAR_boundary in different levels of recombination
# rate (CM_per_bp).
#***********************************************************************************************
# For PAR
#***********************************************************************************************
given.CM_per_bp <- co.intervals(black1000Kb_PAR$CM_per_bp, number = 3)
coplot(GC_density ~ CM_from_PAR_boundary | CM_per_bp, data = black1000Kb_PAR, 
       given.v = given.CM_per_bp, rows = 1)
# For PAR it is not quite clear that higher GC values further away from the PAR boundary coincide
# with higher recombination rate.
#***********************************************************************************************
# For nonPAR
#***********************************************************************************************
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
# GC content in avian coding sequence as shown in various neognathae species is higher than in
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
# (non-PAR Z spends 1/3 of the time in females and 2/3 of the time in male).

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
dwtest(full_Z)
library(DataCombine)
econ_data <- data.frame(black1Mb_stat, resid_mod1=resid(full_Z))
econ_data_1 <- slide(econ_data, Var="resid_mod1", NewVar = "lag1", slideBy = -1)
econ_data_2 <- na.omit(econ_data_1)
full_Z <- lm(scaled_pi_Z ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) + 
            scale(CDS_density) + lag1, na.action = na.omit, data = econ_data_1)
E <- residuals(full_Z)
acf(E)
dwtest(full_Z)
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
black1Mb_stat$obs = as.numeric(rownames(black1Mb_stat))

Z0 <- gls(scaled_pi_Z ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) +
          scale(CDS_density), na.action = "na.exclude", method = "ML", data = black1Mb_stat)

Z1 <- gls(scaled_pi_Z ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) +
            scale(CDS_density), na.action = na.omit, method = "ML", data = black1Mb_stat, 
            correlation = corAR1(form =~ obs))
Z2 <- gls(scaled_pi_Z ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) +
          scale(CDS_density), na.action = na.omit, method = "ML", data = black1Mb_stat, 
          correlation = corSpher(form =~ obs))
Z3 <- gls(scaled_pi_Z ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) +
          scale(CDS_density), na.action = na.omit, method = "ML", data = black1Mb_stat, 
          correlation = corLin(form =~ obs))
Z4 <- gls(scaled_pi_Z ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) +
          scale(CDS_density), na.action = na.omit, method = "ML", data = black1Mb_stat, 
          correlation = corRatio(form =~ obs))
Z5 <- gls(scaled_pi_Z ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) +
          scale(CDS_density), na.action = na.omit, method = "ML", data = black1Mb_stat, 
          correlation = corGaus(form =~ obs))
Z6 <- gls(scaled_pi_Z ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) +
          scale(CDS_density), na.action = na.omit, method = "ML", data = black1Mb_stat, 
          correlation = corExp(form =~ obs))

AIC(Z0, Z1, Z2, Z3, Z4, Z5, Z6)

# All models with spatial structure have a lower AIC, showing that adding a spatial correlation
# structure improves the model. This is particularly the case for models with corAR1, corExp and
# corSpher and corRatio.

# Z chromosome
Z1 <- gls(scaled_pi_Z ~ scale(GC_density) + scale(CM_per_bp) +
            scale(CDS_density), na.action = na.omit, method = "ML", data = black1Mb_stat, 
          correlation = corAR1(form =~ obs))
summary(Z1)

# Full model has the least AIC
#***********************************************************************************************
# PAR
#***********************************************************************************************
black1000Kb_PAR$obs = as.numeric(rownames(black1000Kb_PAR))
scaled_pi_PAR = scale(black1000Kb_PAR$Resampling_mean)
PAR1 <- gls(scaled_pi_PAR ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(CM_per_bp) + 
            scale(CDS_density), na.action = "na.exclude", method = "ML", data = black1000Kb_PAR, 
            correlation = corAR1(form =~ obs))
summary(PAR1)
#***********************************************************************************************
# nonPAR
#***********************************************************************************************
black1000Kb_nonPAR$obs = as.numeric(rownames(black1000Kb_nonPAR))
scaled_pi_nonPAR = scale(black1000Kb_nonPAR$Resampling_mean)
nonPAR1 <- gls(scaled_pi_nonPAR ~ scale(CM_from_PAR_boundary) + scale(GC_density) + scale(male_CM_per_bp) + 
            scale(CDS_density), na.action = "na.exclude", method = "ML", data = black1000Kb_nonPAR, 
            correlation = corAR1(form =~ obs))
summary(nonPAR1)

nonPAR1A <- gls(scaled_pi_nonPAR ~ scale(CM_from_PAR_boundary) + scale(CM_per_bp) + 
                scale(GC_density), na.action = "na.exclude", method = "ML", data = black1000Kb_nonPAR, 
               correlation = corAR1(form =~ obs))
summary(nonPAR1A)
#***********************************************************************************************
# Plot genetic diversity as a function of genetic distance to the pseudo-autosomal boundary
#***********************************************************************************************
#pdf("../result/graphs/cM_from_PAR_boundary.pdf", width = 5.78, height = 7.63)
par(mar=c(6,6,4,4))
layout(matrix(1:3, ncol = 1), widths = 1, heights = c(2.3,2,2.3), respect = FALSE)
par(bty = 'n')

par(mar = c(0, 4.1, 4.1, 2.1))
plot(x=black1000Kb_PAR$CM_from_PAR_boundary, y=black1000Kb_PAR$Resampling_mean, ylab = "", 
     xlab = "", tck = 0, axes = FALSE, pch = 20, col = "black")
loesspi <- predict(loess(black1000Kb_PAR$Resampling_mean~black1000Kb_PAR$CM_from_PAR_boundary, 
                   span=0.5), se=T)
lines(black1000Kb_PAR$CM_from_PAR_boundary,loesspi$fit)
box(which = "plot", bty = "l")
axis(2, c(0.0004, 0.0008, 0.0012), col = NA, col.ticks = 1)

par(mar = c(0, 4.1, 0, 2.1))
plot(x=bluePAR1Mb$CM_from_PAR_boundary, y=bluePAR1Mb$Resampling_mean, ylab = "Diversity (pi)", 
     xlab = "", ylim = c(0.00038, 0.0014), tck = 0, axes = FALSE, pch = 20, col = "blue")
loesspi <- predict(loess(bluePAR1Mb$Resampling_mean~bluePAR1Mb$CM_from_PAR_boundary, span=0.5), se=T)
lines(bluePAR1Mb$CM_from_PAR_boundary,loesspi$fit)
box(which = "plot", bty = "l")
axis(2, c(0.0004, 0.0008, 0.0012), col = NA, col.ticks = 1)

par(mar = c(4.1, 4.1, 0, 2.1))
plot(x=redPAR1Mb$CM_from_PAR_boundary, y=redPAR1Mb$Resampling_mean, 
               , ylab = "", xlab = "cM from PAR boundary", pch = 20, col = "red", tck = 0, axes = FALSE)
loesspi <- predict(loess(redPAR1Mb$Resampling_mean~redPAR1Mb$CM_from_PAR_boundary, span=0.5), se=T)
lines(redPAR1Mb$CM_from_PAR_boundary,loesspi$fit)
box(which = "plot", bty = "l")
axis(1, at = seq(0, 150, 50), col = NA, col.ticks = 1)
axis(2, at = c(0.0002, 0.0006, 0.001), labels = formatC(c(0.0002, 0.0006, 0.001)), col = NA, col.ticks = 1)
#dev.off()

#***********************************************************************************************
# FST between males and females in the pseudo-autosomal boundary
#***********************************************************************************************
superscaffold26 =  25310599
superscaffold54 = superscaffold26 + 5000 + 16379243
superscaffold35 = superscaffold54 + 5000 + 4625539
superscaffold36 = superscaffold35 + 5000 + 9394175
superscaffold62 = superscaffold36 + 5000 + 2917291
superscaffold67 = superscaffold62 + 5000 + 5300260
superscaffold69_1 = superscaffold67 + 5000 + 5978518
superscaffold93 = superscaffold69_1 + 5000 + 4983591
superscaffold63 = superscaffold93 + 5000 + 1692925
superscaffold88 = superscaffold63 + 5000 + 624114
superscaffold83 = superscaffold88 + 5000 + 782506
superscaffold92 = superscaffold83 + 5000 + 2882843

# V8 is the weigthed Fst
black_100 <- read.table("../data/black_male_female_100Kb.windowed.weir.Z.fst")
black_100_0 <- black_100[which(black_100$V8<0),]
black_100_0$V8 = 0
black_100_no0 <- black_100[black_100$V8>0,]
black_100_modified <- rbind(black_100_0, black_100_no0)

black_200 <- read.table("../data/black_male_female_200Kb.windowed.weir.Z.fst")
black_200_0 <- black_200[which(black_200$V8<0),]
black_200_0$V8 = 0
black_200_no0 <- black_200[black_200$V8>0,]
black_200_modified <- rbind(black_200_0, black_200_no0)

# Black, male-female
black_200_modified_PAR <- black_200_modified[black_200_modified$V5 < 53065700,]
black_200_modified_nonPAR <- black_200_modified[black_200_modified$V5 > 53065700,]

#*************************************************************************************
# Z
mean(black_200_modified$V8, na.rm = T)
sd(black_200_modified$V8, na.rm = T)
# PAR
mean(black_200_modified_PAR$V8, na.rm = T)
sd(black_200_modified_PAR$V8, na.rm = T)
# non-PAR
mean(black_200_modified_nonPAR$V8, na.rm = T)
sd(black_200_modified_nonPAR$V8, na.rm = T)
#*************************************************************************************

Mb = 1000000
par(bty = 'n')
plot(black_100_modified$V5/1000000, black_100_modified$V8, pch = 20, 
     ylab = "Fst-100Kb", xlab = "Position (Mb)")
abline(v = 53.0657, col="red", lty = 2)
abline(v = superscaffold26/Mb, col="blue", lty = 3)
abline(v = superscaffold54/Mb, col="blue", lty = 3)
abline(v = superscaffold35/Mb, col="blue", lty = 3)
abline(v = superscaffold36/Mb, col="blue", lty = 3)
abline(v = superscaffold62/Mb, col="blue", lty = 3)
abline(v = superscaffold67/Mb, col="blue", lty = 3)
abline(v = superscaffold69_1/Mb, col="blue", lty = 3)
abline(v = superscaffold93/Mb, col="blue", lty = 3)
abline(v = superscaffold63/Mb, col="blue", lty = 3)
abline(v = superscaffold88/Mb, col="blue", lty = 3)
abline(v = superscaffold83/Mb, col="blue", lty = 3)
abline(v = superscaffold92/Mb, col="blue", lty = 3)

fst_resample <- rep(0, 1000)
for(i in seq(1, 1000)){
  fst_resample[i] <-mean(sample(black_100_modified$V8, length(black_100_modified$V8), replace = TRUE))
}

hist(fst_resample)

par(mfrow=c(1,1))
par(bty = 'n')
plot(black_200_modified$V5/1000000, black_200_modified$V8, pch = 20, 
     ylab = "Male_female_Fst-200Kb", xlab = "Position (Mb)")
axis(1, at = seq(0, 80, 10), col = NA, col.ticks = 1)
axis(2, seq(0, 1, 0.2), col = NA, col.ticks = 1)
box(which = "plot", bty = "l")


abline(v = 53.0657, col="red", lty = 2)
abline(v = superscaffold26/Mb, col="blue", lty = 3)
abline(v = superscaffold54/Mb, col="blue", lty = 3)
abline(v = superscaffold35/Mb, col="blue", lty = 3)
abline(v = superscaffold36/Mb, col="blue", lty = 3)
abline(v = superscaffold62/Mb, col="blue", lty = 3)
abline(v = superscaffold67/Mb, col="blue", lty = 3)
abline(v = superscaffold69_1/Mb, col="blue", lty = 3)
abline(v = superscaffold93/Mb, col="blue", lty = 3)
abline(v = superscaffold63/Mb, col="blue", lty = 3)
abline(v = superscaffold88/Mb, col="blue", lty = 3)
abline(v = superscaffold83/Mb, col="blue", lty = 3)
abline(v = superscaffold92/Mb, col="blue", lty = 3)

blue_100 <- read.table("../data/blue_male_female_100Kb.windowed.weir.Z.fst")
blue_100_0 <- blue_100[which(blue_100$V8<0),]
blue_100_0$V8 = 0
blue_100_no0 <- blue_100[blue_100$V8>0,]
blue_100_modified <- rbind(blue_100_0, blue_100_no0)

blue_200 <- read.table("../data/blue_male_female_200Kb.windowed.weir.Z.fst")
blue_200_0 <- blue_200[which(blue_200$V8<0),]
blue_200_0$V8 = 0
blue_200_no0 <- blue_200[blue_200$V8>0,]
blue_200_modified <- rbind(blue_200_0, blue_200_no0)


# Blue, male-female
blue_200_modified_PAR <- blue_200_modified[blue_200_modified$V5 < 53065700,]
blue_200_modified_nonPAR <- blue_200_modified[blue_200_modified$V5 > 53065700,]
#*************************************************************************************
# Z
mean(blue_200_modified$V8, na.rm = T)
# PAR
mean(blue_200_modified_PAR$V8, na.rm = T)
# non-PAR
mean(blue_200_modified_nonPAR$V8, na.rm = T)
#*************************************************************************************


Mb = 1000000
plot(blue_100_modified$V5/1000000, blue_100_modified$V8, pch = 20, 
     ylab = "Fst-100Kb", xlab = "Position (Mb)")
abline(v = 53.0657, col="red", lty = 2)
abline(v = superscaffold26/Mb, col="blue", lty = 3)
abline(v = superscaffold54/Mb, col="blue", lty = 3)
abline(v = superscaffold35/Mb, col="blue", lty = 3)
abline(v = superscaffold36/Mb, col="blue", lty = 3)
abline(v = superscaffold62/Mb, col="blue", lty = 3)
abline(v = superscaffold67/Mb, col="blue", lty = 3)
abline(v = superscaffold69_1/Mb, col="blue", lty = 3)
abline(v = superscaffold93/Mb, col="blue", lty = 3)
abline(v = superscaffold63/Mb, col="blue", lty = 3)
abline(v = superscaffold88/Mb, col="blue", lty = 3)
abline(v = superscaffold83/Mb, col="blue", lty = 3)
abline(v = superscaffold92/Mb, col="blue", lty = 3)

plot(blue_200_modified$V5/1000000, blue_200_modified$V8, pch = 20, 
     ylab = "Fst-200Kb", xlab = "Position (Mb)")
abline(v = 53.0657, col="red", lty = 2)
abline(v = superscaffold26/Mb, col="blue", lty = 3)
abline(v = superscaffold54/Mb, col="blue", lty = 3)
abline(v = superscaffold35/Mb, col="blue", lty = 3)
abline(v = superscaffold36/Mb, col="blue", lty = 3)
abline(v = superscaffold62/Mb, col="blue", lty = 3)
abline(v = superscaffold67/Mb, col="blue", lty = 3)
abline(v = superscaffold69_1/Mb, col="blue", lty = 3)
abline(v = superscaffold93/Mb, col="blue", lty = 3)
abline(v = superscaffold63/Mb, col="blue", lty = 3)
abline(v = superscaffold88/Mb, col="blue", lty = 3)
abline(v = superscaffold83/Mb, col="blue", lty = 3)
abline(v = superscaffold92/Mb, col="blue", lty = 3)


red_100 <- read.table("../data/red_male_female_100Kb.windowed.weir.Z.fst")
red_100_0 <- red_100[which(red_100$V8<0),]
red_100_0$V8 = 0
red_100_no0 <- red_100[red_100$V8>0,]
red_100_modified <- rbind(red_100_0, red_100_no0)

red_200 <- read.table("../data/red_male_female_200Kb.windowed.weir.Z.fst")
red_200_0 <- red_200[which(red_200$V8<0),]
red_200_0$V8 = 0
red_200_no0 <- red_200[red_200$V8>0,]
red_200_modified <- rbind(red_200_0, red_200_no0)

# Blue, male-female
red_200_modified_PAR <- red_200_modified[red_200_modified$V5 < 53065700,]
red_200_modified_nonPAR <- red_200_modified[blue_200_modified$V5 > 53065700,]
#*************************************************************************************
# Z
mean(red_200_modified$V8, na.rm = T)
# PAR
mean(red_200_modified_PAR$V8, na.rm = T)
# non-PAR
mean(red_200_modified_nonPAR$V8, na.rm = T)
#*************************************************************************************


Mb = 1000000
plot(red_100_modified$V5/1000000, red_100_modified$V8, pch = 20, 
     ylab = "Fst-100Kb", xlab = "Position (Mb)")
abline(v = 53.0657, col="red", lty = 2)
abline(v = superscaffold26/Mb, col="blue", lty = 3)
abline(v = superscaffold54/Mb, col="blue", lty = 3)
abline(v = superscaffold35/Mb, col="blue", lty = 3)
abline(v = superscaffold36/Mb, col="blue", lty = 3)
abline(v = superscaffold62/Mb, col="blue", lty = 3)
abline(v = superscaffold67/Mb, col="blue", lty = 3)
abline(v = superscaffold69_1/Mb, col="blue", lty = 3)
abline(v = superscaffold93/Mb, col="blue", lty = 3)
abline(v = superscaffold63/Mb, col="blue", lty = 3)
abline(v = superscaffold88/Mb, col="blue", lty = 3)
abline(v = superscaffold83/Mb, col="blue", lty = 3)
abline(v = superscaffold92/Mb, col="blue", lty = 3)

plot(red_200_modified$V5/1000000, red_200_modified$V8, pch = 20, 
     ylab = "Fst-200Kb", xlab = "Position (Mb)")
abline(v = 53.0657, col="red", lty = 2)
abline(v = superscaffold26/Mb, col="blue", lty = 3)
abline(v = superscaffold54/Mb, col="blue", lty = 3)
abline(v = superscaffold35/Mb, col="blue", lty = 3)
abline(v = superscaffold36/Mb, col="blue", lty = 3)
abline(v = superscaffold62/Mb, col="blue", lty = 3)
abline(v = superscaffold67/Mb, col="blue", lty = 3)
abline(v = superscaffold69_1/Mb, col="blue", lty = 3)
abline(v = superscaffold93/Mb, col="blue", lty = 3)
abline(v = superscaffold63/Mb, col="blue", lty = 3)
abline(v = superscaffold88/Mb, col="blue", lty = 3)
abline(v = superscaffold83/Mb, col="blue", lty = 3)
abline(v = superscaffold92/Mb, col="blue", lty = 3)

par(mfrow=c(3,1))
black_blue_200 <- read.table("../data/black_blue_male_200Kb.Z")
black_blue_200_0 <- black_blue_200[which(black_blue_200$V8<0),]
black_blue_200_0$V8 = 0
black_blue_200_no0 <- black_blue_200[black_blue_200$V8>0,]
black_blue_200_modified <- rbind(black_blue_200_0, black_blue_200_no0)

black_red_200 <- read.table("../data/black_red_male_200Kb.Z")
black_red_200_0 <- black_red_200[which(black_red_200$V8<0),]
black_red_200_0$V8 = 0
black_red_200_no0 <- black_red_200[black_red_200$V8>0,]
black_red_200_modified <- rbind(black_red_200_0, black_red_200_no0)

blue_red_200 <- read.table("../data/blue_red_male_200Kb.Z")
#blue_red_200_0 <- blue_red_200[which(blue_red_200$V8<0),]
#blue_red_200_0$V8 = 0
#blue_red_200_no0 <- blue_red_200[blue_red_200$V8>0,]
#blue_red_200_modified <- rbind(blue_red_200_0, blue_red_200_no0)

plot(black_blue_200_modified$V5/1000000, black_blue_200_modified$V8, pch = 20, 
     ylab = "Black-Blue-Fst-200Kb", xlab = "Position (Mb)")
mean(black_blue_200$V8)
abline(v = 53.0657, col="red", lty = 2)
abline(v = superscaffold26/Mb, col="blue", lty = 3)
abline(v = superscaffold54/Mb, col="blue", lty = 3)
abline(v = superscaffold35/Mb, col="blue", lty = 3)
abline(v = superscaffold36/Mb, col="blue", lty = 3)
abline(v = superscaffold62/Mb, col="blue", lty = 3)
abline(v = superscaffold67/Mb, col="blue", lty = 3)
abline(v = superscaffold69_1/Mb, col="blue", lty = 3)
abline(v = superscaffold93/Mb, col="blue", lty = 3)
abline(v = superscaffold63/Mb, col="blue", lty = 3)
abline(v = superscaffold88/Mb, col="blue", lty = 3)
abline(v = superscaffold83/Mb, col="blue", lty = 3)
abline(v = superscaffold92/Mb, col="blue", lty = 3)
plot(black_red_200$V5/1000000, black_red_200$V8, pch = 20, 
     ylab = "Black-Red-Fst-200Kb", xlab = "Position (Mb)")
mean(black_red_200$V8)
abline(v = 53.0657, col="red", lty = 2)
abline(v = superscaffold26/Mb, col="blue", lty = 3)
abline(v = superscaffold54/Mb, col="blue", lty = 3)
abline(v = superscaffold35/Mb, col="blue", lty = 3)
abline(v = superscaffold36/Mb, col="blue", lty = 3)
abline(v = superscaffold62/Mb, col="blue", lty = 3)
abline(v = superscaffold67/Mb, col="blue", lty = 3)
abline(v = superscaffold69_1/Mb, col="blue", lty = 3)
abline(v = superscaffold93/Mb, col="blue", lty = 3)
abline(v = superscaffold63/Mb, col="blue", lty = 3)
abline(v = superscaffold88/Mb, col="blue", lty = 3)
abline(v = superscaffold83/Mb, col="blue", lty = 3)
abline(v = superscaffold92/Mb, col="blue", lty = 3)
plot(blue_red_200$V5/1000000, blue_red_200$V8, pch = 20, 
     ylab = "Blue-Red-Fst-200Kb", xlab = "Position (Mb)")
mean(blue_red_200$V8)
abline(v = 53.0657, col="red", lty = 2)
abline(v = superscaffold26/Mb, col="blue", lty = 3)
abline(v = superscaffold54/Mb, col="blue", lty = 3)
abline(v = superscaffold35/Mb, col="blue", lty = 3)
abline(v = superscaffold36/Mb, col="blue", lty = 3)
abline(v = superscaffold62/Mb, col="blue", lty = 3)
abline(v = superscaffold67/Mb, col="blue", lty = 3)
abline(v = superscaffold69_1/Mb, col="blue", lty = 3)
abline(v = superscaffold93/Mb, col="blue", lty = 3)
abline(v = superscaffold63/Mb, col="blue", lty = 3)
abline(v = superscaffold88/Mb, col="blue", lty = 3)
abline(v = superscaffold83/Mb, col="blue", lty = 3)
abline(v = superscaffold92/Mb, col="blue", lty = 3)

par(mfrow=c(1,1))
boxplot(black_blue_200$V9, black_red_200$V9, blue_red_200$V9, 
        names = c("Black-Blue", "Black-Red", "Blue-Red"), ylab = "Fst-200Kb")

# WEIGHTED_FST, MEAN_FST

# Black-Blue Fst
black_blue_200_PAR <- black_blue_200[black_blue_200$V5 < 53065700,]
black_blue_200_nonPAR <- black_blue_200[black_blue_200$V5 > 53065700,]

# Z
mean(black_blue_200$V8, na.rm = T)
# PAR
mean(black_blue_200_PAR$V8, na.rm = T)
# non-PAR 
mean(black_blue_200_nonPAR$V8, na.rm = T)

# Black-Red Fst
black_red_200_PAR <- black_red_200[black_red_200$V5 < 53065700,]
black_red_200_nonPAR <- black_red_200[black_red_200$V5 > 53065700,]

# Z
mean(black_red_200$V8, na.rm = T)
# PAR
mean(black_red_200_PAR$V8, na.rm = T)
sd(black_red_200_PAR$V8, na.rm = T)
# non-PAR 
mean(black_red_200_nonPAR$V8, na.rm = T)
sd(black_red_200_nonPAR$V8, na.rm = T)

t.test(black_red_200_nonPAR$V8, black_red_200_PAR$V8)

# Blue-Red Fst
blue_red_200_PAR <- blue_red_200[blue_red_200$V5 < 53065700,]
blue_red_200_nonPAR <- blue_red_200[blue_red_200$V5 > 53065700,]

# Z
mean(blue_red_200$V8, na.rm = T)
# PAR
mean(blue_red_200_PAR$V8, na.rm = T)
# non-PAR 
mean(blue_red_200_nonPAR$V8, na.rm = T)

t.test(blue_red_200_PAR$V8, blue_red_200_nonPAR$V8)
#***********************************************************************************************
# Tajima's D
#***********************************************************************************************
black_tajimaD <- read.table("../data/black_male_200Kb.allsites.Tajima.D.chrZ.txt")
blue_tajimaD <- read.table("../data/blue_male_200Kb.allsites.Tajima.D.chrZ.txt")
red_tajimaD <- read.table("../data/red_male_200Kb.allsites.Tajima.D.chrZ.txt")

black_tajimaD_PAR <- black_tajimaD[black_tajimaD$V5 < 53065700,]
black_tajimaD_nonPAR <- black_tajimaD[black_tajimaD$V5 > 53065700,]

blue_tajimaD_PAR <- blue_tajimaD[blue_tajimaD$V5 < 53065700,]
blue_tajimaD_nonPAR <- blue_tajimaD[blue_tajimaD$V5 > 53065700,]

red_tajimaD_PAR <- red_tajimaD[red_tajimaD$V5 < 53065700,]
red_tajimaD_nonPAR <- red_tajimaD[red_tajimaD$V5 > 53065700,]

mean(black_tajimaD_PAR$V7, na.rm = T)
mean(black_tajimaD_nonPAR$V7, na.rm = T)

mean(blue_tajimaD_PAR$V7, na.rm = T)
mean(blue_tajimaD_nonPAR$V7, na.rm = T)

mean(red_tajimaD_PAR$V7, na.rm = T)
mean(red_tajimaD_nonPAR$V7, na.rm = T)
#***********************************************************************************************
# Linkage Disequilibrium
#***********************************************************************************************
black_LD <- read.table("../../OS_Z_301018/result/LD_chromosome_plot/black.Z.coordinate.LD.05-50.200kbBin.50.kbStep.Int.v2.out")
blue_LD <- read.table("../../OS_Z_301018/result/LD_chromosome_plot/blue.Z.coordinate.LD.05-50.200kbBin.50.kbStep.Int.v2.out")
red_LD <- read.table("../../OS_Z_301018/result/LD_chromosome_plot/red.Z.coordinate.LD.05-50.200kbBin.50.kbStep.Int.v2.out")


colnames(x = black_LD) = c("Chromosome", "Window_start", "Window_end", "Dprime_LD",
                           "LOD_Dprime", "r_squared", "CIlow_rsq", "CIHi_rsq",
                           "mean_SNP_distance", "NUM_snp_pairs")

colnames(x = blue_LD) = c("Chromosome", "Window_start", "Window_end", "Dprime_LD",
                          "LOD_Dprime", "r_squared", "CIlow_rsq", "CIHi_rsq",
                          "mean_SNP_distance", "NUM_snp_pairs")

colnames(x = red_LD) = c("Chromosome", "Window_start", "Window_end", "Dprime_LD",
                         "LOD_Dprime", "r_squared", "CIlow_rsq", "CIHi_rsq",
                         "mean_SNP_distance", "NUM_snp_pairs")

# LD for PAR and non-PAR
black_LDPAR <- black_LD[black_LD$Window_start < 53065700,]
black_LDnonPAR <- black_LD[black_LD$Window_start > 53065700,]

blue_LDPAR <- blue_LD[blue_LD$Window_start < 53065700,]
blue_LDnonPAR <- blue_LD[blue_LD$Window_start > 53065700,]

red_LDPAR <- red_LD[red_LD$Window_start < 53065700,]
red_LDnonPAR <- red_LD[red_LD$Window_start > 53065700,]

# Black mean and sd
mean(black_LD$r_squared, na.rm = T)
sd(black_LD$r_squared, na.rm = T)

mean(black_LDPAR$r_squared)
sd(black_LDPAR$r_squared)

mean(black_LDnonPAR$r_squared, na.rm = T)
sd(black_LDnonPAR$r_squared, na.rm = T)

t.test(black_LDPAR$r_squared, black_LDnonPAR$r_squared)

# Blue mean and sd
mean(blue_LD$r_squared, na.rm = T)
sd(blue_LD$r_squared, na.rm = T)

mean(blue_LDPAR$r_squared)
sd(blue_LDPAR$r_squared)

mean(blue_LDnonPAR$r_squared, na.rm = T)
sd(blue_LDnonPAR$r_squared, na.rm = T)

t.test(blue_LDPAR$r_squared, blue_LDnonPAR$r_squared)

# Red mean and sd
mean(red_LD$r_squared, na.rm = T)
sd(red_LD$r_squared, na.rm = T)

mean(red_LDPAR$r_squared)
sd(red_LDPAR$r_squared)

mean(red_LDnonPAR$r_squared, na.rm = T)
sd(red_LDnonPAR$r_squared, na.rm = T)

t.test(red_LDPAR$r_squared, red_LDnonPAR$r_squared)


Mb = 1000000
par(mfrow=c(3,1))
superscaffold26 =  25310599
superscaffold54 = superscaffold26 + 5000 + 16379243
superscaffold35 = superscaffold54 + 5000 + 4625539
superscaffold36 = superscaffold35 + 5000 + 9394175
superscaffold62 = superscaffold36 + 5000 + 2917291
superscaffold67 = superscaffold62 + 5000 + 5300260
superscaffold69_1 = superscaffold67 + 5000 + 5978518
superscaffold93 = superscaffold69_1 + 5000 + 4983591
superscaffold63 = superscaffold93 + 5000 + 1692925
superscaffold88 = superscaffold63 + 5000 + 624114
superscaffold83 = superscaffold88 + 5000 + 782506
superscaffold92 = superscaffold83 + 5000 + 2882843
par(bty = 'n')
plot((black_LD$Window_start+black_LD$Window_end)/2000000, black_LD$r_squared, 
     xlab="Position (Mb)", ylab="LD (r^2)", pch = 20, tck = 0, axes = FALSE)
axis(1, at = seq(0, 80, 10), col = NA, col.ticks = 1)
axis(2, seq(0.1, 0.85, 0.1), col = NA, col.ticks = 1)
box(which = "plot", bty = "l")

abline(v = 53.0657, col="red", lty = 2)
abline(v = superscaffold26/Mb, col="blue", lty = 3)
abline(v = superscaffold54/Mb, col="blue", lty = 3)
abline(v = superscaffold35/Mb, col="blue", lty = 3)
abline(v = superscaffold36/Mb, col="blue", lty = 3)
abline(v = superscaffold62/Mb, col="blue", lty = 3)
abline(v = superscaffold67/Mb, col="blue", lty = 3)
abline(v = superscaffold69_1/Mb, col="blue", lty = 3)
abline(v = superscaffold93/Mb, col="blue", lty = 3)
abline(v = superscaffold63/Mb, col="blue", lty = 3)
abline(v = superscaffold88/Mb, col="blue", lty = 3)
abline(v = superscaffold83/Mb, col="blue", lty = 3)
abline(v = superscaffold92/Mb, col="blue", lty = 3)
plot((blue_LD$Window_start+blue_LD$Window_end)/2000000, blue_LD$r_squared, 
     xlab="Position (Mb)", ylab="LD (r^2)", pch = 20, tck = 0, axes = FALSE, col = "blue")
axis(1, at = seq(0, 80, 10), col = NA, col.ticks = 1)
axis(2, seq(0.1, 0.85, 0.1), col = NA, col.ticks = 1)
box(which = "plot", bty = "l")
abline(v = 53.0657, col="red", lty = 2)
abline(v = superscaffold26/Mb, col="blue", lty = 3)
abline(v = superscaffold54/Mb, col="blue", lty = 3)
abline(v = superscaffold35/Mb, col="blue", lty = 3)
abline(v = superscaffold36/Mb, col="blue", lty = 3)
abline(v = superscaffold62/Mb, col="blue", lty = 3)
abline(v = superscaffold67/Mb, col="blue", lty = 3)
abline(v = superscaffold69_1/Mb, col="blue", lty = 3)
abline(v = superscaffold93/Mb, col="blue", lty = 3)
abline(v = superscaffold63/Mb, col="blue", lty = 3)
abline(v = superscaffold88/Mb, col="blue", lty = 3)
abline(v = superscaffold83/Mb, col="blue", lty = 3)
abline(v = superscaffold92/Mb, col="blue", lty = 3)
plot((red_LD$Window_start+red_LD$Window_end)/2000000, red_LD$r_squared, 
     xlab="Position (Mb)", ylab="LD (r^2)", pch = 20, tck = 0, axes = FALSE, col = "red")
axis(1, at = seq(0, 80, 10), col = NA, col.ticks = 1)
axis(2, seq(0.1, 0.95, 0.1), col = NA, col.ticks = 1)
box(which = "plot", bty = "l")
abline(v = 53.0657, col="red", lty = 2)
abline(v = superscaffold26/Mb, col="blue", lty = 3)
abline(v = superscaffold54/Mb, col="blue", lty = 3)
abline(v = superscaffold35/Mb, col="blue", lty = 3)
abline(v = superscaffold36/Mb, col="blue", lty = 3)
abline(v = superscaffold62/Mb, col="blue", lty = 3)
abline(v = superscaffold67/Mb, col="blue", lty = 3)
abline(v = superscaffold69_1/Mb, col="blue", lty = 3)
abline(v = superscaffold93/Mb, col="blue", lty = 3)
abline(v = superscaffold63/Mb, col="blue", lty = 3)
abline(v = superscaffold88/Mb, col="blue", lty = 3)
abline(v = superscaffold83/Mb, col="blue", lty = 3)
abline(v = superscaffold92/Mb, col="blue", lty = 3)

