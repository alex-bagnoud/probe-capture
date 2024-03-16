#!/bin/sh

#  Script_Fig_3.sh
#  
#
#  Created by Henri Siljanen2 on 9.2.2024.
#  

#superpose function
superpose.eb <-
function (x, y, ebl, ebu = ebl, length = 0.04, ...)
    arrows(x, y + ebu, x, y - ebl, angle = 90, code = 3,
    length = length, ...)



setwd("/Users/siljanen2/Documents/Projektit_Kasikirjoitukset/Probe_capture2_MG-documents/MetCap-N_1st_Validation_Outcome_data+results+Manuscript")

#Fig. 3
#Belfontaine:
Reads_avg <- matrix(c(
617.333333, 7.33333333, 43, 908, 273.333333, 0.33333333, 897, 1279.33333, 108.333333, 557.333333, 613, 527.666667, 50.3333333, 29.3333333, 179.666667,
30803.6667,353,3004,12861.3333,1068.66667,2.66666667,14391.6667,21687.3333,2578,23747.3333,1856.33333,9156.33333,2902,301.333333,4159
),nrow = 2,ncol = 15,byrow = TRUE)

Values_sd <- matrix(c(
97.6057831, 5.73488351, 22.9056034, 317.754622, 5.79271573, 0.47140452, 195.545391, 123.321621, 25.9529489, 12.1197726, 88.0567999, 150.380259, 21.5612822, 13.2245563, 130.527988,
8261.49338,234.355855,1348.43045,3643.93106,366.147815,2.05480467,3124.73842,5432.02365,737.130925,3665.31239,489.10417,3387.29059,1284.95292,62.8985073,3274.41384
),nrow = 2,ncol = 15,byrow = TRUE)

colnames(Reads_avg) = c("nifH",    "TamoA",    "amoA",    "nxrB",    "nrfA",    "hzoA",    "napA",    "narG",    "nirK",    "nirS",    "norB",    "nosZ",    "pmoA",    "mmoX",    "mcrA")
rownames(Reads_avg) = c("Shotgun ", "Captured  ")

setEPS()
postscript("belfointaine_N_reads_fig4.eps")  # create eps -file
#png(file = "belfointaine_N_reads_fig4.png")   # Give the chart file a name.
belfontaine_fig2 <- barplot(Reads_avg, beside=FALSE,
 log= "TRUE",     ylim=c(0,100000.0), scale_y_log10(breaks=c(1,10,100,1000,10000,100000), main="Belfontaine",
    xlab="Genes", ylab = "Number of reads", legend.text = c("Shotgun metagenomics", "Captured metagenomics")
)
dev.off()  # Save the file.

superpose.eb(belfontaine_fig2, Reads_avg, Values_sd, col="black", lwd=1.5)


#Hungarian:
Values_avg <- matrix(c(
1.976184514, 4.950860694, 0.556781523,  1.100563387, 1.538320223, 0, 37.6885202, 17.92418326, 9.995615811, 1.19427028, 5.1459641, 9.203134248, 0.114219414, 8.610971834, 0.000410519,
0.65870656, 2.57957066, 0.287746176,  6.2148364, 6.54494727, 0, 26.3627441, 26.2020598, 8.99831344, 1.089098868, 4.741509139, 14.34173513, 0.310447747, 1.445692189, 0),nrow = 2,ncol = 15,byrow = TRUE)

Values_sd <- matrix(c(
0.171393046, 1.074144765, 0.197682724, 0.074547684, 0.031503633, 0, 0.734114453, 0.087407806, 0.31870049, 0.151706566, 0.195100733, 0.55277326, 0.01130385, 1.057703178, 0.00071104,
0.342774882, 2.144038334, 0.115799971, 0.550684377, 0.610072558, 0, 1.142455965, 0.600751255, 0.920879301, 0.162092609, 0.166995678, 0.085737433, 0.063295503, 0.434609049, 0),nrow = 2,ncol = 15,byrow = TRUE)

colnames(Values_avg) = c("nifH",    "TamoA",    "amoA",    "nxrB",    "nrfA",    "hzoA",    "napA",    "narG",    "nirK",    "nirS",    "norB",    "nosZ",    "pmoA",    "mmoX",    "mcrA")
rownames(Reads_avg) = c("Shotgun ", "Captured  ")

setEPS()
postscript("hungarian_fig2.eps")  # create eps -file
#png(file = "hungarian_fig2.png")   # Give the chart file a name.
hungarian_fig2 <- barplot(Values_avg, beside=TRUE,
    ylim=c(-0.2000,40.0), main="Hungarian",
    xlab="Genes", ylab = "Relative abundance (%)", legend.text = c("Shotgun metagenomics", "Captured metagenomics"),
    args.legend = list(x = "topleft"))
abline(h=0)

superpose.eb(hungarian_fig2, Values_avg, Values_sd, col="black", lwd=1.5)
dev.off()  # Save the file.


######################
# ANOVA

#Fig. 3
#Belfontaine
reads_shotgun_captured_belfontaine <-read.csv("reads_shotgun_captured_belfontaine.csv",header=T, sep=";")
##ANOVA + tukeyHSD for number of reads with capture or shotgun:

reads_shotgun_captured_belfontaine.aov <- aov(number_reads~factor(Method)*factor(Gene), reads_shotgun_captured_belfontaine)
summary(reads_shotgun_captured_belfontaine.aov)

                           Df    Sum Sq   Mean Sq F value   Pr(>F)
factor(Method)               1 1.508e+09 1.508e+09  187.17  < 2e-16 ***
factor(Gene)                14 2.191e+09 1.565e+08   19.43  < 2e-16 ***
factor(Method):factor(Gene) 14 1.952e+09 1.394e+08   17.31 4.63e-16 ***
Residuals                   60 4.833e+08 8.054e+06
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


TukeyHSD(reads_shotgun_captured_belfontaine.aov, ordered=FALSE)

$`factor(Method):factor(Gene)`
                                      diff           lwr          upr     p adj
shotgun:napA-captured:napA   -13494.666667 -22615.343210  -4373.99012 0.0000907


shotgun:narG-captured:narG   -20408.000000 -29528.676544 -11287.32346 0.0000000
shotgun:nifH-captured:nifH   -30186.333333 -39307.009877 -21065.65679 0.0000000
shotgun:nirS-captured:nirS   -23190.000000 -32310.676544 -14069.32346 0.0000000
shotgun:nxrB-captured:nxrB   -11953.333333 -21074.009877  -2832.65679 0.0010095

#Hungary
reads_shotgun_captured_hungary <-read.csv("reads_shotgun_captured_hungary.csv",header=T, sep=";")
##ANOVA + tukeyHSD for number of reads with capture or shotgun:

reads_shotgun_captured_hungary.aov <- aov(number_reads~factor(Method)*factor(Gene), reads_shotgun_captured_hungary)
summary(reads_shotgun_captured_hungary.aov)

                                  Df    Sum Sq   Mean Sq F value Pr(>F)
factor(Method)               1 1.771e+09 1.771e+09  441.94 <2e-16 ***
factor(Gene)                14 4.062e+09 2.901e+08   72.38 <2e-16 ***
factor(Method):factor(Gene) 14 3.749e+09 2.678e+08   66.80 <2e-16 ***
Residuals                   60 2.405e+08 4.008e+06
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


TukeyHSD(reads_shotgun_captured_hungary.aov, ordered=FALSE)

$`factor(Method):factor(Gene)`
                                      diff         lwr          upr     p adj
shotgun:TamoA-captured:TamoA -6.727667e+03 -13161.9234   -293.40990 0.0299949

shotgun:mmoX-captured:mmoX   -1.160833e+04 -18042.5901  -5174.07657 0.0000007
shotgun:napA-captured:napA   -5.059400e+04 -57028.2568 -44159.74324 0.0000000
shotgun:narG-captured:narG   -2.357333e+04 -30007.5901 -17139.07657 0.0000000

shotgun:nirK-captured:nirK   -1.333500e+04 -19769.2568  -6900.74324 0.0000000
shotgun:nirS-captured:nirS   -1.595667e+03  -8029.9234   4838.59010 0.9999999
shotgun:norB-captured:norB   -6.870000e+03 -13304.2568   -435.74324 0.0231908
shotgun:nosZ-captured:nosZ   -1.199433e+04 -18428.5901  -5560.07657 0.0000003
shotgun:nrfA-captured:nrfA   -1.878333e+03  -8312.5901   4555.92343 0.9999964
shotgun:nxrB-captured:nxrB   -1.294667e+03  -7728.9234   5139.59010 1.0000000
shotgun:pmoA-captured:pmoA   -1.456667e+02  -6579.9234   6288.59010 1.0000000




