#!/bin/sh

#  Script_Fig_3.sh
#  
#
# Created by Henri Siljanen on 28.9.2024, after pci-microbiology revision


setwd("/Users/siljanen/Library/Mobile\ Documents/com\~apple\~CloudDocs/Documents/Probe_capture2_MG-documents/MetCap-N_1st_Validation_Outcome_data+results+Manuscript/")

#superpose function
superpose.eb <-
function (x, y, ebl, ebu = ebl, length = 0.04, ...)
    arrows(x, y + ebu, x, y - ebl, angle = 90, code = 3,
    length = length, ...)

#Belfontaine:
Total_read_RA <- matrix(c(
0.001911033, 1.54089E-05, 0.000102607, 0.002340381, 0.000821016, 7.64167E-07, 0.00244101, 0.003668652, 0.000292022, 0.001648217, 0.001723354, 0.001398577, 0.000128771, 0.000103205, 0.000689391, 
13.34051243, 0.139450074, 0.733383736, 2.59240322, 0.58008041, 0.00130378, 10.03299, 9.35380698, 1.454513908, 10.255586, 1.070992322, 4.059711994, 0.981120149, 0.171831185, 1.898530578
),nrow = 2,ncol = 15,byrow = TRUE)

Values_sd <- matrix(c(
0.000997664, 1.15718E-05, 3.02979E-05, 0.000274423, 0.000346357, 1.0807E-06, 0.000335242, 0.001070409, 3.16418E-05, 0.000621598, 0.000404019, 9.31457E-05, 2.22262E-05, 8.52551E-05, 0.000711681, 
3.949794036, 0.078180128, 0.49125469, 2.894437223, 0.123749205, 0.001029063, 2.504484603, 1.117919135, 0.323352616, 0.327699814, 0.355213166, 1.107116742, 0.176496915, 0.061601251, 1.526503534
),nrow = 2,ncol = 15,byrow = TRUE)

colnames(Total_read_RA ) = c("nifH",	"TamoA",	"amoA",	"nxrB",	"nrfA",	"hzoA",	"napA",	"narG",	"nirK",	"nirS",	"norB",	"nosZ",	"pmoA",	"mmoX",	"mcrA")
rownames(Total_read_RA ) = c("Shotgun ", "Captured  ")

setEPS()
postscript("belfointaine_N_reads_fig3pci.eps")  # create eps -file
#png(file = "belfointaine_N_reads_fig3pci.png")   # Give the chart file a name.
belfontaine_fig3 <- barplot(Total_read_RA , beside=TRUE,  
 ylim=c(0,20), main="Belfontaine",
    xlab="Genes", ylab = "Relative abundance, total reads (%)", legend.text = c("Shotgun metagenomics", "Captured metagenomics")
)
dev.off()  # Save the file.

setEPS()
postscript("belfointaine_N_reads_fig3pci.eps")  # create eps -file
#png(file = "belfointaine_N_reads_fig3pci.png")   # Give the chart file a name.
belfontaine_fig3 <- barplot(Total_read_RA, beside=TRUE,  
 ylim=c(-1.0,30.0), main="Belfontaine",
    xlab="Genes", ylab = "Relative abundance, total reads (%)", legend.text = c("Shotgun metagenomics", "Captured metagenomics"),
    args.legend = list(x = "topleft"))

superpose.eb(belfontaine_fig3, Total_read_RA , Values_sd, col="black", lwd=1.5)
dev.off()  # Save the file.


setEPS()
postscript("belfointaine_N_reads_fig3pci_s.eps")  # create eps -file
#png(file = "belfointaine_N_reads_fig3pci.png")   # Give the chart file a name.
belfontaine_fig3_s <- barplot(Total_read_RA, beside=TRUE,  
 ylim=c(-0.005,0.005), main="Belfontaine",
    xlab="Genes", ylab = "Relative abundance, total reads (%)")
abline(h=0)

superpose.eb(belfontaine_fig3_s, Total_read_RA , Values_sd, col="black", lwd=1.5)
dev.off()  # Save the file.





#Hungarian:
Total_read_RA_H <- matrix(c(
7.34066E-05, 0.000292496, 3.11156E-05, 0.00057338, 0.000603068, 0, 0.002462011, 0.002462252, 0.000829175, 0.00010485, 0.000442137, 0.001355456, 3.00713E-05, 0.000151054, 0,

1.042072243, 2.611045776, 0.293029832, 0.580382904, 0.811505869, 0, 19.88344074, 9.456528627, 5.273714839, 0.62998497, 2.714774935, 4.857191124, 0.060250088, 4.54600597, 0.000213994

),nrow = 2,ncol = 15,byrow = TRUE)

Values_sd_H <- matrix(c(
6.30724E-05, 0.000239411, 2.31809E-05, 0.000236922, 0.000258183, 0, 0.001115166, 0.001131136, 0.000356943, 5.68066E-05, 0.000197648, 0.000644432, 1.70702E-05, 0.000107277, 0,

0.065579758, 0.460207367, 0.082249559, 0.02696476, 0.006875067, 0, 0.295429452, 0.049540031, 0.145910308, 0.064743462, 0.079618242, 0.272345285, 0.004760226, 0.483157914, 0.000302633
),nrow = 2,ncol = 15,byrow = TRUE)

colnames(Total_read_RA_H) = c("nifH",	"TamoA",	"amoA",	"nxrB",	"nrfA",	"hzoA",	"napA",	"narG",	"nirK",	"nirS",	"norB",	"nosZ",	"pmoA",	"mmoX",	"mcrA")
rownames(Total_read_RA_H) = c("Shotgun ", "Captured  ")

setEPS()
postscript("hungarian_RA_reads_fig3pci.eps")  # create eps -file
#png(file = "hungarian_RA_reads_fig3pci.png")   # Give the chart file a name.
hungarian_fig3 <- barplot(Total_read_RA_H, beside=TRUE,  
    ylim=c(-1.0,25.0), main="Hungarian",
    xlab="Genes", ylab = "Relative abundance, total reads (%)", legend.text = c("Shotgun metagenomics", "Captured metagenomics"),
    args.legend = list(x = "topleft"))

superpose.eb(hungarian_fig3, Total_read_RA_H, Values_sd_H, col="black", lwd=1.5)
dev.off()  # Save the file.

setEPS()
postscript("hungarian_RA_reads_fig3pci_s.eps")  # create eps -file
#png(file = "hungarian_RA_reads_fig3pci_s.png")   # Give the chart file a name.
hungarian_fig3s <- barplot(Total_read_RA_H, beside=TRUE,  
    ylim=c(-0.005,0.005), main="Hungarian",
    xlab="Genes", ylab = "Relative abundance, total reads (%)" )
abline(h=0)

superpose.eb(hungarian_fig3s, Total_read_RA_H, Values_sd_H, col="black", lwd=1.5)
dev.off()  # Save the file.


#Fig. 3_Revision_after_PCI_microbiology


getwd()

setwd("/Users/siljanen/Library/Mobile Documents/com~apple~CloudDocs/Documents/Probe_capture2_MG-documents/MetCap-N_1st_Validation_Outcome_data+results+Manuscript")

#Belfontaine
Total_reads_RA_shotgun_captured_belfontaine <-read.csv("total_reads_RA_shotgun_captured_belfontaine.csv",header=T, sep=";")
##ANOVA + tukeyHSD for number of reads with capture or shotgun:

Total_reads_RA_shotgun_captured_belfontaine.aov <- aov(Total_read_RA~factor(Method)*factor(Gene), Total_reads_RA_shotgun_captured_belfontaine)
summary(Total_reads_RA_shotgun_captured_belfontaine.aov)

         Df Sum Sq Mean Sq F value   Pr(>F)    
factor(Method)               1  320.9   320.9  179.82  < 2e-16 ***
factor(Gene)                14  433.9    31.0   17.37 4.27e-16 ***
factor(Method):factor(Gene) 14  433.6    31.0   17.36 4.34e-16 ***
Residuals                   60  107.1     1.8                     
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
> 


TukeyHSD(Total_reads_RA_shotgun_captured_belfontaine.aov, ordered=FALSE)

$`factor(Method):factor(Gene)`
                                      diff         lwr        upr     p adj
shotgun:TamoA-captured:TamoA -1.394347e-01  -4.4326442  4.1537749 1.0000000
shotgun:amoA-captured:amoA   -7.332811e-01  -5.0264907  3.5599284 1.0000000
shotgun:hzoA-captured:hzoA   -1.303016e-03  -4.2945126  4.2919065 1.0000000
shotgun:mcrA-captured:mcrA   -1.897841e+00  -6.1910507  2.3953684 0.9945745
shotgun:mmoX-captured:mmoX   -1.717280e-01  -4.4649375  4.1214816 1.0000000
shotgun:napA-captured:napA   -1.003055e+01 -14.3237585 -5.7373395 0.0000000
shotgun:narG-captured:narG   -9.350138e+00 -13.6433479 -5.0569288 0.0000000
shotgun:nifH-captured:nifH   -1.333860e+01 -17.6318109 -9.0453919 0.0000000
shotgun:nirK-captured:nirK   -1.454222e+00  -5.7474314  2.8389877 0.9999322
shotgun:nirS-captured:nirS   -1.025394e+01 -14.5471473 -5.9607282 0.0000000
shotgun:norB-captured:norB   -1.069269e+00  -5.3624785  3.2239406 0.9999999
shotgun:nosZ-captured:nosZ   -4.058313e+00  -8.3515230  0.2348961 0.0888384
shotgun:nrfA-captured:nrfA   -5.792594e-01  -4.8724689  3.7139501 1.0000000
shotgun:nxrB-captured:nxrB   -2.590063e+00  -6.8832724  1.7031467 0.8382351
shotgun:pmoA-captured:pmoA   -9.809914e-01  -5.2742009  3.3122182 1.0000000





#Hungary
Total_reads_RA_shotgun_captured_hungary <-read.csv("total_reads_RA_shotgun_captured_hungary.csv",header=T, sep=";")
##ANOVA + tukeyHSD for number of reads with capture or shotgun:

Total_reads_RA_shotgun_captured_hungary.aov <- aov(number_reads~factor(Method)*factor(Gene), Total_reads_RA_shotgun_captured_hungary)
summary(Total_reads_RA_shotgun_captured_hungary.aov)

                    Df Sum Sq Mean Sq F value Pr(>F)    
factor(Method)               1  278.3  278.26    8525 <2e-16 ***
factor(Gene)                14  582.2   41.59    1274 <2e-16 ***
factor(Method):factor(Gene) 14  581.9   41.56    1273 <2e-16 ***
Residuals                   60    2.0    0.03                   
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


TukeyHSD(Total_reads_RA_shotgun_captured_hungary.aov, ordered=FALSE) 

$`factor(Method):factor(Gene)`
                                      diff           lwr           upr     p adj
shotgun:TamoA-captured:TamoA -2.610753e+00 -3.191372e+00 -2.030134e+00 0.0000000
shotgun:amoA-captured:amoA   -2.929987e-01 -8.736176e-01  2.876202e-01 0.9707278
shotgun:hzoA-captured:hzoA    6.217249e-15 -5.806189e-01  5.806189e-01 1.0000000
shotgun:mcrA-captured:mcrA   -2.139940e-04 -5.808329e-01  5.804049e-01 1.0000000
shotgun:mmoX-captured:mmoX   -4.545855e+00 -5.126474e+00 -3.965236e+00 0.0000000
shotgun:napA-captured:napA   -1.988098e+01 -2.046160e+01 -1.930036e+01 0.0000000
shotgun:narG-captured:narG   -9.454066e+00 -1.003469e+01 -8.873448e+00 0.0000000
shotgun:nifH-captured:nifH   -1.041999e+00 -1.622618e+00 -4.613800e-01 0.0000008
shotgun:nirK-captured:nirK   -5.272886e+00 -5.853505e+00 -4.692267e+00 0.0000000
shotgun:nirS-captured:nirS   -6.298801e-01 -1.210499e+00 -4.926125e-02 0.0189298
shotgun:norB-captured:norB   -2.714333e+00 -3.294952e+00 -2.133714e+00 0.0000000
shotgun:nosZ-captured:nosZ   -4.855836e+00 -5.436455e+00 -4.275217e+00 0.0000000
shotgun:nrfA-captured:nrfA   -8.109028e-01 -1.391522e+00 -2.302839e-01 0.0003009
shotgun:nxrB-captured:nxrB   -5.798095e-01 -1.160428e+00  8.093489e-04 0.0507650
shotgun:pmoA-captured:pmoA   -6.022002e-02 -6.408389e-01  5.203989e-01 1.0000000




#  Created by Henri Siljanen2 on 9.2.2024. (earlier version)
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




