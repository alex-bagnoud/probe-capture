


#  Script_Fig_2.sh
#  
#
#  Created by Henri Siljanen2 on 9.2.2024.
#  

# R-script

#superpose function
superpose.eb <-
function (x, y, ebl, ebu = ebl, length = 0.04, ...)
    arrows(x, y + ebu, x, y - ebl, angle = 90, code = 3,
    length = length, ...)



setwd("/Users/siljanen2/Documents/Projektit_Kasikirjoitukset/Probe_capture2_MG-documents/MetCap-N_1st_Validation_Outcome_data+results+Manuscript")

#Fig. 2
#Belfontaine:
Values_avg <- matrix(c(10.3441966, 0.109424147, 0.664492183, 14.39209094, 4.539866608, 0.00518081, 14.43998378, 20.9303173, 1.738539335, 9.211585595, 9.965618453, 8.419408721, 0.793126123, 0.514654926, 3.253568927,
23.42597185, 0.248137156, 1.312049869, 4.668767264, 1.024495521, 0.002316736, 17.65993275, 16.54639777, 2.571267121, 18.11658892, 1.878197085, 7.200984347, 1.738314019, 0.302832082, 3.303747508
),nrow = 2,ncol = 15,byrow = TRUE)

Values_sd <- matrix(c(3.293084462, 0.1027814, 0.373856334, 4.812941784, 0.843358131, 0.008973427, 1.984129157, 0.735582029, 0.283538461, 1.253614875, 0.284579262, 1.845266111, 0.321203189, 0.369892126, 3.362029308,
7.919045244, 0.170961233, 1.089212775, 6.423335833, 0.270439206, 0.002275213, 5.243197789, 2.731354432, 0.71512939, 1.144139515, 0.721808655, 2.51656519, 0.424580645, 0.13258371, 3.186479547
),nrow = 2,ncol = 15,byrow = TRUE)

colnames(Values_avg) = c("nifH",    "TamoA",    "amoA",    "nxrB",    "nrfA",    "hzoA",    "napA",    "narG",    "nirK",    "nirS",    "norB",    "nosZ",    "pmoA",    "mmoX",    "mcrA")
rownames(Values_avg) = c( "Shotgun metagenomics ", "Targeted metagenomics")

setEPS()
postscript("belfointaine_fig3.eps")  # create eps -file
#png(file = "belfointaine_fig3.png")   # Give the chart file a name.
belfontaine_fig3 <- barplot(Values_avg, beside=TRUE,
    ylim=c(-0.2000,40.0), main="Belfontaine",
    xlab="Genes", ylab = "Relative abundance (%)", legend.text = c( "Shotgun metagenomics ", "Targeted metagenomics "),
    args.legend = list(x = "topleft"))
abline(h=0)

superpose.eb(belfontaine_fig3, Values_avg, Values_sd, col="black", lwd=1.0)
dev.off()  # Save the file.

#Hungarian:
Values_avg <- matrix(c(
0.65870656, 2.57957066, 0.287746176,  6.2148364, 6.54494727, 0, 26.3627441, 26.2020598, 8.99831344, 1.089098868, 4.741509139, 14.34173513, 0.310447747, 1.445692189, 0,
1.976184514, 4.950860694, 0.556781523,  1.100563387, 1.538320223, 0, 37.6885202, 17.92418326, 9.995615811, 1.19427028, 5.1459641, 9.203134248, 0.114219414, 8.610971834, 0.000410519
),nrow = 2,ncol = 15,byrow = TRUE)

Values_sd <- matrix(c(
 0.342774882, 2.144038334, 0.115799971, 0.550684377, 0.610072558, 0, 1.142455965, 0.600751255, 0.920879301, 0.162092609, 0.166995678, 0.085737433, 0.063295503, 0.434609049, 0,
0.171393046, 1.074144765, 0.197682724, 0.074547684, 0.031503633, 0, 0.734114453, 0.087407806, 0.31870049, 0.151706566, 0.195100733, 0.55277326, 0.01130385, 1.057703178, 0.00071104),nrow = 2,ncol = 15,byrow = TRUE)

colnames(Values_avg) = c("nifH",    "TamoA",    "amoA",    "nxrB",    "nrfA",    "hzoA",    "napA",    "narG",    "nirK",    "nirS",    "norB",    "nosZ",    "pmoA",    "mmoX",    "mcrA")
rownames(Values_avg) = c("Shotgun metagenomics ", "targeted metagenomics ")

setEPS()
postscript("hungarian_fig3.eps")  # create eps -file
#png(file = "hungarian_fig3.png")   # Give the chart file a name.
hungarian_fig3 <- barplot(Values_avg, beside=TRUE,
    ylim=c(-0.2000,40.0), main="Hungarian",
    xlab="Genes", ylab = "Relative abundance (%)", legend.text = c("Shotgun metagenomics ", "Targeted metagenomics "),
    args.legend = list(x = "topleft"))
abline(h=0)

superpose.eb(hungarian_fig3, Values_avg, Values_sd, col="black", lwd=1.0)
dev.off()  # Save the file.


############
# ANOVA Fig 2.

#Fig.2
#BF
shotgun_BF_comparison <-read.csv("shotgun_BF_comparison2.csv",header=T, sep=";")
##ANOVA + tukeyHSD for relative abundance conc in the inhibitor experiment:

shotgun_BF_comparison.aov <- aov(RelativeAbun~factor(DNA)*factor(gene), shotgun_BF_comparison)
summary(shotgun_BF_comparison.aov)

         Df Sum Sq Mean Sq F value   Pr(>F)
factor(DNA)               1      0    0.05   0.006    0.936
factor(gene)             14   3787  270.51  37.699  < 2e-16 ***
factor(DNA):factor(gene) 14    684   48.84   6.807 4.81e-08 ***
Residuals                60    431    7.18
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


TukeyHSD(shotgun_BF_comparison.aov, ordered=FALSE)

$`factor(DNA):factor(gene)`
 $`factor(DNA):factor(gene)`
                                      diff          lwr           upr     p adj
shotgun:TamoA-captured:TamoA  -0.138713009  -8.74737205   8.469946032 1.0000000
shotgun:amoA-captured:amoA    -0.647557686  -9.25621673   7.961101356 1.0000000
shotgun:hzoA-captured:hzoA     0.002864074  -8.60579497   8.611523115 1.0000000
shotgun:mcrA-captured:mcrA    -0.050178581  -8.65883762   8.558480460 1.0000000
shotgun:mmoX-captured:mmoX     0.211822844  -8.39683620   8.820481885 1.0000000
shotgun:napA-captured:napA    -3.219948967 -11.82860801   5.388710075 0.9995994
shotgun:narG-captured:narG     4.383919527  -4.22473951  12.992578568 0.9675545
shotgun:nifH-captured:nifH   -13.081775259 -21.69043430  -4.473116217 0.0000503  ***
shotgun:nirK-captured:nirK    -0.832727786  -9.44138683   7.775931256 1.0000000
shotgun:nirS-captured:nirS    -8.905003324 -17.51366237  -0.296344283 0.0340817 *
shotgun:norB-captured:norB     8.087421367  -0.52123767  16.696080409 0.0941996
shotgun:nosZ-captured:nosZ     1.218424376  -7.39023467   9.827083417 1.0000000
shotgun:nrfA-captured:nrfA     3.515371087  -5.09328795  12.124030128 0.9983213
shotgun:nxrB-captured:nxrB     9.723323675   1.11466463  18.331982716 0.0109847 *
shotgun:pmoA-captured:pmoA    -0.945187895  -9.55384694   7.663471146 1.0000000

#HUN
shotgun_HUN_comparison <-read.csv("shotgun_HUN_comparison2.csv",header=T, sep=";")
##ANOVA + tukeyHSD for relative abundance conc in the inhibitor experiment:

shotgun_HUN_comparison.aov <- aov(RelativeAbun~factor(DNA)*factor(gene), shotgun_HUN_comparison)
summary(shotgun_HUN_comparison.aov)

                       Df Sum Sq Mean Sq  F value Pr(>F)
factor(DNA)               1      0     0.0    0.013   0.91
factor(gene)             14   7037   502.7 1307.427 <2e-16 ***
factor(DNA):factor(gene) 14    502    35.8   93.189 <2e-16 ***
Residuals                60     23     0.4
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
>


TukeyHSD(shotgun_HUN_comparison.aov, ordered=FALSE)

$`factor(DNA):factor(gene)`
                                      diff          lwr          upr     p adj
shotgun:TamoA-captured:TamoA -2.371290e+00  -4.36398728  -0.37859279 0.0050933  *

shotgun:amoA-captured:amoA   -2.690353e-01  -2.26173259   1.72366190 1.0000000

shotgun:hzoA-captured:hzoA    1.022619e-15  -1.99269725   1.99269725 1.0000000

shotgun:mcrA-captured:mcrA   -4.105190e-04  -1.99310776   1.99228673 1.0000000

shotgun:mmoX-captured:mmoX   -7.165280e+00  -9.15797689  -5.17258240 0.0000000 *

shotgun:napA-captured:napA   -1.132578e+01 -13.31847334  -9.33307885 0.0000000 *

shotgun:narG-captured:narG    8.277877e+00   6.28517930  10.27057379 0.0000000  *

shotgun:nifH-captured:nifH   -1.317478e+00  -3.31017520   0.67521929 0.6973344

shotgun:nirK-captured:nirK   -9.973024e-01  -2.98999962   0.99539487 0.9733785

shotgun:nirS-captured:nirS   -1.051714e-01  -2.09786866   1.88752583 1.0000000

shotgun:norB-captured:norB   -4.044550e-01  -2.39715221   1.58824229 1.0000000

shotgun:nosZ-captured:nosZ    5.138601e+00   3.14590364   7.13129813 0.0000000 *

shotgun:nrfA-captured:nrfA    5.006627e+00   3.01392980   6.99932429 0.0000000 *

shotgun:nxrB-captured:nxrB    5.114273e+00   3.12157577   7.10697026 0.0000000  *
shotgun:pmoA-captured:pmoA    1.962283e-01  -1.79646891   2.18892558 1.0000000


###########################################################################################
# Fig. 4. Statistical test for 
# Dropbox link to the data:
https://www.dropbox.com/scl/fi/37lml9f5fx3geknnr2f4h/4b-bubble_plot_table_kruskall_test_input2_transp.csv?rlkey=9chj3f23fgvcsubxplfisna84&dl=0

#Analysis of Wilcox comparisons at done in R 3.5.3. 

setwd("/Users/siljanen 1/Library/Mobile Documents/com~apple~CloudDocs/Documents/TamoA_Hungary_shotgun_dada_usearch/usearch_captured_shotgun_TamoA")

bubble_krus_c_t <-read.csv("4b-bubble_plot_table_kruskall_test_input2_transp.csv",header=T, sep=";")

attach(bubble_krus_c_t)
pairwise.wilcox.test(NS, group, , alternative = "less", p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3, group, , alternative = "less", p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.1, group, , alternative = "less", p.adjust.method = "BH")
#Example outcome
#Pairwise comparisons using Wilcoxon rank sum test 

#data:  NS.Alpha.3.1 and group 

#  1     2    
#  2 0.048 -    
#  3 0.048 1.000

#P value adjustment method: BH 

pairwise.wilcox.test(NS.Alpha.3.1.1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.1.1.1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.1.1.1_OTU1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.1.1.1_OTU2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.1.1.1_OTU3, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.1.1.2_OTU3, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.1.1.3, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.1.1.Incertae_sedis_OTU1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.1.1_OTU4, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.1.2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.1.2_OTU1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.1.2_OTU2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.2.1.1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.2.1.1.1.2_OTU2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.2.2_OTU1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.2.2_OTU7, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.2.3, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.2.3.1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.2.3.1.3_OTU1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.2.3.1.4, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.2.3.1.4_OTU3, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.2.3.1.7, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.3.1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.3.1_OTU1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Alpha.3.3.1_OTU3, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Beta.1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Beta.1_OTU10, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Beta.1_OTU4, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Beta.1_OTU6, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Beta.1_OTU9, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Beta.2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Beta.2.1_OTU2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Beta.2.2_OTU2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Beta.2_OTU1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.1.1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.1.1.2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.1.1.2.1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.1.1.2.1_OTU12, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.1.1.2.1_OTU13, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.1.1.2.2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.1.1.2.2_OTU2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.1.1.2_OTU2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.1.1.3, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.1.1.3_OTU1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.1.1.3_OTU2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.1.Incertae_sedis.3_OTU2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.1_OTU9, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.2.1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.2.1_OTU1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.2.1_OTU2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.2.2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.2.2.1_OTU2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.2.2.1_OTU4, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Delta.2.2.2_OTU3, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Epsilon.2.1_OTU5, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Epsilon.2.2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Epsilon.2.2_OTU5, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.1.1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.1.1_OTU2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.1.2_OTU1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.1.2_OTU3, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.1.Incertae_sedis_OTU2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.1.1_OTU2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.1.1_OTU6, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.1.2.2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.1.2.2_OTU3, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.1.2.2_OTU4, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.1_OTU1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.2.1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.2.1.1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.2.1_OTU3, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.2.1_OTU4, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.2.1_OTU5, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.2.2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.2.2_OTU6, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.2.3_OTU1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.2.3_OTU3, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.2.3_OTU4, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.3.2.3_OTU3, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.3.2_OTU2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.3.2_OTU5, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.3.Incertae_sedis_OTU2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.3_OTU1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Gamma.2.3_OTU2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Zeta, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Zeta.1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Zeta.1.1_OTU1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Zeta.1.1_OTU3, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Zeta.1.2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Zeta.1.2_OTU4, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Zeta.1.2_OTU8, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Zeta.2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Zeta_OTU2, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS.Zeta_OTU3, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(NS_OTU1, group, , alternative = "less",  p.adjust.method = "BH")
pairwise.wilcox.test(unknown, group, , alternative = "less",  p.adjust.method = "BH")



