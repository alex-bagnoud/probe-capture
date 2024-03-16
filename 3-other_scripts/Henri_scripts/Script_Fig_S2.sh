#!/bin/sh

#  Script_Fig_S2.sh
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

#Figure typing: Fig. S2a.

Values_avg <- matrix(c(

3.179626083,    6.116078497,    3.391259425,    1.43302037,    2.086387763,    2.949354487,    14.78407557,    7.858650332,    11.5467057,    10.48134228,    16.43394986,    18.22244882,    0.855507499,    0.415457419,    0.246135895,

1.308616263,    8.720021012,    4.668370225,    0.462870465,    5.883429453,    0,    15.86614427,    12.19713962,    4.135848376,    12.05277245,    13.12913278,    18.9456358,    2.431374648,    0.198070681,    0.000573956

),nrow = 2,ncol = 15,byrow = TRUE)

Values_sd <- matrix(c(
1.907785398,    2.468375967,    1.585962337,    0.759938791,    1.217293747,    1.085429017,    1.984894097,    2.256081587,    2.644964871,    1.847367375,    1.419556898,    1.635882655,    0.34711287,    0.173475441,    0.383946649,
0.609848736,    3.430733493,    2.260109803,    0.20996556,    3.312497916,    0,    6.519199618,    1.052038435,    1.779413507,    5.592749984,    0.640208738,    6.219980888,    1.53003224,    0.112656434,    0.000623349
),nrow = 2,ncol = 15,byrow = TRUE)

colnames(Values_avg) = c("nifH",    "TamoA",    "amoA",    "nxrB",    "nrfA",    "hzoA",    "napA",    "narG",    "nirK",    "nirS",    "norB",    "nosZ",    "pmoA",    "mmoX",    "mcrA")
rownames(Values_avg) = c("Original DNA ", "Captured sequences ")

setEPS()
postscript("mock-community_DNA_captured3.eps")  # create eps -file
#png(file = "mock-community_DNA_captured.png")   # Give the chart file a name.
mock_community <- barplot(Values_avg, beside=TRUE,
    ylim=c(-0.2000,30.0), main="mock-community comparison",
    xlab="Sites", ylab = "Relative abundance (%)", legend.text = c("Original DNA ", "Captured sequences "),
    args.legend = list(x = "topleft"))
abline(h=0)

superpose.eb(mock_community, Values_avg, Values_sd/sqrt(6), col="black", lwd=1.5)
dev.off()  # Save the file.

#Fig. S2b :
#

setwd("/Users/siljanen/Documents/Projektit_Kasikirjoitukset/Probe_capture2/MetCap-N_1st_Validation_Outcome_data+results+Manuscript")

Values_avg <- matrix(c(1.06679986, 1.17121851, 2.27958349, 3.93658389, 6.29969237, 4.76221198, 5.63762848, 7.51003856, 9.64507717, 8.62094584, 4.54789559, 2.01688518, 2.88328791, 3.55307684, 5.70398209, 5.59583918, 1.89540013, 1.43131383, 0.86538508, 2.99702931, 1.333217, 1.076501, 1.85024986, 0.72032239, 4.39610722, 3.26583111, 2.20922934, 1.13953859, 1.31056039, 0.69902583, 5.24031281, 2.82723607, 3.26789521, 3.53118631, 1.80202276, 1.70837654, 18.7861081, 16.6724583, 13.6088024, 12.3152289, 14.1066947, 16.0365435, 6.04732888, 5.99492141, 7.03992275, 7.17288745, 10.7400191, 11.4199652, 8.32467573, 9.95600043, 14.7963899, 16.0517008, 12.00208, 10.2688786, 13.5493849, 13.0090566, 9.76056127, 7.89216656, 9.64519453, 11.0426691, 16.9600904, 16.4267259, 15.639848, 15.8806442, 17.2821479, 19.3816492, 15.2234876, 14.4205492, 13.2061362, 15.3363579, 16.0579025, 19.7798146, 0.84255393, 1.5147575, 1.03768391, 0.89210745, 0.68565609, 0.33171107, 0.138299, 0.60878337, 0.33969213, 0.46603676, 0.65517352, 0.34848996, 0.03855017, 0.07231693, 0.13197911, 0.0922752, 1.11931059, 0.05214314, 1.47262408, 2.29607319, 1.65691349, 1.2473477, 0.77622069, 0.40251842, 10.8253599, 14.7125234, 8.74373252, 5.70273476, 8.18601925, 4.14975621, 3.51455128, 5.15520691, 7.9067164, 7.11138625, 2.32405057, 1.99830995, 0.35227986, 0.85281551, 0.56030075, 0.29795519, 0.50027409, 0.21359739, 10.9627373, 9.09489373, 6.62371215, 3.37973968, 3.51508372, 1.72441017, 0, 0, 0, 0, 0, 0, 27.7516415, 15.5568588, 15.3851265, 19.0649579, 9.43186982, 8.0064111, 13.0232799, 12.6887341, 11.5871321, 10.1119675, 12.731254, 13.04047, 1.81901298, 2.15068485, 3.61464373, 4.69293209, 6.20812919, 6.32968742, 4.3579022, 7.83727949, 10.2784498, 11.6522297, 17.0722437, 21.1185299, 12.6442551, 13.1471722, 12.9024592, 12.2050574, 13.8665772, 14.0092755, 11.0316616, 12.3530025, 17.5140441, 19.6944994, 24.5319507, 28.5486565, 2.03179473, 3.75741501, 2.96863667, 4.65898583, 0.76469193, 0.40672372, 0.21289957, 0.39611426, 0.25756941, 0.17863134, 0.0916351, 0.05157441, 0, 0.00122604, 0.00056312, 0.00157523, 0, 7.9345E-05
),nrow = 2,ncol = 90,byrow = TRUE)

Values_avg_sided <- matrix(c(1.06679986,1.17121851,2.27958349,3.93658389,6.29969237,4.76221198,1.47262408,2.29607319,1.65691349,1.2473477,0.77622069,0.40251842,5.63762848,7.51003856,9.64507717,8.62094584,4.54789559,2.01688518,10.8253599,14.7125234,8.74373252,5.70273476,8.18601925,4.14975621,2.88328791,3.55307684,5.70398209,5.59583918,1.89540013,1.43131383,3.51455128,5.15520691,7.9067164,7.11138625,2.32405057,1.99830995,0.86538508,2.99702931,1.333217,1.076501,1.85024986,0.72032239,0.35227986,0.85281551,0.56030075,0.29795519,0.50027409,0.21359739,4.39610722,3.26583111,2.20922934,1.13953859,1.31056039,0.69902583,10.9627373,9.09489373,6.62371215,3.37973968,3.51508372,1.72441017,5.24031281,2.82723607,3.26789521,3.53118631,1.80202276,1.70837654,0,0,0,0,0,0,18.7861081,16.6724583,13.6088024,12.3152289,14.1066947,16.0365435,27.7516415,15.5568588,15.3851265,19.0649579,9.43186982,8.0064111,6.04732888,5.99492141,7.03992275,7.17288745,10.7400191,11.4199652,13.0232799,12.6887341,11.5871321,10.1119675,12.731254,13.04047,8.32467573,9.95600043,14.7963899,16.0517008,12.00208,10.2688786,1.81901298,2.15068485,3.61464373,4.69293209,6.20812919,6.32968742,13.5493849,13.0090566,9.76056127,7.89216656,9.64519453,11.0426691,4.3579022,7.83727949,10.2784498,11.6522297,17.0722437,21.1185299,16.9600904,16.4267259,15.639848,15.8806442,17.2821479,19.3816492,12.6442551,13.1471722,12.9024592,12.2050574,13.8665772,14.0092755,10.9627373,9.09489373,6.62371215,3.37973968,3.51508372,1.72441017,11.0316616,12.3530025,17.5140441,19.6944994,24.5319507,28.5486565,0.84255393,1.5147575,1.03768391,0.89210745,0.68565609,0.33171107,2.03179473,3.75741501,2.96863667,4.65898583,0.76469193,0.40672372,0.138299,0.60878337,0.33969213,0.46603676,0.65517352,0.34848996,0.21289957,0.39611426,0.25756941,0.17863134,0.0916351,0.05157441,0.03855017,0.07231693,0.13197911,0.0922752,1.11931059,0.05214314,0,0.00122604,0.00056312,0.00157523,0,7.9345E-05
),nrow = 30,ncol = 6,byrow = TRUE)


Values_avg_sided2 <- matrix(c(
0.99949425, 1.14431643, 2.189446174, 3.802846374, 6.203749154, 4.737904116, 1.472624084, 2.296073187, 1.656913492, 1.247347703, 0.776220692, 0.402518418,
5.28194411, 7.337538172, 9.263699867, 8.328066552, 4.478631934, 2.006590345, 10.82535991, 14.71252342, 8.743732515, 5.702734761, 8.186019251, 4.149756212,
2.701378008, 3.471465124, 5.478440156, 5.405731806, 1.866533516, 1.424007938, 3.514551279, 5.155206905, 7.906716401, 7.111386252, 2.324050569, 1.998309946,
0.810786955, 2.92818962, 1.280500088, 1.039929048, 1.822070879, 0.716645631, 0.352279863, 0.852815507, 0.560300749, 0.297955191, 0.500274087, 0.213597394,
4.11875183, 3.190817228, 2.121873906, 1.100825057, 1.290600784, 0.695457775, 10.96273727, 9.094893733, 6.623712153, 3.379739677, 3.515083717, 1.724410167,
4.909695539, 2.762296409, 3.138678921, 3.411221357, 1.774578267, 1.699656431, 0, 0, 0, 0, 0, 0,
17.60087128, 16.28950343, 13.0706949, 11.89684377, 13.89185222, 15.95468783, 27.75164152, 15.55685879, 15.38512652, 19.06495785, 9.431869825, 8.006411096,
5.665796046, 5.85722221, 6.761556205, 6.929203034, 10.57645051, 11.36167398, 13.02327993, 12.68873409, 11.58713212, 10.11196751, 12.73125404, 13.04047004,
7.799462502, 9.727317985, 14.21132386, 15.50637658, 11.8192904, 10.21646289, 1.81901298, 2.150684848, 3.614643728, 4.69293209, 6.208129194, 6.329687419,
12.69453881, 12.71024761, 9.374617588, 7.624046085, 9.498299896, 10.98630372, 4.357902203, 7.837279491, 10.27844976, 11.65222966, 17.07224371, 21.11852989,
15.89005905, 16.0494154, 15.02143064, 15.34113125, 17.01894377, 19.28271903, 12.64425507, 13.14717223, 12.90245924, 12.20505744, 13.86657721, 14.00927546,
20.57213382, 16.38624992, 16.63806417, 18.21263476, 17.33632598, 20.18928426, 11.03166159, 12.35300248, 17.51404412, 19.69449945, 24.53195067, 28.54865649,
0.789396244, 1.47996457, 0.996652707, 0.861799891, 0.675213665, 0.330017913, 2.03179473, 3.75741501, 2.96863667, 4.65898583, 0.76469193, 0.40672372,
0.129573561, 0.594800041, 0.326260319, 0.450204099, 0.645195338, 0.346711156, 0.21289957, 0.396114259, 0.25756941, 0.178631344, 0.091635099, 0.051574408,
0.036117994, 0.070655859, 0.126760508, 0.089140336, 1.102263689, 0.051876981, 0, 0.001226044, 0.000563116, 0.001575232, 0, 7.93452E-05
),nrow = 30,ncol = 6,byrow = TRUE)



rownames(Values_avg_sided2) = c("nifH",     "nifH_C",
"TamoA",    "TamoA_C",
"amoA",    "amoA_C",
"nxrB",    "nxrB_C",
"nrfA",    "nrfA_C",
"hzoA",    "hzoA_C",
"napA", "napA_C",
"narG", "narG_C",
"nirK",    "nirK_C",
"nirS",    "nirS_C",
"norB",    "norB_C",
"nosZ",    "nosZ_C",
"pmoA",    "pmoA_C",
"mmoX",    "mmoX_C",
"mcrA", "mcrA_C")

colnames(Values_avg_sided2) = c("GC 48%", "GC 50%", "GC 53%", "GC 58%", "GC 60%", "GC 63%")


setEPS()
postscript("mock_community_DNA_captured_sided2.eps")  # create eps -file
#png(file = "mock_community_DNA_captured_sided2.png")   # Give the chart file a name.
mock_community_DNA_captured_sided2 <- barplot(t(Values_avg_sided2), beside=TRUE,
    ylim=c(-0.2000,30.0), main="mock community comparison",
    xlab="Functional Gene", ylab = "Relative abundance (%)", legend.text = c("GC 48%", "GC 50%", "GC 53%", "GC 58%", "GC 60%", "GC 63%"),
    args.legend = list(x = "topleft"))
abline(h=0)
dev.off()  # Save the file.



#nirK comparison; Original DNA, Origical DNA without nitrifiers and Capruted nirK (GC% compiled together)  Inserted figure in Fig. S2a

Values_avg_nirKsepS <- matrix(c(
11.547,    4.703,    4.128),nrow = 1,ncol = 3,byrow = TRUE)

colnames(Values_avg_nirKsepS) = c("Orig. nirK with nitrif", "Orig. nirK no nitrif", "Capture nirK ")

Values_sd <- matrix(c(2.644964871,    2.54935184,    1.781201997
),nrow = 1,ncol = 3,byrow = TRUE)


setEPS()
postscript("mock_community_DNA_nirKsepS.eps")  # create eps -file
#png(file = "mock_community_DNA_nirKsepS.png")   # Give the chart file a name.
mock_community_DNA_nirKsepS <- barplot((Values_avg_nirKsepS), beside=TRUE,
    ylim=c(-0.2000,20.0), main="nirK comparison",
    xlab="Functional Gene", ylab = "Relative abundance (%)")
abline(h=0)

superpose.eb(mock_community_DNA_nirKsepS, Values_avg_nirKsepS, Values_sd, col="black", lwd=1.0)
dev.off()  # Save the file.

#nirK comparison; Original DNA, Origincal DNA without nitrifiers and Capruted nirK (GC% separated)  Inserted figure in Fig. S2b

Values_avg_nirKsep <- matrix(c(
7.799,    9.727317985,    14.21132386,    15.50637658,    11.8192904,    10.21646289,
1.848,    1.382,    3.67,    6.217,    7.2,    7.9,
1.8,    2.15,    3.6,    4.7,    6.2,    6.32
),nrow = 3,ncol = 6,byrow = TRUE)

rownames(Values_avg_nirKsep) = c("Orig. nirK with nitrif", "Orig. nirK no nitrif", "Capture nirK ")

colnames(Values_avg_nirKsep) = c("O 48%", "O 50%", "O 53%", "O 58%", "O 60%", "O 63%")

setEPS()
postscript("mock_community_DNA_nirKsepC.eps")  # create eps -file
#png(file = "mock_community_DNA_nirKsepC.png")   # Give the chart file a name.
mock_community_DNA_nirKsepC <- barplot(t(Values_avg_nirKsep), beside=TRUE,
    ylim=c(-0.2000,20.0), main="mock community comparison",
    xlab="Functional Gene", ylab = "Relative abundance (%)", legend.text = c("GC 48%", "GC 50%", "GC 53%", "GC 58%", "GC 60%", "GC 63%"),
    args.legend = list(x = "topright"))
abline(h=0)
dev.off()  # Save the file.

###########################
# Anova Fig. S2a

#anova
#Fig.S2
mock_community_comparison_simple <-read.csv("mock_community_comparison.csv",header=T, sep=";")
mock_community_comparison <-read.csv("mock_community_comparison2.csv",header=T, sep=";")

##ANOVA + tukeyHSD for relative abundance conc in the inhibitor experiment:

mock_community_comparison.aov <- aov(RelativeAbun~factor(DNA)*factor(gene), mock_community_comparison)
summary(mock_community_comparison.aov)
mock_community_comparison2.aov <- aov(RelativeAbun~factor(DNA), mock_community_comparison)
summary(mock_community_comparison2.aov)

                          Df Sum Sq Mean Sq F value   Pr(>F)
factor(DNA)                1      0     0.0   0.000        1
factor(gene)              14   5994   428.2  55.951  < 2e-16 ***
factor(DNA):factor(gene)  14    426    30.4   3.974 8.47e-06 ***
Residuals                150   1148     7.7
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


TukeyHSD(mock_community_comparison.aov, ordered=FALSE)

$`factor(DNA):factor(gene)`
                                       diff           lwr          upr     p adj
Original:TamoA-Captured:TamoA -2.390276e+00  -8.496733140   3.71618139 0.9996997
Original:amoA-Captured:amoA   -1.157887e+00  -7.264344162   4.94857037 1.0000000
Original:hzoA-Captured:hzoA    3.062838e+00  -3.043618984   9.16929555 0.9852004
Original:mcrA-Captured:mcrA    2.505219e-01  -5.855935366   6.35697917 1.0000000
Original:mmoX-Captured:mmoX    2.280084e-01  -5.878448823   6.33446571 1.0000000
Original:napA-Captured:napA   -6.118383e-01  -6.718295546   5.49461899 1.0000000
Original:narG-Captured:narG   -4.127965e+00 -10.234422757   1.97849177 0.7153201
Original:nifH-Captured:nifH    1.944065e+00  -4.162391844   8.05052269 0.9999949
Original:nirK-Captured:nirK    7.764106e+00   1.657648595  13.87056313 0.0010815
Original:nirS-Captured:nirS   -1.236267e+00  -7.342724233   4.87019030 1.0000000
Original:norB-Captured:norB    3.799385e+00  -2.307072449   9.90584208 0.8470638
Original:nosZ-Captured:nosZ   -3.274928e+00  -9.381385096   2.83152944 0.9658426
Original:nrfA-Captured:nrfA   -3.713381e+00  -9.819837973   2.39307656 0.8746122
Original:nxrB-Captured:nxrB    1.010914e+00  -5.095543623   7.11737091 1.0000000
Original:pmoA-Captured:pmoA   -1.547296e+00  -7.653753590   4.55916094 1.0000000

attach(mock_community_comparison_simple)
mock_community_comparison_simple.aov <- aov(Original~Captured)
summary(mock_community_comparison_simple.aov)
t.test(RelativeAbun~factor(DNA), mock_community_comparison)

Welch Two Sample t-test

data:  RelativeAbun by factor(DNA)
t = 1.8518e-10, df = 174.83, p-value = 1
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -1.91838  1.91838
sample estimates:
mean in group Captured mean in group Original
              6.666667               6.666667

