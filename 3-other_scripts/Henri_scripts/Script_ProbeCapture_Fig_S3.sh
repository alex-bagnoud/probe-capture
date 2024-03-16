
#  Script_ProbeCapture_Fig_S3.sh
#  
#
#  Created by Henri Siljanen2 on 9.2.2024.
#  

# To plot the figures in R.

### 9.7.2019 / make it again for all together 18.11.2021
##Plots for different GC%'s.

setwd("/Users/siljanen1/Documents/Probe_capture2/MetCap-N_1st_Validation_Outcome_data+results")

Data_t <-read.csv("correlation_table_cultures_transpose.csv",header=T,sep=";",row.names=1)

Data_all <-read.csv("mock_community_comparison4_for_scatter.csv",header=T,sep=";",row.names=1)

##
# _DNA = Original DNA relative abudnance of genes, which was taken in to the Mock-community by pipetting this abundance too pooled DNA.
# _bx = Targted probe capture relative abundance, after sequencing and after compiling numbers of reads with hmmer-profile search.

#GC47%
ggscatter(Data_t, x = "X47_DNA", y = "X47_bx",
          
          add = "reg.line", conf.int = TRUE,
          add.params = list(fill = "lightgray"),
          ggtheme = theme_minimal()
          )+
  stat_cor(method = "pearson",
           label.x = 8, label.y = 0.002) +

geom_point(aes(color = as.factor(Group)), size = 3, alpha = 0.6) +
  scale_color_manual(values = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))


#GC50%
ggscatter(Data_t, x = "X50_DNA", y = "X50_bx",
          add = "reg.line", conf.int = TRUE,
          add.params = list(fill = "lightgray"),
          ggtheme = theme_minimal()
          )+
  stat_cor(method = "pearson",
           label.x = 8, label.y = 0.002) +

geom_point(aes(color = as.factor(Group)), size = 3, alpha = 0.6) +
  scale_color_manual(values = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#GC53%
ggscatter(Data_t, x = "X53_DNA", y = "X53_bx",
          add = "reg.line", conf.int = TRUE,
          add.params = list(fill = "lightgray"),
          ggtheme = theme_minimal()
          )+
  stat_cor(method = "pearson",
           label.x = 8, label.y = 0.002) +

geom_point(aes(color = as.factor(Group)), size = 3, alpha = 0.6) +
  scale_color_manual(values = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#GC57%
ggscatter(Data_t, x = "X57_DNA", y = "X57_bx",
          add = "reg.line", conf.int = TRUE,
          add.params = list(fill = "lightgray"),
          ggtheme = theme_minimal()
          )+
  stat_cor(method = "pearson",
           label.x = 8, label.y = 0.002) +

geom_point(aes(color = as.factor(Group)), size = 3, alpha = 0.6) +
  scale_color_manual(values = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#GC60%
ggscatter(Data_t, x = "X60_DNA", y = "X60_bx",
          add = "reg.line", conf.int = TRUE,
          add.params = list(fill = "lightgray"),
          ggtheme = theme_minimal()
          )+
  stat_cor(method = "pearson",
           label.x = 8, label.y = 0.002) +

geom_point(aes(color = as.factor(Group)), size = 3, alpha = 0.6) +
  scale_color_manual(values = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#GC63%
ggscatter(Data_t, x = "X63_DNA", y = "X63_bx",
          add = "reg.line", conf.int = TRUE,
          add.params = list(fill = "lightgray"),
          ggtheme = theme_minimal()
          )+
  stat_cor(method = "pearson",
           label.x = 8, label.y = 0.002) +

geom_point(aes(color = as.factor(Group)), size = 3, alpha = 0.6) +
  scale_color_manual(values = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))




#Figures taken to eps-files, then axes are modified in Illustrator.
