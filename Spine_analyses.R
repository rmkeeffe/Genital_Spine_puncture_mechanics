#article: Evolutionary Morphology of Genital Spines Informed by Puncture Mechanics
#published in: Proceedings of the Royal Society B
#doi: 10.1098/rspb.2025.1698

# 1 - Load Packages & Data ------------------------------------------------------------------

rm(list=ls(all=T)) 
library(geomorph) #v. 4.0.8
library(ape) #v. 5.8
library(phytools) #v. 2.3-0
library(geiger) #v. 2.0.11
library(calibrate) #v. 1.7.7
library(RevGadgets) #v. 1.2.1
library(dplyr) #v. 1.1.4
library(Morpho) #v. 2.12
library(nlme) #v. 3.1-164
library(plotrix) #v. 3.8-4
library(ggplot2) #v. 3.5.1

#Set working directory

#Read data files
SpineShape <- read.morphologika("spine_morphologika_unscaled.txt") #Load in TPS file
SpineGPA <- gpagen(SpineShape, ProcD = FALSE) #Run Generalized Procrustes Analysis using bending energy criterion
SpineData <- read.csv("Spine_Data.csv") #Load in csv file with SpecimenID, SVL (mm)

# 2 - Principal Component Analysis ------------------------------------------------------------------

#Run spine PCA
SpinePCA <- gm.prcomp(SpineGPA$coords, phy = NULL)
SpinePCASum <- summary(SpinePCA) #Prints variance and eigenvalues

#Set tested character to factor
SpineData$Tested_for_Puncture_Performance <- factor(SpineData$Tested_for_Puncture_Performance)

#Set family character to factor
SpineData$Family <- factor(SpineData$Family)

#Set species character to factor
SpineData$Species <- factor(SpineData$Species)

#write plot
SpinePCAPlot <- plot(SpinePCA,
                     xlab = paste("Principal Component 1 ", "(", sep = "", 
                                  paste(round(SpinePCASum$PC.summary$Comp1[2]*100, digits = 2), "%", ")", sep = "")),
                     ylab = paste("Principal Component 2 ", "(", sep = "", 
                                  paste(round(SpinePCASum$PC.summary$Comp2[2]*100, digits = 2), "%", ")", sep = "")),
                     pch = c(21,22,23)[SpineData$Tested_for_Puncture_Performance],
                     col = "black",
                     bg = c("blue","green", "orange","red","purple")[SpineData$Family],
                     cex = 2,
                     cex.lab = 1.2)
text(SpinePCAPlot$PC.points[ , 2] ~ SpinePCAPlot$PC.points[ , 1], 
     labels = SpineData$Specimen_ID, cex= 0.5, col = c("black"))

# 3 - Phylomorphospace ------------------------------------------------------------------

#Read Zaher et al. 2019 tree
snake.tree <- read.nexus("pone.0216148.s012.tre")

#Read species list and prune tree
species_list <- read.csv("phylomorph_pruning.csv", header = F)
tips <- species_list$V1
namecheck <- name.check(snake.tree, data.names=tips)
treenotdata <- as.vector(namecheck$tree_not_data)
datanottree <- as.vector(namecheck$data_not_tree)
droppedtip <- drop.tip(snake.tree,namecheck$tree_not_data)

#Plotting trimmed tree to check for errors
plot(droppedtip) 

#Preparing tree for phylomorphospace
tt <- multi2di(droppedtip)
tt$edge.length[tt$edge.length==0]<-1e-8

#Reading trimmed datasets (largest spine only) for phylomorphospace
SpineShapeTrim <- read.morphologika("trimmed_spine_morphologika_unscaled2.txt") #Load in TPS file
SpineGPATrim <- gpagen(SpineShapeTrim, ProcD = FALSE) #Run Generalized Procrustes Analysis using bending energy criterion
SpineDataTrim <- read.csv("trimmed_Spine_Data.csv") #Load in csv file with SpecimenID, SVL (mm)

#Run PCA on trimmed dataset
SpinePCATrim <- gm.prcomp(SpineGPATrim$coords, phy = NULL)
SpinePCATrimSum <- summary(SpinePCATrim) #Prints variance and eigenvalues
SpineDataTrim$Tested_for_Puncture_Performance <- factor(SpineDataTrim$Tested_for_Puncture_Performance)
SpineDataTrim$Family <- factor(SpineDataTrim$Family)

SpinePCAPlotTrim <- plot(SpinePCATrim)
SpinePCAPlotTrim$PC.points

#Read aligned dataset with PC coords
aligned_data <- read.csv("Aligned_data.csv", header = T, row.names = 1)

#Plot phylomorphospace
phylomorphospace(tt, aligned_data)

#Calculate the degree of phylogenetic signal from Procrustes shape variables
#This might take a few minutes to run
physig <- physignal(SpineGPATrim$coords, phy = tt, iter = 10000, print.progress = T)
#Observed Phylogenetic Signal (K): 0.49056
#P-value: 0.12449

# 4 - Statistics ------------------------------------------------------------------

fittest <- lm(SpineData$SVL_mm ~ SpineData$Volume)
summary(fittest)


#Make a geomorph dataframe for using functions in geomorph package
SpineGDF <- geomorph.data.frame(shape = SpineGPA$coords,
                                SVL = SpineData$SVL_mm,
                                VOL = SpineData$Volume,
                                FAM = SpineData$Family,
                                STRCURV = SpineData$Structural_Curvature_degrees,
                                SPE = SpineData$Species,
                                TIP = SpineData$Tip_angle_degrees,
                                ROC = SpineData$Radius_of_Curvature_at_Tip)

#Remove NAs from dataset
ndf <- na.omit(SpineGDF)

#Run a Procrustes ANOVA of shape and SVL, remove NAs
SpineFit1 <- procD.lm(shape ~ SVL,
                   data = ndf,
                   iter = 1000)
summary(SpineFit1)
###p = 0.6304

#Run a Procrustes ANOVA of shape and Family
SpineFit2 <- procD.lm(shape ~ FAM,
                      data = SpineGDF,
                      iter = 1000)
summary(SpineFit2)
###p = 0.6713

#Run a Procrustes ANOVA of shape and structural curvature
SpineFit3 <- procD.lm(shape ~ STRCURV,
                      data = SpineGDF,
                      iter = 1000)
summary(SpineFit3)
###p = 0.1518

#Run a Procrustes ANOVA of shape and species
SpineFit4 <- procD.lm(shape ~ SPE,
                      data = SpineGDF,
                      iter = 1000)
summary(SpineFit4)
###p = 0.985

#Run a Procrustes ANOVA of shape and tip angle
SpineFit5 <- procD.lm(shape ~ TIP,
                      data = SpineGDF,
                      iter = 1000)
summary(SpineFit5)
###p = 0.7253

#Run a Procrustes ANOVA of shape and radius of curvature
SpineFit6 <- procD.lm(shape ~ ROC,
                      data = SpineGDF,
                      iter = 1000)
summary(SpineFit6)
###p = 0.5215

#Subset Performance Data
PERF_DATA <- subset(SpineData, Tested_for_Puncture_Performance %in% c("Y"))

#remove NERH_Spine5 which broke during testing
PERF_DATA <- PERF_DATA[which(PERF_DATA$Specimen_ID != "NR30M_Spine5"), ]

#Run a linear model between average puncture force and ROC
Perf_Model1 <- lm(Average_Puncture_Force ~ Radius_of_Curvature_at_Tip, data=PERF_DATA)
summary(Perf_Model1)
###p = 0.9648

#Run a linear model between average puncture force and TA
Perf_Model2 <- lm(Average_Puncture_Force ~ Tip_angle_degrees, data=PERF_DATA)
summary(Perf_Model2)
###p = 0.05204

#Run a linear model between average puncture force and PC1 scores
Perf_Model3 <- lm(Average_Puncture_Force ~ PC1_value,
                  data=PERF_DATA)
summary(Perf_Model3)
###p = 0.5152

#Run a linear model between average puncture force and volume
Perf_Model4 <- lm(Average_Puncture_Force ~ Volume, data=PERF_DATA)
summary(Perf_Model4)
###p = 0.3999

#Run a linear model between average puncture force and structural curvature
Perf_Model5 <- lm(Average_Puncture_Force ~ Structural_Curvature_degrees, data=PERF_DATA)
summary(Perf_Model5)
###p = 0.5954

#Run a linear model between puncture angle range and structural curvature
AnglePerf_Model1 <- lm(Puncture_range_degrees ~ Structural_Curvature_degrees,
                  data=PERF_DATA)
summary(AnglePerf_Model1)
###p = 0.04498 *

#Run a linear model between puncture angle range and PC1 scores
AnglePerf_Model2 <- lm(Puncture_range_degrees ~ PC1_value,
                  data=PERF_DATA)
summary(AnglePerf_Model2)
###p = 0.3102

#Run a linear model between puncture angle range and ROC
AnglePerf_Model3 <- lm(Puncture_range_degrees ~ Radius_of_Curvature_at_Tip,
                       data=PERF_DATA)
summary(AnglePerf_Model3)
###p = 0.004972 *

#Run a linear model between puncture angle range and TA
AnglePerf_Model4 <- lm(Puncture_range_degrees ~ Tip_angle_degrees,
                       data=PERF_DATA)
summary(AnglePerf_Model4)
###p = 0.02451 *

#Run a linear model between puncture angle range and volume
AnglePerf_Model5 <- lm(Puncture_range_degrees ~ Volume,
                        data=PERF_DATA)
summary(AnglePerf_Model5)
###p = 0.3816

#Run a linear model between lowest force to puncture and TA
ForcePerf_Model6 <- lm(Lowest_Force_to_Puncture ~ Tip_angle_degrees,
                        data=PERF_DATA)
summary(ForcePerf_Model6)
###p = 0.07609 *

#Run a linear model between lowest force to puncture and ROC
ForcePerf_Model13 <- lm(Lowest_Force_to_Puncture ~ Radius_of_Curvature_at_Tip,
                        data=PERF_DATA)
summary(ForcePerf_Model13)
###p = 0.8406

#Run a linear model between lowest force to puncture and volume
ForcePerf_Model14 <- lm(Lowest_Force_to_Puncture ~ Volume,
                        data=PERF_DATA)
summary(ForcePerf_Model14)
###p = 0.3195

#Run a linear model between lowest force to puncture and PC1 score
ForcePerf_Model14 <- lm(Lowest_Force_to_Puncture ~ PC1_value,
                        data=PERF_DATA)
summary(ForcePerf_Model14)
###p = 0.4311

#Run a linear model between lowest force to puncture and structural curvature
ForcePerf_Model15 <- lm(Lowest_Force_to_Puncture ~ Structural_Curvature_degrees,
                        data=PERF_DATA)
summary(ForcePerf_Model15)
###p = 0.8175

# 5 - Regression Plots ------------------------------------------------------------------

#Write plot for puncture angle range vs. structural curvature
Plot_1 <- plot(PERF_DATA$Puncture_range_degrees ~ PERF_DATA$Structural_Curvature_degrees,
                          data = PERF_DATA,
                          pch = 21,
                          col = "black",
                          bg = "blue",
                          cex = 2,
                          cex.lab = 1.2,
                          xlab = "Structural Curvature (degrees)",
                          ylab = "Puncture Angle Range (degrees)")
abline(a = 30.6931,
       b =  0.7832,
       lw = 3)

#Write plot for puncture angle range vs. tip angle
Plot_2 <- plot(PERF_DATA$Puncture_range_degrees ~ PERF_DATA$Tip_angle_degrees,
                          data = PERF_DATA,
                          pch = 21,
                          col = "black",
                          bg = "blue",
                          cex = 2,
                          cex.lab = 1.2,
                          xlab = "Tip Angle (degrees)",
                          ylab = "Puncture Angle Range (degrees)")
abline(a = 88.0525,
       b =  -0.7702,
       lw = 3)

#Write plot for ROC vs. puncture performance score
Plot_3 <- plot(PERF_DATA$Puncture_range_degrees ~ PERF_DATA$Radius_of_Curvature_at_Tip,
                          data = PERF_DATA,
                          pch = 21,
                          col = "black",
                          bg = "blue",
                          cex = 2,
                          cex.lab = 1.2,
                          xlab = "Radius of Curvature at Tip",
                          ylab = "Puncture Angle Range (degrees)")
abline(a = 77.264,
       b = -364.765,
       lw = 3)

Plot_4 <- plot(PERF_DATA$Average_Puncture_Force ~ PERF_DATA$Tip_angle_degrees,
               data = PERF_DATA,
               pch = 21,
               col = "black",
               bg = "blue",
               cex = 2,
               cex.lab = 1.2,
               xlab = "Tip Angle (degrees)",
               ylab = "Average Puncture Force (N)")
abline(a = -7.2639,
       b = 0.7404,
       lw = 3)
