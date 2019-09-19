## Examining digitizing error
# Digitizer: Vera Weisbecker

# August 2015

library(geomorph)
library(abind)

# Rep1 

setwd("E:/Students/Cruise/Paper/rep/R1")

##Cranial

WomCrData <- read.csv("Crania ORIGINAL.csv", header = T, row.names = 1) # Read data
WomCrData <- t(WomCrData) # Transpose rows and columns
WomCrData <- WomCrData[complete.cases(WomCrData),] # Remove any empty rows
Rep1_Cranial <- arrayspecs(WomCrData, k=3, p=ncol(WomCrData)/3) # Make 3D Array

##Mandible
WomMaData <- read.csv("Mandibles ORIGINAL.csv", header = T, row.names = 1) 
WomMaData <- t(WomMaData)
WomMaData <- WomMaData[complete.cases(WomMaData),] 
Rep1_Mandible <- arrayspecs(WomMaData, k=3, p=ncol(WomMaData)/3) 

# Rep2 

setwd("E:/Students/Cruise/Paper/rep/R2")

##Cranial

WomCrData <- read.csv("Crania RELANDMARKED.csv", header = T, row.names = 1) # Read data
WomCrData <- t(WomCrData) # Transpose rows and columns
WomCrData <- WomCrData[complete.cases(WomCrData),] # Remove any empty rows
Rep2_Cranial <- arrayspecs(WomCrData, k=3, p=ncol(WomCrData)/3) # Make 3D Array

##Mandible
WomMaData <- read.csv("Mandibles RELANDMARKED.csv", header = T, row.names = 1) 
WomMaData <- t(WomMaData)
WomMaData <- WomMaData[complete.cases(WomMaData),] 
Rep2_Mandible <- arrayspecs(WomMaData, k=3, p=ncol(WomMaData)/3) 

#Check match of names
nCran <- length(dimnames(Rep2_Cranial)[[3]])
nMand <- length(dimnames(Rep2_Mandible)[[3]])
  
dimnames(Rep2_Cranial)[[3]] == dimnames(Rep1_Cranial)[[3]]
dimnames(Rep2_Mandible)[[3]] == dimnames(Rep1_Mandible)[[3]]

#Concatenate arrays using abind

CranRep <- abind(Rep1_Cranial,Rep2_Cranial)  
MandRep <- abind(Rep1_Mandible,Rep2_Mandible)

#Check all is well
plot3d(CranRep[,,3],aspect = FALSE)

# Procrustes of the two sets - there's a warning that row names are not used but that's OK because they have to be duplicates

CranRep_GPA <- gpagen(CranRep) 

MandRep_GPA <- gpagen(MandRep)

individuals <- dimnames(CranRep_GPA$coords)[[3]] # name for each pair shouldbe same, so not _1, _2


# create a vector designating specimens to a replicate, 12 each for both datasets
replicate <- c(rep(1, nCran), rep(2, nCran)) 

# to check your replicates, examine the replicates in a PCA
plotTangentSpace(CranRep_GPA$coords, groups = replicate) # PCA coloured by replicate; looking for pairs in PCA; black is rep2
plotTangentSpace(MandRep_GPA$coords, groups = replicate, label = individuals)

# Look at the outliers
out <- plotOutliers(gpa$coords) # Which specimens are way off?
plotRefToTarget(mshape(gpa$coords), gpa$coords[,,out[1]], method="vector", label = T)

## TEST Examining replicate error

rep.er.Cran <- procD.lm(CranRep_GPA$coords ~ factor(individuals))
((rep.er.Cran$aov.table$MS[1] - rep.er.Cran$aov.table$MS[2])/2) / (rep.er.Cran$aov.table$MS[2] + ((rep.er.Cran$aov.table$MS[1] - rep.er.Cran$aov.table$MS[2])/2))
# Number here should be above 0.9 for good repeatability, although in practice 

rep.er.Mand <- procD.lm(MandRep_GPA$coords ~ factor(individuals))
((rep.er.Mand$aov.table$MS[1] - rep.er.Mand$aov.table$MS[2])/2) / (rep.er.Mand$aov.table$MS[2] + ((rep.er.Mand$aov.table$MS[1] - rep.er.Mand$aov.table$MS[2])/2))
