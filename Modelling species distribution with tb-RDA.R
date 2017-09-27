####################################################################################################
### Modelling patterns of species distribution with environmental and spatial data: a tentative code
####################################################################################################
### By: Pedro V. Eisenlohr
### After Borcard et al. (2011):
### Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, Springer, 2011
### Released on: 26 September, 2017

### Please refer to this code as "Eisenlohr (2014) with improvements available at http://github.com/pedroeisenlohr/variancepartition".
### Eisenlohr, P.V. (2014) Modelling species distribution patterns with environmental and spatial data: a tentative code. 
### doi: 10.1007/s40415-014-0064-3.


### Acknowledgments:
### Daniel Borcard
### Danilo Mesquita Neves
### Mario Jose Marques-Azevedo

# Define the working directory.
# Each user should adjust this!
setwd("C:/Users/pedro/Dropbox/Ferramentas LabEc/Matrizes")
getwd()

# If you do not have the packages installed, please use the following commands:
install.packages("ape") 
install.packages("packfor", repos="http://R-Forge.R-project.org")
install.packages("spacemakeR", repos="http://R-Forge.R-project.org")
install.packages("ade4")
install.packages("spdep") 
install.packages("vegan")
install.packages("AEM", repos="http://R-Forge.R-project.org")
install.packages("PCNM", repos="http://R-Forge.R-project.org")
# PCNM and/or AEM package have/has experienced problems with this link. If it happens, please 
# install it manually (.zip files) from another source.
install.packages("usdm")
install.packages("ClustOfVar")
install.packages("FactoMineR")

# Load the required packages
# (vegan must be loaded after ade4 to avoid some conflicts)
library(ape)
library(packfor)
library(spacemakeR)
library(ade4)
library(spdep)
library(vegan)
library(AEM)
library(PCNM)
library(usdm)
library(ClustOfVar)
library(FactoMineR)

# Import the data from CSV files
spp<-read.table(file.choose(),row.names=1,header=T,sep=",") ###community data
dim(spp)
env1<-read.table(file.choose(),row.names=1,header=T,sep=",") ###environmental data
dim(env1)
ll<-read.table(file.choose(),row.names=1,header=T,sep=",") ###coordinate data (long lat)
dim(ll)
# Examples of such matrices can be found assessing Supplementary Material of Eisenlohr (2014).

# Be careful: the following step is recommended only if you are working with
# multivariate data community.
# Transform the data
unicates<-apply(spp,2,sum) ###removing unicates
spp.u<-spp[,which(unicates>1)]
dim(spp.u)
spp.h <- decostand (spp.u, "hellinger") ###See Legendre & Gallagher (2001)
write.table(spp.h,"Spp_Hellinger.csv")

### If you have categorical environmental variables, you should run the
### next step only for quantitative ones. In this example, the columns 1 and 2 represent
### quantitative variables, and variables 3 to 5 represent categorical ones. However, if
### you are working with PCA axes as predictors, you SHOULD NOT follow the next step.
env<-decostand(env1[,1:2],"standardize") ###standardization of environmental data
### Run the example below adapting it to your data. Here, 'substrate', 'shrubs' and 
### 'topography' are categorical variables, and are placed at the columns 3, 4 and 5,
### respectively.
# Recode environmental variables 3 to 5 into dummy binary variables
substrate <- model.matrix(~env1[,3])[,-1]
shrubs <- model.matrix(~env1[,4])[,-1]
topography <- model.matrix(~env1[,5])[,-1]

env2 <- cbind(env[,1:2], substrate, shrubs, topography)
edit (env2)

### If you don't have categorical variables and wish to standardize your
### quantitative ones (which is not the case with PCA axes as predictors), 
### run the following command:
env1<-decostand(env1,"standardize")

### If you do not have categorical variables, please replace env2 by env1 below:

# Checking collinearities
env.rda <- rda(spp.h, env2)
vif.cca(env.rda) ### Checking for collinearity. VIF should be < 10.
### According to Blanchet et al. (2008): "It is not recommended to use a stepwise
### procedure in situations in which there are collinear variables (Chatterjee 
### and Price 1977, Freedman et al. 1992)".
### So, if you have found strong collinearities, I recommend you remove them
### before proceeding with the next step. An interesting way to address this
### issue is applying vifcor (usdm package) or the functions of 'ClustOfVar' and
### 'FactoMineR' packages.

# Forward selection of the environmental variables
anova(env.rda, step=1000)
### According to Blanchet et al. (2008): "If, and only if, the global 
### test is significant, one can proceed with forward selection"
(env.R2a <- RsquareAdj(env.rda)$adj.r.squared)
env.fwd <- forward.sel(spp.h, env2, adjR2thresh=env.R2a)
env.fwd #List of selected variables.
env.sign <- sort(env.fwd$order)
env.red <- env2[,c(env.sign)]

## SPATIAL ANALYSIS
# Is there a linear trend in the species distribution data?
### "As it is the case in most techniques of spatial analysis, the first step
### is to detrend the data. We recommend doing it whenever a significant linear
### trend is detected even though our method is able to model linear trends,
### because this preliminary step allows a separate modelling of the trend while
### retaining all the potential of the principal coordinates to model more 
### complex features" (Borcard & Legendre 2002).
spp.XY.rda<-rda(spp.h,ll)
anova<-anova.cca(spp.XY.rda,step=1000)
anova
# If ANOVA is significant, proceed with linearly detrended community data
spp.h.det <- resid(lm(as.matrix(spp.h) ~ ., data=ll))
write.table(spp.h.det,"Detrended_Data.csv")

## If ANOVA is significant, forward selection of coordinates
(spp.XY.R2a <- RsquareAdj(spp.XY.rda)$adj.r.squared)
(spp.XY.fwd <- forward.sel(spp.h, as.matrix(ll), 
	adjR2thresh=spp.XY.R2a))
XY.sign <- sort(spp.XY.fwd$order)
# Write the significant coordinates to a new object
XY.red <- ll[,c(XY.sign)]

### One can proceed further analysis with undetrended data, and partitioning out
### the variance with the selected linear trend variables. This is the case with
### the current code. If you decide to work with detrended data, please change
### 'spp.h" by 'spp.h.det' from now on.

################ MEM TYPE 1: Classic PCNM #############################
## Construct the PCNM variables
xy.d1 <- dist(ll)
spp.PCNM.auto <- PCNM(xy.d1)
summary(spp.PCNM.auto)
# Truncation distance
(dmin <- spp.PCNM.auto$thresh)
# Number of eigenvalues
(nb.ev <- length(spp.PCNM.auto$values))
# Expected value of I, no spatial correlation
spp.PCNM.auto$expected_Moran
# Moran's I of the PCNM variables
spp.PCNM.auto$Moran_I
# Eigenfunctions with positive spatial correlation
(select <- which(spp.PCNM.auto$Moran_I$Positive == TRUE))
# Number of PCNM with I > E(I)
length(select)
spp.PCNM.pos <- as.data.frame(spp.PCNM.auto$vectors)[,select]

## Run the global PCNM analysis on the *undetrended* data
(spp.PCNM.rda1 <- rda(spp.h, spp.PCNM.pos))
anova(spp.PCNM.rda1, step=1000)

## Since the analysis is significant, compute the adjusted R2 and run a forward
## selection of the PCNM variables. Below, one can perform the analysis with the
## *undetrended* data. 
(spp.R2a <- RsquareAdj(spp.PCNM.rda1)$adj.r.squared)
(spp.PCNM.fwd <- forward.sel(spp.h, as.matrix(spp.PCNM.pos), 
	adjR2thresh=spp.R2a))
(nb.sig.PCNM <- nrow(spp.PCNM.fwd)) #Number of signif. PCNM
# Identity of significant PCNMs in increasing order
(PCNM.sign <- sort(spp.PCNM.fwd[,2]))
# Write the significant PCNMs to a new object
PCNM.red <- spp.PCNM.pos[,c(PCNM.sign)]
write.table(PCNM.red,"PCNMsel.csv")
### Here it's important to register the R2a after PCNM selection. This value will
### be compared with the R2a obtained for alternative spatial models (MEMs 2-16;
### see below), which should also be recorded.

################ MEM TYPE 2 #############################
# Delaunay triangulation. No weighting matrix (binary weights).
ll.temp2 <- tri2nb(ll)
ll.transf2 = nb2listw(ll.temp2)
mem2<-scores.listw(ll.transf2) 
colnames(mem2$vectors)<-paste("MEMT2",1:ncol(mem2$vectors))
colnames(mem2$vectors)
mem.teste2<-test.scores(mem2,ll.transf2,1000)
mem.teste2
sel.mem2<-mem2
sel.mem2$vectors<-sel.mem2$vectors[,mem.teste2[,2]<0.05]
dim(sel.mem2$vectors)
#MEM with positive spatial autocorrelation
MEM.Moran.pos2<-which(mem.teste2[,1]>-1/(nrow(mem2$vectors)-1))
MEM.pos2<-mem2$vectors[,MEM.Moran.pos2]
#MEM with positive and significant spatial autocorrelation
MEM.Moran.pos.sig2<-MEM.Moran.pos2[which(mem.teste2[MEM.Moran.pos2,2]<=0.05)]
MEM.pos.sig2<-mem2$vectors[,MEM.Moran.pos.sig2]
dim(MEM.pos.sig2)
rda2<-rda(spp.h,MEM.pos.sig2)
r2.rda2<-RsquareAdj(rda2)$adj.r.squared
r2.rda2
anova.cca(rda2,step=1000) #Global test of RDA result. 
### According to Blanchet et al. (2008): "If, and only if, the global test is
### significant, one can proceed with forward selection"
rda2.sel<-forward.sel(spp.h,MEM.pos.sig2,adjR2thresh=r2.rda2)
rda2.sel #List of selected variables.
espaciais<-colnames(mem2$vectors)%in%rda2.sel$variables
sel.values<-sel.mem2$values[rda2.sel$order]
sel.vectors2<-sel.mem2$vectors[,rda2.sel$order]
write.table(sel.vectors2,"final.sel_MEM-T2.csv")

################ MEM TYPE 3 #############################
# Delaunay triangulation. "Minmax" matrix. 
### The ?minmax? style is based on Kelejian and Prucha (2010),
### and divides the weights by the minimum of the maximum row sums
### and maximum column sums of the input weights.
ll.temp3<-tri2nb(ll) 
ll.transf3 = nb2listw(ll.temp3, style = "minmax")
mem3<-scores.listw(ll.transf3) 
colnames(mem3$vectors)<-paste("MEMT3",1:ncol(mem3$vectors))
colnames(mem3$vectors)
mem.teste3<-test.scores(mem3,ll.transf3,1000)
mem.teste3
sel.mem3<-mem3
sel.mem3$vectors<-sel.mem3$vectors[,mem.teste3[,2]<0.05]
dim(sel.mem3$vectors)
#MEM with positive spatial autocorrelation
MEM.Moran.pos3<-which(mem.teste3[,1]>-1/(nrow(mem3$vectors)-1))
MEM.pos3<-mem3$vectors[,MEM.Moran.pos3]
dim(MEM.pos3)
#MEM with positive and significant spatial autocorrelation
MEM.Moran.pos.sig3<-MEM.Moran.pos3[which(mem.teste3[MEM.Moran.pos3,2]<=0.05)]
MEM.pos.sig3<-mem3$vectors[,MEM.Moran.pos.sig3]
dim(MEM.pos.sig3)
rda3<-rda(spp.h,MEM.pos.sig3)
r2.rda3<-RsquareAdj(rda3)$adj.r.squared
r2.rda3
anova.cca(rda3,step=1000) #Global test of RDA result. 
### According to Blanchet et al. (2008): "If, and only if, the global test is
### significant, one can proceed with forward selection"
rda3.sel<-forward.sel(spp.h,MEM.pos.sig3,adjR2thresh=r2.rda3)
rda3.sel #List of selected variables.
espaciais<-colnames(mem3$vectors)%in%rda3.sel$variables
sel.values<-sel.mem3$values[rda3.sel$order]
sel.vectors3<-sel.mem3$vectors[,rda3.sel$order]
write.table(sel.vectors3,"final.sel_MEM-T3.csv")

################ MEM TYPE 4 #############################
# Delaunay triangulation. "B" matrix.
ll.temp4<-tri2nb(ll) 
ll.transf4 = nb2listw(ll.temp4, style = "B")
mem4<-scores.listw(ll.transf4) 
colnames(mem4$vectors)<-paste("MEMT4",1:ncol(mem4$vectors))
colnames(mem4$vectors)
mem.teste4<-test.scores(mem4,ll.transf4,1000)
sel.mem4<-mem4
sel.mem4$vectors<-sel.mem4$vectors[,mem.teste4[,2]<0.05]
dim(sel.mem4$vectors)
#MEM with positive spatial autocorrelation
MEM.Moran.pos4<-which(mem.teste4[,1]>-1/(nrow(mem4$vectors)-1))
MEM.pos4<-mem4$vectors[,MEM.Moran.pos4]
dim(MEM.pos4)
#MEM with positive and significant spatial autocorrelation
MEM.Moran.pos.sig4<-MEM.Moran.pos4[which(mem.teste4[MEM.Moran.pos4,2]<=0.05)]
MEM.pos.sig4<-mem4$vectors[,MEM.Moran.pos.sig4]
dim(MEM.pos.sig4)
rda4<-rda(spp.h,MEM.pos.sig4)
r2.rda4<-RsquareAdj(rda4)$adj.r.squared
r2.rda4
anova.cca(rda4,step=1000) #Global test of RDA result.
### According to Blanchet et al. (2008): "If, and only if, the global test is
### significant, one can proceed with forward selection"
rda4.sel<-forward.sel(spp.h,sel.mem4$vectors,adjR2thresh=r2.rda4)
rda4.sel #List of selected variables.
espaciais<-colnames(mem4$vectors)%in%rda4.sel$variables
sel.values<-sel.mem4$values[rda4.sel$order]
sel.vectors4<-sel.mem4$vectors[,rda4.sel$order]
write.table(sel.vectors4,"final.sel_MEM-T4.csv")

################ MEM TYPE 5 #############################
# Delaunay triangulation. "C" matrix.
ll.temp5<-tri2nb(ll) 
ll.transf5 = nb2listw(ll.temp5, style = "C")
mem5<-scores.listw(ll.transf5) 
colnames(mem5$vectors)<-paste("MEMT5",1:ncol(mem5$vectors))
colnames(mem5$vectors)
mem.teste5<-test.scores(mem5,ll.transf5,1000)
sel.mem5<-mem5
sel.mem5$vectors<-sel.mem5$vectors[,mem.teste5[,2]<0.05]
dim(sel.mem5$vectors)
#MEM with positive spatial autocorrelation
MEM.Moran.pos5<-which(mem.teste5[,1]>-1/(nrow(mem5$vectors)-1))
MEM.pos5<-mem5$vectors[,MEM.Moran.pos5]
dim(MEM.pos5)
#MEM with positive and significant spatial autocorrelation
MEM.Moran.pos.sig5<-MEM.Moran.pos5[which(mem.teste5[MEM.Moran.pos5,2]<=0.05)]
MEM.pos.sig5<-mem5$vectors[,MEM.Moran.pos.sig5]
dim(MEM.pos.sig5)
rda5<-rda(spp.h,MEM.pos.sig5)
r2.rda5<-RsquareAdj(rda5)$adj.r.squared
r2.rda5
anova.cca(rda5,step=1000) #Global test of RDA result.
### According to Blanchet et al. (2008): "If, and only if, the global test is
### significant, one can proceed with forward selection"
rda5.sel<-forward.sel(spp.h,sel.mem5$vectors,adjR2thresh=r2.rda5)
rda5.sel #List of selected variables.
espaciais<-colnames(mem5$vectors)%in%rda5.sel$variables
sel.values<-sel.mem5$values[rda5.sel$order]
sel.vectors5<-sel.mem5$vectors[,rda5.sel$order]
write.table(sel.vectors5,"final.sel_MEM-T5.csv")

################ MEM TYPE 6 #############################
# Delaunay triangulation. "U" matrix.
ll.temp6<-tri2nb(ll) 
ll.transf6 = nb2listw(ll.temp6, style = "U")
mem6<-scores.listw(ll.transf6) 
colnames(mem6$vectors)<-paste("MEMT6",1:ncol(mem6$vectors))
colnames(mem6$vectors)
mem.teste6<-test.scores(mem6,ll.transf6,1000)
sel.mem6<-mem6
sel.mem6$vectors<-sel.mem6$vectors[,mem.teste6[,2]<0.05]
dim(sel.mem6$vectors)
#MEM with positive spatial autocorrelation
MEM.Moran.pos6<-which(mem.teste6[,1]>-1/(nrow(mem6$vectors)-1))
MEM.pos6<-mem6$vectors[,MEM.Moran.pos6]
dim(MEM.pos6)
#MEM with positive and significant spatial autocorrelation
MEM.Moran.pos.sig6<-MEM.Moran.pos6[which(mem.teste6[MEM.Moran.pos6,2]<=0.05)]
MEM.pos.sig6<-mem6$vectors[,MEM.Moran.pos.sig6]
dim(MEM.pos.sig6)
rda6<-rda(spp.h,MEM.pos.sig6)
r2.rda6<-RsquareAdj(rda6)$adj.r.squared
r2.rda6
anova.cca(rda6,step=1000) #Global test of RDA result.
### According to Blanchet et al. (2008): "If, and only if, the global test is
### significant, one can proceed with forward selection"
rda6.sel<-forward.sel(spp.h,sel.mem6$vectors,adjR2thresh=r2.rda6)
rda6.sel #List of selected variables.
espaciais<-colnames(mem6$vectors)%in%rda6.sel$variables
sel.values<-sel.mem6$values[rda6.sel$order]
sel.vectors6<-sel.mem6$vectors[,rda6.sel$order]
write.table(sel.vectors6,"final.sel_MEM-T6.csv")

################ MEM TYPE 7 #############################
# Delaunay triangulation. "S" matrix.
ll.temp7<-tri2nb(ll) 
ll.transf7 = nb2listw(ll.temp7, style = "S")
mem7<-scores.listw(ll.transf7) 
colnames(mem7$vectors)<-paste("MEMT7",1:ncol(mem7$vectors))
colnames(mem7$vectors)
mem.teste7<-test.scores(mem7,ll.transf7,1000)
sel.mem7<-mem7
sel.mem7$vectors<-sel.mem7$vectors[,mem.teste7[,2]<0.05]
dim(sel.mem7$vectors)
#MEM with positive spatial autocorrelation
MEM.Moran.pos7<-which(mem.teste7[,1]>-1/(nrow(mem7$vectors)-1))
MEM.pos7<-mem7$vectors[,MEM.Moran.pos7]
dim(MEM.pos7)
#MEM with positive and significant spatial autocorrelation
MEM.Moran.pos.sig7<-MEM.Moran.pos7[which(mem.teste7[MEM.Moran.pos7,2]<=0.05)]
MEM.pos.sig7<-mem7$vectors[,MEM.Moran.pos.sig7]
dim(MEM.pos.sig7)
rda7<-rda(spp.h,MEM.pos.sig7)
r2.rda7<-RsquareAdj(rda7)$adj.r.squared
r2.rda7
anova.cca(rda7,step=1000) #Global test of RDA result.
### According to Blanchet et al. (2008): "If, and only if, the global test is
### significant, one can proceed with forward selection"
rda7.sel<-forward.sel(spp.h,sel.mem7$vectors,adjR2thresh=r2.rda7)
rda7.sel #List of selected variables.
espaciais<-colnames(mem7$vectors)%in%rda7.sel$variables
sel.values<-sel.mem7$values[rda7.sel$order]
sel.vectors7<-sel.mem7$vectors[,rda7.sel$order]
write.table(sel.vectors7,"final.sel_MEM-T7.csv")

################ MEM TYPE 8 #############################
# Gabriel graph
ll.temp8 <- graph2nb(gabrielneigh(as.matrix(ll)), sym=TRUE)
ll.transf8 = nb2listw(ll.temp8)
mem8<-scores.listw(ll.transf8) 
colnames(mem8$vectors)<-paste("MEMT8",1:ncol(mem8$vectors))
colnames(mem8$vectors)
mem.teste8<-test.scores(mem8,ll.transf8,1000)
sel.mem8<-mem8
sel.mem8$vectors<-sel.mem8$vectors[,mem.teste8[,2]<0.05]
dim(sel.mem8$vectors)
#MEM with positive spatial autocorrelation
MEM.Moran.pos8<-which(mem.teste8[,1]>-1/(nrow(mem8$vectors)-1))
MEM.pos8<-mem8$vectors[,MEM.Moran.pos8]
dim(MEM.pos8)
#MEM with positive and significant spatial autocorrelation
MEM.Moran.pos.sig8<-MEM.Moran.pos8[which(mem.teste8[MEM.Moran.pos8,2]<=0.05)]
MEM.pos.sig8<-mem8$vectors[,MEM.Moran.pos.sig8]
dim(MEM.pos.sig8)
rda8<-rda(spp.h,MEM.pos.sig8)
r2.rda8<-RsquareAdj(rda8)$adj.r.squared
r2.rda8
anova.cca(rda8,step=1000) #Global test of RDA result.
### According to Blanchet et al. (2008): "If, and only if, the global test is
### significant, one can proceed with forward selection"
rda8.sel<-forward.sel(spp.h,sel.mem8$vectors,adjR2thresh=r2.rda8)
rda8.sel #List of selected variables
espaciais<-colnames(mem8$vectors)%in%rda8.sel$variables
sel.values<-sel.mem8$values[rda8.sel$order]
sel.vectors8<-sel.mem8$vectors[,rda8.sel$order]
write.table(sel.vectors8,"final.sel_MEM-T8.csv")

################ MEM TYPE 9 #############################
# Relative neighbourhood
ll.temp9 <- graph2nb(relativeneigh(as.matrix(ll)), sym=TRUE)
ll.transf9 = nb2listw(ll.temp9)
mem9<-scores.listw(ll.transf9) 
colnames(mem9$vectors)<-paste("MEMT9",1:ncol(mem9$vectors))
colnames(mem9$vectors)
mem.teste9<-test.scores(mem9,ll.transf9,1000)
sel.mem9<-mem9
sel.mem9$vectors<-sel.mem9$vectors[,mem.teste9[,2]<0.05]
dim(sel.mem9$vectors)
#MEM with positive spatial autocorrelation
MEM.Moran.pos9<-which(mem.teste9[,1]>-1/(nrow(mem9$vectors)-1))
MEM.pos9<-mem9$vectors[,MEM.Moran.pos9]
dim(MEM.pos9)
#MEM with positive and significant spatial autocorrelation
MEM.Moran.pos.sig9<-MEM.Moran.pos9[which(mem.teste9[MEM.Moran.pos9,2]<=0.05)]
MEM.pos.sig9<-mem9$vectors[,MEM.Moran.pos.sig9]
dim(MEM.pos.sig9)
rda9<-rda(spp.h,MEM.pos.sig9)
r2.rda9<-RsquareAdj(rda9)$adj.r.squared
r2.rda9
anova.cca(rda9,step=1000) #Global test of RDA result.
### According to Blanchet et al. (2008): "If, and only if, the global test is
### significant, one can proceed with forward selection"
rda9.sel<-forward.sel(spp.h,sel.mem9$vectors,adjR2thresh=r2.rda9)
rda9.sel #List of selected variables
espaciais<-colnames(mem9$vectors)%in%rda9.sel$variables
sel.values<-sel.mem9$values[rda9.sel$order]
sel.vectors9<-sel.mem9$vectors[,rda9.sel$order]
write.table(sel.vectors9,"final.sel_MEM-T9.csv")

################ MEM TYPE 10 #############################
# Minimum spanning tree
ll.temp10 <- mst.nb(dist(ll))
ll.transf10 = nb2listw(ll.temp10)
mem10<-scores.listw(ll.transf10) 
colnames(mem10$vectors)<-paste("MEMT10",1:ncol(mem10$vectors))
colnames(mem10$vectors)
mem.teste10<-test.scores(mem10,ll.transf10,1000)
sel.mem10<-mem10
sel.mem10$vectors<-sel.mem10$vectors[,mem.teste10[,2]<0.05]
dim(sel.mem10$vectors)
#MEM with positive spatial autocorrelation
MEM.Moran.pos10<-which(mem.teste10[,1]>-1/(nrow(mem10$vectors)-1))
MEM.pos10<-mem10$vectors[,MEM.Moran.pos10]
dim(MEM.pos10)
#MEM with positive and significant spatial autocorrelation
MEM.Moran.pos.sig10<-MEM.Moran.pos10[which(mem.teste10[MEM.Moran.pos10,2]<=0.05)]
MEM.pos.sig10<-mem10$vectors[,MEM.Moran.pos.sig10]
dim(MEM.pos.sig10)
rda10<-rda(spp.h,MEM.pos.sig10)
r2.rda10<-RsquareAdj(rda10)$adj.r.squared
r2.rda10
anova.cca(rda10,step=1000) #Global test of RDA result.
### According to Blanchet et al. (2008): "If, and only if, the global test is
### significant, one can proceed with forward selection"
rda10.sel<-forward.sel(spp.h,sel.mem10$vectors,adjR2thresh=r2.rda10)
rda10.sel #List of selected variables
espaciais<-colnames(mem10$vectors)%in%rda10.sel$variables
sel.values<-sel.mem10$values[rda10.sel$order]
sel.vectors10<-sel.mem10$vectors[,rda10.sel$order]
write.table(sel.vectors10,"final.sel_MEM-T10.csv")

################ MEM TYPE 11-16 #############################
# Connectivity matrix based on a distance (radius around points)
## Construct the PCNM variables
xy.d1 <- dist(ll)
PCNM.auto <- PCNM(xy.d1)
# Truncation distance
(dmin <- PCNM.auto$thresh)
# Using the same truncation distance dmin as in the classical PCNM.
thresh4 <- dnearneigh(as.matrix(ll), 0, dmin*4)
thresh4.d1 <- nbdists(thresh4, as.matrix(ll))
# Weights as function of inverse distance
inv.dist <- lapply(thresh4.d1, function(x) 1-x/max(dist(ll)))

################ MEM TYPE 11 #############################
### A - Creation of spatial weighting matrix W. Argument "B".
invdist1 = nb2listw(thresh4, glist=inv.dist, style="B")
mem11 <- scores.listw(invdist1, echo=TRUE)
colnames(mem11$vectors)<-paste("MEMT11",1:ncol(mem11$vectors))
colnames(mem11$vectors)
MEM.Moran11 <- test.scores(mem11,invdist1, 1000)
MEM.Moran11
# MEM with significant spatial correlation
which(MEM.Moran11[,2] <= 0.05)
length(which(MEM.Moran11[,2] <= 0.05))
# Store the MEM vectors in new objects
# All MEM
invdist.MEM.vec11 <-mem11$vectors
# MEM with positive spatial correlation
MEM.Moran.pos11 <- which(MEM.Moran11[,1] > -1/(nrow(mem11$vectors)-1))
invdist.MEM.pos11 <- invdist.MEM.vec11[,MEM.Moran.pos11]
# MEM with positive *and significant* spatial correlation
MEM.Moran.pos.sig11 <- MEM.Moran.pos11[which(MEM.Moran11[MEM.Moran.pos11,2] <= 0.05)]
invdist.MEM.pos.sig11 <- invdist.MEM.vec11[,MEM.Moran.pos.sig11]
dim(invdist.MEM.pos.sig11)
rda11<-rda(spp.h,invdist.MEM.pos.sig11)
(undet.PCNM.R2a <- RsquareAdj(rda11)$adj.r.squared)
anova.cca(rda11,step=1000) #Global test of RDA result. 
### According to Blanchet et al. (2008): "If, and only if, the global test is
### significant, one can proceed with forward selection"
rda11.sel<-forward.sel(spp.h, as.matrix(invdist.MEM.pos.sig11), 
	adjR2thresh=undet.PCNM.R2a)
rda11.sel #List of selected variables.
espaciais11<-colnames(mem11$vectors)%in%rda11.sel$variables
sel.values<-mem11$values[rda11.sel$order]
sel.vectors11<-mem11$vectors[,rda11.sel$order]
write.table(sel.vectors11,"final.sel_MEM-T11.csv")

################ MEM TYPE 12 #############################
### B - Creation of spatial weighting matrix W. Argument "C"
invdist2 = nb2listw(thresh4, glist=inv.dist, style="C")
mem12 <- scores.listw(invdist2, echo=TRUE)
colnames(mem12$vectors)<-paste("MEMT12",1:ncol(mem12$vectors))
colnames(mem12$vectors)
MEM.Moran12 <- test.scores(mem12,invdist2, 1000)
MEM.Moran12
# MEM with significant spatial correlation
which(MEM.Moran12[,2] <= 0.05)
length(which(MEM.Moran12[,2] <= 0.05))
# Store the MEM vectors in new objects
# All MEM
invdist.MEM.vec12 <-mem12$vectors
# MEM with positive spatial correlation
MEM.Moran.pos12 <- which(MEM.Moran12[,1] > -1/(nrow(mem12$vectors)-1))
invdist.MEM.pos12 <- invdist.MEM.vec12[,MEM.Moran.pos12]
# MEM with positive *and significant* spatial correlation
MEM.Moran.pos.sig12 <- MEM.Moran.pos12[which(MEM.Moran12[MEM.Moran.pos12,2] <= 0.05)]
invdist.MEM.pos.sig12 <- invdist.MEM.vec12[,MEM.Moran.pos.sig12]
dim(invdist.MEM.pos.sig12)
rda12<-rda(spp.h,invdist.MEM.pos.sig12)
(undet.PCNM.R2a <- RsquareAdj(rda12)$adj.r.squared)
anova.cca(rda12,step=1000) #Global test of RDA result. 
### According to Blanchet et al. (2008): "If, and only if, the global test is
### significant, one can proceed with forward selection"
rda12.sel<-forward.sel(spp.h, as.matrix(invdist.MEM.pos.sig12), 
	adjR2thresh=undet.PCNM.R2a)
rda12.sel #List of selected variables.
espaciais12<-colnames(mem12$vectors)%in%rda12.sel$variables
sel.values<-mem12$values[rda12.sel$order]
sel.vectors12<-mem12$vectors[,rda12.sel$order]
write.table(sel.vectors12,"final.sel_MEM-T12.csv")

################ MEM TYPE 13 #############################
### C - Creation of spatial weighting matrix W. Argument "minmax"
invdist3 = nb2listw(thresh4, glist=inv.dist, style="minmax")
mem13 <- scores.listw(invdist3, echo=TRUE)
colnames(mem13$vectors)<-paste("MEMT13",1:ncol(mem13$vectors))
colnames(mem13$vectors)
MEM.Moran13 <- test.scores(mem13,invdist3, 1000)
MEM.Moran13
# MEM with significant spatial correlation
which(MEM.Moran13[,2] <= 0.05)
length(which(MEM.Moran13[,2] <= 0.05))
# Store the MEM vectors in new objects
# All MEM
invdist.MEM.vec13 <-mem13$vectors
# MEM with positive spatial correlation
MEM.Moran.pos13 <- which(MEM.Moran13[,1] > -1/(nrow(mem13$vectors)-1))
invdist.MEM.pos13 <- invdist.MEM.vec13[,MEM.Moran.pos13]
# MEM with positive *and significant* spatial correlation
MEM.Moran.pos.sig13 <- MEM.Moran.pos13[which(MEM.Moran13[MEM.Moran.pos13,2] <= 0.05)]
invdist.MEM.pos.sig13 <- invdist.MEM.vec13[,MEM.Moran.pos.sig13]
dim(invdist.MEM.pos.sig13)
rda13<-rda(spp.h,invdist.MEM.pos.sig13)
(undet.PCNM.R2a <- RsquareAdj(rda13)$adj.r.squared)
anova.cca(rda13,step=1000) #Global test of RDA result. 
### According to Blanchet et al. (2008): "If, and only if, the global test is
### significant, one can proceed with forward selection"
rda13.sel<-forward.sel(spp.h, as.matrix(invdist.MEM.pos.sig13), 
	adjR2thresh=undet.PCNM.R2a)
rda13.sel #List of selected variables.
espaciais13<-colnames(mem13$vectors)%in%rda13.sel$variables
sel.values<-mem13$values[rda13.sel$order]
sel.vectors13<-mem13$vectors[,rda13.sel$order]
write.table(sel.vectors13,"final.sel_MEM-T13.csv")

################ MEM TYPE 14 #############################
### D - Creation of spatial weighting matrix W. Argument "W".
invdist4 = nb2listw(thresh4, glist=inv.dist, style="minmax")
mem14 <- scores.listw(invdist4, echo=TRUE)
colnames(mem14$vectors)<-paste("MEMT14",1:ncol(mem14$vectors))
colnames(mem14$vectors)
MEM.Moran14 <- test.scores(mem14,invdist4, 1000)
MEM.Moran14
# MEM with significant spatial correlation
which(MEM.Moran14[,2] <= 0.05)
length(which(MEM.Moran14[,2] <= 0.05))
# Store the MEM vectors in new objects
# All MEM
invdist.MEM.vec14 <-mem14$vectors
# MEM with positive spatial correlation
MEM.Moran.pos14 <- which(MEM.Moran14[,1] > -1/(nrow(mem14$vectors)-1))
invdist.MEM.pos14 <- invdist.MEM.vec14[,MEM.Moran.pos14]
# MEM with positive *and significant* spatial correlation
MEM.Moran.pos.sig14 <- MEM.Moran.pos14[which(MEM.Moran14[MEM.Moran.pos14,2] <= 0.05)]
invdist.MEM.pos.sig14 <- invdist.MEM.vec14[,MEM.Moran.pos.sig14]
dim(invdist.MEM.pos.sig14)
rda14<-rda(spp.h,invdist.MEM.pos.sig14)
(undet.PCNM.R2a <- RsquareAdj(rda14)$adj.r.squared)
anova.cca(rda14,step=1000) #Global test of RDA result. 
### According to Blanchet et al. (2008): "If, and only if, the global test is
### significant, one can proceed with forward selection"
rda14.sel<-forward.sel(spp.h, as.matrix(invdist.MEM.pos.sig14), 
	adjR2thresh=undet.PCNM.R2a)
rda14.sel #List of selected variables.
espaciais14<-colnames(mem14$vectors)%in%rda14.sel$variables
sel.values<-mem14$values[rda14.sel$order]
sel.vectors14<-mem14$vectors[,rda14.sel$order]
write.table(sel.vectors14,"final.sel_MEM-T14.csv")

################ MEM TYPE 15 #############################
### E - Creation of spatial weighting matrix W. Argument "S".
invdist5 = nb2listw(thresh4, glist=inv.dist, style="S")
mem15 <- scores.listw(invdist5, echo=TRUE)
colnames(mem15$vectors)<-paste("MEMT15",1:ncol(mem15$vectors))
colnames(mem15$vectors)
MEM.Moran15 <- test.scores(mem15,invdist5, 1000)
MEM.Moran15
# MEM with significant spatial correlation
which(MEM.Moran15[,2] <= 0.05)
length(which(MEM.Moran15[,2] <= 0.05))
# Store the MEM vectors in new objects
# All MEM
invdist.MEM.vec15 <-mem15$vectors
# MEM with positive spatial correlation
MEM.Moran.pos15 <- which(MEM.Moran15[,1] > -1/(nrow(mem15$vectors)-1))
invdist.MEM.pos15 <- invdist.MEM.vec15[,MEM.Moran.pos15]
# MEM with positive *and significant* spatial correlation
MEM.Moran.pos.sig15 <- MEM.Moran.pos15[which(MEM.Moran15[MEM.Moran.pos15,2] <= 0.05)]
invdist.MEM.pos.sig15 <- invdist.MEM.vec15[,MEM.Moran.pos.sig15]
dim(invdist.MEM.pos.sig15)
rda15<-rda(spp.h,invdist.MEM.pos.sig15)
(undet.PCNM.R2a <- RsquareAdj(rda15)$adj.r.squared)
anova.cca(rda15,step=1000) #Global test of RDA result. 
### According to Blanchet et al. (2008): "If, and only if, the global test is
### significant, one can proceed with forward selection"
rda15.sel<-forward.sel(spp.h, as.matrix(invdist.MEM.pos.sig15), 
	adjR2thresh=undet.PCNM.R2a)
rda15.sel #List of selected variables.
espaciais15<-colnames(mem15$vectors)%in%rda15.sel$variables
sel.values<-mem15$values[rda15.sel$order]
sel.vectors15<-mem15$vectors[,rda15.sel$order]
write.table(sel.vectors15,"final.sel_MEM-T15.csv")

################ MEM TYPE 16 #############################
### F - Creation of spatial weighting matrix W. Argument "U".
invdist6 = nb2listw(thresh4, glist=inv.dist, style="U")
mem16 <- scores.listw(invdist6, echo=TRUE)
colnames(mem16$vectors)<-paste("MEMT16",1:ncol(mem16$vectors))
colnames(mem16$vectors)
MEM.Moran16 <- test.scores(mem16,invdist6, 1000)
MEM.Moran16
# MEM with significant spatial correlation
which(MEM.Moran16[,2] <= 0.05)
length(which(MEM.Moran16[,2] <= 0.05))
# Store the MEM vectors in new objects
# All MEM
invdist.MEM.vec16 <-mem16$vectors
# MEM with positive spatial correlation
MEM.Moran.pos16 <- which(MEM.Moran16[,1] > -1/(nrow(mem16$vectors)-1))
invdist.MEM.pos16 <- invdist.MEM.vec16[,MEM.Moran.pos16]
# MEM with positive *and significant* spatial correlation
MEM.Moran.pos.sig16 <- MEM.Moran.pos16[which(MEM.Moran16[MEM.Moran.pos16,2] <= 0.05)]
invdist.MEM.pos.sig16 <- invdist.MEM.vec16[,MEM.Moran.pos.sig16]
dim(invdist.MEM.pos.sig16)
rda16<-rda(spp.h,invdist.MEM.pos.sig16)
(undet.PCNM.R2a <- RsquareAdj(rda16)$adj.r.squared)
anova.cca(rda16,step=1000) #Global test of RDA result. 
### According to Blanchet et al. (2008): "If, and only if, the global test is 
### significant, one can proceed with forward selection"
rda16.sel<-forward.sel(spp.h, as.matrix(invdist.MEM.pos.sig16), 
	adjR2thresh=undet.PCNM.R2a)
rda16.sel #List of selected variables.
espaciais16<-colnames(mem16$vectors)%in%rda16.sel$variables
sel.values<-mem16$values[rda16.sel$order]
sel.vectors16<-mem16$vectors[,rda16.sel$order]
write.table(sel.vectors16,"final.sel_MEM-T16.csv")


---------------------------------------------------------------------------------

# Preparing selected variables for variance partitioning:
### Please change * by the number of the best spatial model.
### In case of the classical PCNM was the best MEM type, you should skip this step.
sel.vectors<-sel.mem*$vectors[,rda*.sel$order]
write.table(sel.vectors,"best.mem.csv")

### I recommend that you check whether your data indeed fits a linear model.
### A possibility is to model the relationship between the selected environmental
### -spatial variables and the axes of a tb-PCA (transformation-based PCA,
### constructed by means of Hellinger-transformation) yielded only by species data. 
### You can also use ordiresids {vegan} (see details in http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/ordiresids.html).

# Variance partitioning:
### In case of detrended data as response matrix along the analysis above, one should remove 
### the 'XY.red' component, and change 'spp.h' by 'spp.h.det' from now on.
### Additionally, in case of classical PCNM (MEM Type 1) was the best MEM type, please change 
### sel.vectors by PCNM.red from now on.
### In case of using detrended response data (spp.h.det instead of spp.h), you should not use
### XY.red during the steps below.

all.varpart<-varpart(spp.h,env.red,sel.vectors,XY.red)
all.varpart
plot(all.varpart)

### Testing the environmental significance, after considering the effect of 
### selected MEMs and XY coordinates:
env.rda<-rda(spp.h,env.red,cbind(sel.vectors, XY.red))
summary(env.rda) ### Please observe the explanation of each axis.
teste.env<-anova(env.rda, step=1000)
teste.env
plot(env.rda)

### Testing the MEMs significance, after considering the effect of selected 
### environmental variables and XY coordinates:
esp.rda<-rda(spp.h,sel.vectors,cbind(XY.red,env.red)) 
summary(esp.rda) ### Please observe the explanation of each axis.
teste.esp<-anova(esp.rda, step=1000)
teste.esp
plot(esp.rda)

### Testing the XY signficance, after considering the effect of selected 
### environmental variables and MEMs:
esp.rda2<-rda(spp.h,XY.red,cbind(env.red,sel.vectors)) 
summary(esp.rda2) ### Please observe the explanation of each axis.
teste.esp<-anova(esp.rda2, step=1000)
teste.esp
plot(esp.rda2)

### Testing the significance of all predictors:
all<-rda(spp.h,cbind(env.red,sel.vectors,XY.red))
summary(all) ### Please observe the explanation of each axis.
teste.all<-anova(all, step=1000)
teste.all
plot(all)

### If the "pure" spatial fraction is significant, I recommend you proceed
### with testing if this fraction represents neutral processes. Mario Marques-
### Azevedo has implemented the method of Diniz-Filho et al. (2012) (see
### "Checking for missing environmental factor in spatial component" at
### https://github.com/MarioJose/scripts/tree/master/VarPartAndSpatialChecking).

###############################################################
### If you would like to clean all the steps above, please run:
### PLEASE NOTE: THE FOLLOWING COMMAND WILL ERASE ALL THE MEMORY OF THIS ROUTINE!
rm(list=ls(all=TRUE))
