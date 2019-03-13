### SELECION OF ENVIRONMENTAL AND SPATIAL DATA, AND VARIATION PARTITIONING ###
#By: Pedro V. Eisenlohr - UNEMAT, Alta Floresta (pedro.eisenlohr@unemat.br)
#Based on Bauman et al. (2018), Clappe et al. (2018), Blanchet et al. (2008), Chavent et al. (2012) and others

# Acknowledgments:
# Sylvie Clappe
# João Carlos Pires de Oliveira

library(ade4)
library(adespatial)
library(spdep)
library(vegan)
library(usdm)
library(ape)

##Species matrix
species<-read.table(file.choose(),row.names=1,header=T,sep=",") ###community species data (abundance or binary)
View(species)
dim(species)
unicates<-apply(species,2,sum) 
spp.u<-species[,which(unicates>1)] ###removing unicates
dim(spp.u)
spp.h <- decostand(spp.u, method = "hell")
View(spp.h)

##Environmental matrix
environment <- read.table(file.choose(),row.names=1,header=T,sep=",") #environmental variables
View(environment)
dim(environment)
#environment<-decostand(environment,"standardize") #caso queira padronizar a escala das variáveis ambientais
#View(environment)

##Spatial matrix
ll<-read.table(file.choose(),row.names=1,header=T,sep=",") #longitude and latitude
View(ll)
dim(ll)


#### CONDENSING SPECIES MATRIX WITH A PCoA (OPTIONAL) ##########
spp.dist <- vegdist(spp.u, method="bray")
pcoa.species <- cmdscale(spp.dist, k=(nrow(spp.u)-1), eig=TRUE)
scores.pcoa<-pcoa.species$points
head(scores.pcoa,30)
# Eigenvalues
(ev<-pcoa.species$eig)

# You may retain axes that capture 95% of all inertia:
(explic <- pcoa.species$eig/sum(pcoa.species$eig))
(cmu <- cumsum(explic))
(length.pcoa <- which.min(cmu <= 0.95))
(scores.pcoa <- scores.pcoa[,1:length.pcoa])
head(scores.pcoa,20)
dim(scores.pcoa)
##################################################################





#############################################
####### SELECTING ENVIRONMENTAL DATA ########
#############################################

### According to Blanchet et al. (2008): "It is not recommended to use a stepwise
### procedure in situations in which there are collinear variables (Chatterjee 
### and Price 1977, Freedman et al. 1992)".
### So, if you have found strong collinearities, I recommend you remove them
### before proceeding with the next step. An interesting way to address this
### issue is applying Hierarchical Clustering of Variables (Chavent et al. 2012).
# Checking collinearities
env.rda <- rda(scores.pcoa, environment) #If you prefer to work with the whole response matrix, please change 'scores.pcoa' by 'spp.h'
vif.cca(env.rda) ### Checking for collinearity. VIF should be < 10.
# If VIF > 10, consider using Hierarchical Clustering of Variables, as follows:


----------------------------------------------------------------------

######################################################################
#### HIERARCHICAL CLUSTERING OF VARIABLES (Chavent et al. 2012) ######
######################################################################
#library(ClustOfVar)
#tree <- hclustvar(environment)
#plot(tree)
#stab <- stability(tree,B=999) #To help in the selection of the number of partitions.

### Here, you need to adjust the routine to the number of selected clusters.
### For example, if this number is 4:
#P4<-cutreevar(tree,4,matsim=TRUE)
#cluster <- P4$cluster
#X <- environment
#princomp(X[,which(cluster==1)],cor=TRUE)$sdev^2
#princomp(X[,which(cluster==2)],cor=TRUE)$sdev^2
#princomp(X[,which(cluster==3)],cor=TRUE)$sdev^2
#princomp(X[,which(cluster==4)],cor=TRUE)$sdev^2
#P4$cluster
#clusterID <- P4$var
#clusterID
#write.table(P4$scores,"PCAScores.csv")

#env2 <- read.table("PCAScores.csv",row.names=1,header=T,sep=" ")
#dim(env2)
#View(env2)
#vif(env2) #Be sure that no variable presents VIF>10.

#In case of no colinearities, please run the command below:
#env3=env2
#In case of collinearities, please run the command below:
#(v1<-vifcor(env2,th=0.8)) 
#If necessary, please return to the above command and adjust the threshold.
#env3 <- exclude(env2, v1)

#names(env3)
#write.table(env3,"env_without_collinearities.csv")
#environment = env3
#View(environment)

#################################################################
#### END OF HIERARCHICAL CLUSTERING OF VARIABLES ################
#################################################################

-----------------------------------------------------------------

#################################################################
### FORWARD SELECTION OF PREDICTORS
#################################################################

# Forward selection of the environmental variables
env.rda<-rda(scores.pcoa,environment) #If you prefer to work with the whole response matrix, please change 'scores.pcoa' by 'spp.h'
anova(env.rda, permutations = how(nperm=999))
### According to Blanchet et al. (2008): "If, and only if, the global 
### test is significant, one can proceed with forward selection"
(env.R2a <- RsquareAdj(env.rda)$adj.r.squared)
env.fwd <- packfor::forward.sel(scores.pcoa, environment, adjR2thresh=env.R2a)
env.fwd #List of selected variables
env.sign <- sort(env.fwd$order)
env.red <- environment[,c(env.sign)]
head(env.red)

save.image()

#############################################
######### SELECTING SPATIAL DATA ############
#############################################

##### OPTIMIZING THE SELECTION OF SMW #####
### The function listw.candidates is used to build the spatial weighting matrices that
### we want to test and compare (with the listw.select function).
### I strongly recommend a careful reading on Bauman et al. (2018), mainly with respect
### the trade-off between accuracy and power analysis, since a p-value correction for 
### multiple tests (Sidak correction) is performed.

candidates <- listw.candidates(coord = ll, nb = c("del", "gab", "rel", "mst",
  			"pcnm", "dnear"), weights = c("binary", "flin", "fup", "fdown"))
names(candidates)                              
(nbw<-length(candidates)) ### Number of spatial weighting matrices generated


### Optimization the selection of a subset of SWM among the candidates generated above,
### using the corrected significance threshold calculated ("forward"):
(W_sel_fwd <- listw.select(scores.pcoa, candidates, MEM.autocor = "positive", method = "FWD",
                    p.adjust = TRUE, MEM.all = FALSE, nperm = 999)) 
#If you prefer to work with condensed response matrix, please change 'spp.h' by 'scores' above.

### Some characteristics of the best spatial model:
# Best SWM:
W_sel_fwd$best.id
(names.sel<-names(W_sel_fwd$best.id))
# Retained object for further analysis
SWM.selected<-W_sel_fwd$best$MEM.select
class(SWM.selected)
dim(SWM.selected)
#write.table(SWM.selected,"SWMselected.csv")
# Write the selected MEMs to a new object
spatial.red <- as.matrix(SWM.selected)
class(spatial.red)

# Creating an object with the name of the selected SWM:
x<-data.frame(names(candidates),W_sel_fwd$candidates$R2Adj.select)
names(x)<- c('MEMS', "select")
(ms<- which.max(x$select))
candidates1=array(candidates)
(candidates.sel<-candidates1[[ms]])

######################################################

### Obtaining all spatial eigenvectors from the SWM selected above:
mem.all<-mem(candidates.sel, MEM.autocor = "all")
mem.all
class(mem.all)
dim(mem.all)

### Obtaining positive spatial eigenvectors from the SWM selected above:
mem.positive<-mem(candidates.sel,MEM.autocor = "positive")
mem.positive
class(mem.positive)
dim(mem.positive)


#############################################
########## VARIATION PARTITIONING ###########
#############################################
vprda <- varipart(scores.pcoa, env.red, spatial.red, scale=TRUE) #classic variation partitioning
vprda

vprdaMSR <- msr(vprda, mem.all, nrepet = 999) #new variation partitioning (Clappe et al. 2018)
vprdaMSR

res <- rbind(vprda$R2.adj, vprdaMSR$R2.adj.msr)

rownames(res) <- c("Standard VP", "MSR VP")
res


### Alternatively, one can use the whole response matrix:
#vprda.all <- varipart(ssp.h, env.red, spatial.red, scale=TRUE) #classic variation partitioning
#vprda.all

#vprdaMSR <- msr(vprda.all, mem.all, nrepet = 999) #new variation partitioning 
#vprdaMSR

#res <- rbind(vprda.all$R2.adj, vprdaMSR$R2.adj.msr)
#rownames(res) <- c("Standard VP", "MSR VP")
#res

###########################################################################
####### Testing the environmental significance, after considering #########  
################# the effect of selected MEMs #############################
###########################################################################
env.rda<-rda(scores.pcoa,env.red,spatial.red)
head(summary(env.rda)) ### Please observe the explanation of each axis.
spenvcor(env.rda) # species-environment correlation
intersetcor(env.rda) #"interset" correlation
test.env<-anova(env.rda, permutations = how(nperm=999))
test.env
#plot(env.rda)

##############################################################################
#### Model validation: checking for spatial independence of RDA residuals ####
##############################################################################

##### Pending some adjustments recommended by David Bauman (personal contact) #####

#Generating objects required for the analysis:
X <- ll [,1]
Y <- ll [,2]
res <- residuals (env.rda)
cor.res<-dist(res)

#Mantel correlogram with Pearson correlation
#(cor.mantel<-mantel.correlog(cor.res, XY=ll, nperm=9999))
#summary(cor.mantel)
#plot(cor.mantel)

#Mantel correlogram with Spearman correlation
(mite.correlog2 <- mantel.correlog(cor.res, XY=ll, cutoff=FALSE, 
                                  r.type="spearman", nperm=99))
summary(mite.correlog2)
plot(mite.correlog2)


#Testing significance for each predictor
#Creating objects for each predictor in env.red:
### For example, if your predictors are oxy, slo, flo and pho:
oxy<-data.frame(env.red$oxy)
slo<-data.frame(env.red$slo)
flo<-data.frame(env.red$flo)
pho<-data.frame(env.red$pho)

#Testing significance of each predictor:
env.rda.oxy <-rda(spp.h,oxy,cbind(slo,flo,pho,spatial.red))
summary(env.rda.oxy)
(anova(env.rda.oxy,permutations = how(nperm=999)))

env.rda.slo <-rda(spp.h,slo,cbind(oxy,flo,pho,spatial.red))
summary(env.rda.slo)
(anova(env.rda.slo,permutations = how(nperm=999)))

env.rda.flo <-rda(spp.h,flo,cbind(oxy,slo,pho,spatial.red))
summary(env.rda.flo)
(anova(env.rda.flo,permutations = how(nperm=999)))

env.rda.pho <-rda(spp.h,pho,cbind(oxy,flo,slo,spatial.red))
summary(env.rda.pho)
(anova(env.rda.pho,permutations = how(nperm=999)))


###########################################################################
######## Testing the spatial significance, after considering ##############  
################# the effect of selected env ##############################
###########################################################################
spatial.rda<-rda(scores.pcoa,spatial.red,env.red)
head(summary(spatial.rda)) ### Please observe the explanation of each axis.
spenvcor(spatial.rda) # species-space correlation
intersetcor(spatial.rda) #"interset" correlation
test.spatial<-anova(spatial.rda, permutations = how(nperm=999))
test.spatial
#plot(spatial.rda)


###########################################################################
### Evaluating the whole model ######
###########################################################################
all.predictors <- cbind(env.red,spatial.red)
all<-rda(scores.pcoa,all.predictors)
head(summary(all)) ### Please observe the explanation of each axis.
teste.all<-anova(all, permutations = how(nperm=999))
teste.all
#plot(all)
# To test each axis individually:
rda.formula <- rda(scores.pcoa~., data=all.predictors)
anova(rda.formula, by="axis")

#Testing significance for each predictor.
### In this case, the MEMs selected above were MEM3,MEM6,MEM4,MEM1,MEM9,MEM10,MEM2,MEM7,MEM16,MEM11,MEM8
#Creating objects for each predictor in spatial.red:
spatial.red<-as.data.frame(spatial.red)
MEM3<-data.frame(spatial.red$MEM3)
MEM6<-data.frame(spatial.red$MEM6)
MEM4<-data.frame(spatial.red$MEM4)
MEM1<-data.frame(spatial.red$MEM1)
MEM9<-data.frame(spatial.red$MEM9)
MEM10<-data.frame(spatial.red$MEM10)
MEM2<-data.frame(spatial.red$MEM2)
MEM7<-data.frame(spatial.red$MEM7)
MEM16<-data.frame(spatial.red$MEM16)
MEM11<-data.frame(spatial.red$MEM11)
MEM8<-data.frame(spatial.red$MEM8)
all.spatial<-cbind(MEM3,MEM6,MEM4,MEM1,MEM9,MEM10,MEM2,MEM7,MEM16,MEM11,MEM8)
View(all.spatial)

#Testing significance of each predictor:
spatial.rda.MEM3 <-rda(scores.pcoa,MEM3,cbind(all.spatial[,-c(1)],env.red))
(anova(spatial.rda.MEM3,permutations = how(nperm=999)))

spatial.rda.MEM6 <-rda(scores.pcoa,MEM6,cbind(all.spatial[,-c(2)],env.red))
(anova(spatial.rda.MEM6,permutations = how(nperm=999)))

spatial.rda.MEM4 <-rda(scores.pcoa,MEM4,cbind(all.spatial[,-c(3)],env.red))
(anova(spatial.rda.MEM4,permutations = how(nperm=999)))

spatial.rda.MEM1 <-rda(scores.pcoa,MEM1,cbind(all.spatial[,-c(4)],env.red))
(anova(spatial.rda.MEM1,permutations = how(nperm=999)))

spatial.rda.MEM9 <-rda(scores.pcoa,MEM9,cbind(all.spatial[,-c(5)],env.red))
(anova(spatial.rda.MEM9,permutations = how(nperm=999)))

spatial.rda.MEM10 <-rda(scores.pcoa,MEM10,cbind(all.spatial[,-c(6)],env.red))
(anova(spatial.rda.MEM10,permutations = how(nperm=999)))

spatial.rda.MEM2 <-rda(scores.pcoa,MEM2,cbind(all.spatial[,-c(7)],env.red))
(anova(spatial.rda.MEM2,permutations = how(nperm=999)))

spatial.rda.MEM7 <-rda(scores.pcoa,MEM7,cbind(all.spatial[,-c(8)],env.red))
(anova(spatial.rda.MEM7,permutations = how(nperm=999)))

spatial.rda.MEM16 <-rda(scores.pcoa,MEM16,cbind(all.spatial[,-c(9)],env.red))
(anova(spatial.rda.MEM16,permutations = how(nperm=999)))

spatial.rda.MEM11 <-rda(scores.pcoa,MEM11,cbind(all.spatial[,-c(10)],env.red))
(anova(spatial.rda.MEM11,permutations = how(nperm=999)))

spatial.rda.MEM8 <-rda(scores.pcoa,MEM8,cbind(all.spatial[,-c(11)],env.red))
(anova(spatial.rda.MEM8,permutations = how(nperm=999)))

#END
