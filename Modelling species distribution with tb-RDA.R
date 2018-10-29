### SELECION OF ENVIRONMENTAL AND SPATIAL DATA, AND VARIATION PARTITIONING ###
#By: Pedro V. Eisenlohr - UNEMAT, Alta Floresta
#Based on Bauman et al. (2018), Clappe et al. (2018), Blanchet et al. (2008) and others

# Acknowledgments:
# Sylvie Clappe
# João Carlos Pires de Oliveira


library(ade4)
library(adespatial)
library(spdep)
library(vegan)
library(packfor)
library(usdm)

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
environment<-decostand(environment,"standardize") #caso queira padronizar a escala das variáveis ambientais
View(environment)

##Spatial matrix
ll<-read.table(file.choose(),row.names=1,header=T,sep=",") #longitude and latitude
View(ll)



#### CONDENSING SPECIES MATRIX ##########
(length.pca <- (nrow(spp.h)-1))
pca.species.p <- dudi.pca(spp.h, scannf = F, nf = length.pca)
summary(pca.species.p)
# Eigenvalues
(ev<-pca.species.p$eig)

# Apply Kaiser-Guttman criterion to select axes
#(length.pca2<- length(ev[ev > mean(ev)]))
#pca.species <- dudi.pca(spp.h, scannf = F, nf = length.pca2)
#scores <- pca.species$l1
#head(scores)
#dim(scores)

# Alternatively, you may retain axes that capture 95% of all inertia:
(explic <- pca.species.p$eig/sum(pca.species.p$eig))
(cmu <- cumsum(explic))
(length.pca3 <- which.min(cmu <= 0.95))
pca.species <- dudi.pca(spp.h, scannf = F, nf = length.pca3)
scores <- pca.species$l1
head(scores)
dim(scores)
#################################




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
env.rda <- rda(spp.h, environment)
vif.cca(env.rda) ### Checking for collinearity. VIF should be < 10.
# If VIF > 10, consider using Hierarchical Clustering of Variables, as follows:


#################################################################
#################################################################
#### HIERARCHICAL CLUSTERING OF VARIABLES #######################
#################################################################
#library(ClustOfVar)
#tree <- hclustvar(environment)
#plot(tree)
#stab <- stability(tree,B=100) #To help in the selection of the number of partitions.

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
#vif(env2)#Be sure that no variable presents VIF>10.

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
#################################################################
#### END OF HIERARCHICAL CLUSTERING OF VARIABLES ################
#################################################################



#################################################################
### FORWARD SELECTION OF PREDICTORS
#################################################################

### IMPORTANT: The forward selection procedure suggested above uses 
### "spp.h" as response variable matrix. This can lead to different
### results in comparison to "scores" as response variable matrix.
### One possible strategy would be comparing the results of both models
### and select one of them based on an objective criterion (R² adjusted, 
### AIC etc.).

# Forward selection of the environmental variables
env.rda<-rda(spp.h,environment)
anova(env.rda, permutations = how(nperm=999))
### According to Blanchet et al. (2008): "If, and only if, the global 
### test is significant, one can proceed with forward selection"
(env.R2a <- RsquareAdj(env.rda)$adj.r.squared)
env.fwd <- packfor::forward.sel(spp.h, environment, adjR2thresh=env.R2a)
env.fwd #List of selected variables
env.sign <- sort(env.fwd$order)
(env.red <- environment[,c(env.sign)])


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
(W_sel_fwd <- listw.select(spp.h, candidates, MEM.autocor = "positive", method = "FWD",
                    p.adjust = TRUE, MEM.all = TRUE, nperm = 999))
### Some characteristics of the best spatial model:
# Best SWM:
W_sel_fwd$best.id
(names.sel<-names(W_sel_fwd$best.id))
# Retained object for further analysis
SWM.selected<-W_sel_fwd$best$MEM.select
class(SWM.selected)
dim(SWM.selected)
write.table(SWM.selected,"SWMselected.csv")
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


#############################################
########## VARIATION PARTITIONING ###########
#############################################
scores.dudi <- dudi.pca(spp.h, scannf = F, nf=length.pca3) #nf refers to the number of PCA components retained. Default is 2.

vprda <- varipart(scores.dudi, env.red, spatial.red, nrepet=999)#classic variation partitioning
vprda

vprdaMSR <- msr(vprda, mem.all, nrepet = 999)#new variation partitioning 
vprdaMSR

res <- rbind(vprda$R2.adj, vprdaMSR$R2.adj.msr)
rownames(res) <- c("Standard VP", "MSR VP")
res


###########################################################################
### Testing the environmental significance, after considering #############  
################# the effect of selected MEMs #############################
###########################################################################
env.rda<-rda(spp.h,env.red,spatial.red)
summary(env.rda) ### Please observe the explanation of each axis.
spenvcor(env.rda) # species-environment correlation
intersetcor(env.rda) #"interset" correlation
test.env<-anova(env.rda, permutations = how(nperm=999))
test.env
plot(env.rda)


#### Model validation: checking for spatial independence of RDA residuals ####

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
                                  r.type="spearman", nperm=9999))
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
### Testing the spatial significance, after considering ###################  
################# the effect of selected env ##############################
###########################################################################
spatial.rda<-rda(spp.h,spatial.red,env.red)
summary(spatial.rda) ### Please observe the explanation of each axis.
spenvcor(spatial.rda) # species-environment correlation
intersetcor(spatial.rda) #"interset" correlation
test.spatial<-anova(spatial.rda, permutations = how(nperm=999))
test.spatial
#plot(spatial.rda)

#END
