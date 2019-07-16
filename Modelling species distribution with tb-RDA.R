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

setwd(choose.dir())
dir()

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
environment <- read.table(file.choose(),row.names=1,header=T,sep=" ") #environmental variables
View(environment)
dim(environment)
#environment<-decostand(environment,"standardize") #if you wish to standardize the scale of environmental variables
#View(environment)

##Spatial matrix
ll<-read.table(file.choose(),row.names=1,header=T,sep=",") #longitude and latitude
View(ll)
dim(ll)


#### CONDENSING SPECIES MATRIX WITH A PCoA (OPTIONAL) ##########
spp.dist <- vegdist(spp.u, method="bray") #Bray-Curtis dissimilarity index
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

# Alternatively, you may retain the significant axes (Broken-stick's approach):
ev <- pcoa.species$eig
(ev.r <- ev[1:30])
n <- length(ev)
bsm <- data.frame(j=seq(1:n), p=0)
bsm$p[1] <- 1/n
for (i in 2:n)
{
  bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))
}
bsm$p <- 100*bsm$p/n
bsm
dev.new(title="PCoA eigenvalues")
barplot(t(cbind(100*ev/sum(ev),bsm$p[n:1])), beside=TRUE, 
        main="% variance", col=c("bisque",2), las=2)
legend("topright", c("% eigenvalue", "Broken stick model"), 
       pch=15, col=c("bisque",2), bty="n")
source("evplot.R")
#evplot(ev.r) #In case of using ev.r, please use only the lower graphic.
evplot(ev)  #In case of using ev, you can use both graphics.

# Alternatively, you may retain the significant axes (permutation's approach):
library(PCPS)
pcoa.sig <- pcoa.sig(spp.u, method="bray", axis=6, by=100, iterations=1000)

# If, for instance, the first three axes are retained in the above analysis:
(scores.pcoa.3 <- scores.pcoa[,1:3])
head(scores.pcoa.3,20)
dim(scores.pcoa.3)

# Plot of the sites
dev.new(title="PCoA")
ordiplot(scores(pcoa.species, choices=c(1,2)), type="t", main="PCoA")
abline(h=0, lty=3)
abline(v=0, lty=3)
# Add weighted average projection of species
spe.wa <- wascores(pcoa.species$points[,1:2], spp.u)
text(spe.wa, rownames(spe.wa), cex=0.7, col="red")
write.table(spe.wa,"Species_Scores_PCoA-vazante.csv")
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
vif(environment)
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
env.rda <- rda(scores.pcoa.3,environment) #If you prefer to work with the whole response matrix, please change 'scores.pcoa.3' by 'spp.h'
env.rda
anova(env.rda, permutations = how(nperm=999))
### According to Blanchet et al. (2008): "If, and only if, the global 
### test is significant, one can proceed with forward selection"
(env.R2a <- RsquareAdj(env.rda)$adj.r.squared)
env.fwd <- forward.sel(scores.pcoa.3, environment)
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
(W_sel_fwd <- listw.select(scores.pcoa.3, candidates, MEM.autocor = "positive", method = "FWD",
                    p.adjust = TRUE, MEM.all = FALSE, nperm = 999)) 
#If you prefer to work with the whole response matrix, please change 'scores.pcoa.3' by 'spp.h' above.

save.image()

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

save.image()


#############################################
########## VARIATION PARTITIONING ###########
#############################################
vprda <- varipart(scores.pcoa.3, env.red, spatial.red, scale=TRUE) #classic variation partitioning
vprda

vprdaMSR <- msr(vprda, mem.all, nrepet = 999) #new variation partitioning (Clappe et al. 2018)
vprdaMSR

res <- rbind(vprda$R2.adj, vprdaMSR$R2.adj.msr)

rownames(res) <- c("Standard VP", "MSR VP")
res


### Alternatively, one can use the whole response matrix:
#vprda.all <- varipart(ssp.h, env.red, spatial.red, scale=TRUE) #classic variation partitioning
#vprda.all

#vprdaMSR <- msr(vprda, mem.all, nrepet = 999) #new variation partitioning 
#vprdaMSR

#res <- rbind(vprda.all$R2.adj, vprdaMSR$R2.adj.msr)
#rownames(res) <- c("Standard VP", "MSR VP")
#res

###########################################################################
####### Testing the environmental significance, after considering #########  
################# the effect of selected MEMs #############################
###########################################################################
env.rda<-rda(scores.pcoa.3,env.red,spatial.red)
head(summary(env.rda)) ### Please observe the explanation of each axis.
spenvcor(env.rda) # species-environment correlation
intersetcor(env.rda) #"interset" correlation
envfit(env.rda,env.red) #fits environmental vectors or factors onto an ordination
test.env<-anova(env.rda, permutations = how(nperm=999))
test.env
#plot(env.rda)


###########################################################################
######## Testing the spatial significance, after considering ##############  
################# the effect of selected env ##############################
###########################################################################
spatial.rda<-rda(scores.pcoa.3,spatial.red,env.red)
head(summary(spatial.rda)) ### Please observe the explanation of each axis.
spenvcor(spatial.rda) # species-space correlation
intersetcor(spatial.rda) #"interset" correlation
envfit(spatial.rda, spatial.red) #fits environmental vectors or factors onto an ordination
test.spatial<-anova(spatial.rda, permutations = how(nperm=999))
test.spatial
#plot(spatial.rda)


###########################################################################
### Evaluating the whole model ######
###########################################################################
all.predictors <- cbind(env.red,spatial.red)
all<-rda(scores.pcoa.3,all.predictors)
head(summary(all)) ### Please observe the explanation of each axis.
teste.all<-anova(all, permutations = how(nperm=999))
teste.all
spenvcor(all) # species-predictors correlation
intersetcor(all) #"interset" correlation
envfit(all, all.predictors) #fits environmental vectors or factors onto an ordination
# To test each axis individually:
rda.formula <- rda(scores.pcoa.3~., data=all.predictors)
anova(rda.formula, by="axis")
#plot(all)

save.image()

#END
