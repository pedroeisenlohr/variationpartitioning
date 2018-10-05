### SELECION OF ENVIRONMENTAL AND SPATIAL DATA, AND VARIATION PARTITIONING ###
#By: Pedro V. Eisenlohr - UNEMAT, Alta Floresta (pedro.eisenlohr@unemat.br)
#Based on Bauman et al. (2018), Clappe et al. (2018), Blanchet et al. (2008) and others
### Released on: October 04, 2018

### Please refer to this code as "Eisenlohr (2014) with improvements available at http://github.com/pedroeisenlohr/variancepartition".
### Eisenlohr, P.V. (2014) doi: 10.1007/s40415-014-0064-3.

### Acknowledgments: 
### Sylvie Clappe
### João Carlos Pires de Oliveira


library(ade4)
library(adespatial)
library(spdep)
library(vegan)
library(packfor)

##Species matrix
species<-read.table(file.choose(),row.names=1,header=T,sep=",") ###community data
View(species)
unicates<-apply(spp,2,sum) 
spp.u<-spp[,which(unicates>1)] ###removing unicates
dim(spp.u)
spp.h <- decostand(spp.u, method = "hell")

##Environmental matrix
environment <- read.table(file.choose(),row.names=1,header=T,sep=",")
View(environment)
#environment<-decostand(environment,"standardize")
View(environment)

##Spatial matrix
ll<-read.table(file.choose(),row.names=1,header=T,sep=",")
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

# Alternatively, you may retain axes that capture about 95% of all inertia:
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
### issue is applying vifcor (usdm package).
# Checking collinearities
env.rda <- rda(scores, environment)
vif.cca(env.rda) ### Checking for collinearity. VIF should be < 10.
# If VIF > 10, you should remove one or more variables:
#v1<-vifcor(environment,th=0.8)
#v1
# Verify VIF. VIF should be < 10. If necessary, reduce the th.
#environment<-exclude(environment,v1)
#names(environment)

# Forward selection of the environmental variables
env.rda<-rda(scores,environment)
anova(env.rda, permutations = how(nperm=999))
### According to Blanchet et al. (2008): "If, and only if, the global 
### test is significant, one can proceed with forward selection"
(env.R2a <- RsquareAdj(env.rda)$adj.r.squared)
env.fwd <- forward.sel(scores, environment, adjR2thresh=env.R2a)
env.fwd #List of selected variables
length(env.fwd) #Number of selected variables
env.sign <- sort(env.fwd$order)
(env.red <- environment[,c(env.sign)])


#############################################
######### SELECTING SPATIAL DATA ############
#############################################

##### OPTIMIZING THE SELECTION OF SMW #####
candidates <- listw.candidates(coord = ll, nb = c("del", "gab", "rel", "mst",
  			"pcnm", "dnear"), weights = c("binary", "flin", "fup", "fdown"))
names(candidates)                              
(nbw<-length(candidates)) ### Number of spatial weighting matrices generated

### Optimization of the selection of the SWM among the candidates generated above,
### using the corrected significance threshold calculated above for the global tests (global):
(W_sel_all <- listw.select(scores, candidates, MEM.autocor = "all", method = "global",
                    p.adjust = TRUE, MEM.all = TRUE, nperm = 999))
### Some characteristics of the best spatial model:
# Best SWM:
(names.sel<-names(W_sel_all$best.id))
# Retained object for further analysis
W_sel_all$best$MEM.all
(length(W_sel_all$best$MEM.all))
MEM.all<-W_sel_all$best$MEM.all

### Forward selection of the best SWM selected above:
## Run the global RDA analysis
(spp.spatial.rda <- rda(scores, W_sel_all$best$MEM.select))
anova(spp.spatial.rda, permutations = how(nperm=999))
(spp.R2a <- RsquareAdj(spp.spatial.rda)$adj.r.squared)
(spp.spatial.fwd <- forward.sel(scores, as.matrix(W_sel_all$best$MEM.select), 
	adjR2thresh=spp.R2a))
(nb.sig.spatial <- nrow(spp.spatial.fwd)) #Number of signif. MEM
# Identity of significant MEMs in increasing order
(spatial.sign <- sort(spp.spatial.fwd[,2]))
# Write the significant MEMs to a new object
spatial.red <- W_sel_all$best$MEM.select[,c(spatial.sign)]
spatial.red <- as.matrix(spatial.red)
class(spatial.red)



#############################################
########## VARIATION PARTITIONING ###########
#############################################
scores.dudi <- dudi.pca(spp.h, scannf = F, nf=length.pca3)

vprda <- varipart(scores.dudi, env.red, spatial.red, nrepet=999)#classic variation partitioning
vprda

vprdaMSR <- msr(vprda, MEM.all, nrepet = 999)#new variation partitioning 
vprdaMSR

res <- rbind(vprda$R2.adj, vprdaMSR$R2.adj.msr)
rownames(res) <- c("Standard VP", "MSR VP")
res


###########################################################################
### Testing the environmental significance, after considering #############  
#################### the effect of selected MEMs ##########################
###########################################################################
env.rda<-rda(scores,env.red,spatial.red)
summary(env.rda) ### Please observe the explanation of each axis.
spenvcor(env.rda) # species-environment correlation
intersetcor(env.rda) #"interset" correlation
test.env<-anova(env.rda, permutations = how(nperm=999))
test.env
#plot(env.rda)

###########################################################################
### Testing the spatial significance, after considering ###################  
############### the effect of selected env ################################
###########################################################################
spatial.rda<-rda(scores,spatial.red,env.red)
summary(spatial.rda) ### Please observe the explanation of each axis.
spenvcor(spatial.rda) # species-environment correlation
intersetcor(spatial.rda) #"interset" correlation
test.spatial<-anova(spatial.rda, permutations = how(nperm=999))
test.spatial
#plot(spatial.rda)
