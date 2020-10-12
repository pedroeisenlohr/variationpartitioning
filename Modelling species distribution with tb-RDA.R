### SELECION OF ENVIRONMENTAL AND SPATIAL DATA, AND VARIATION PARTITIONING ###
#By: Pedro V. Eisenlohr - UNEMAT, Alta Floresta (pedro.eisenlohr@unemat.br)
#Based on Bauman et al. (2018), Clappe et al. (2018), Blanchet et al. (2008), Chavent et al. (2012) and others

# Acknowledgments:
# Sylvie Clappe
# Jo?o Carlos Pires de Oliveira

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
environment <- read.table(file.choose(),row.names=1,header=T,sep=",") #environmental variables
View(environment)
dim(environment)
#environment_std<-decostand(environment,"standardize") #caso queira padronizar a escala das vari?veis ambientais
#View(environment_std)

##Spatial matrix
ll<-read.table(file.choose(),row.names=1,header=T,sep=",") #longitude and latitude
View(ll)
dim(ll)




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


#----------------------------------------------------------------------#

######################################################################
#### HIERARCHICAL CLUSTERING OF VARIABLES (Chavent et al. 2012) ######
######################################################################
library(ClustOfVar)
tree <- hclustvar(environment)
plot(tree)
stab <- stability(tree,B=1000, graph = T) #To help in the selection of the number of partitions.
(nc<-which.max(stab$meanCR)+1)

### Here, you need to adjust he routine to the number of selected clusters (if you don't want to use nc).
P<-cutreevar(tree,nc,matsim=TRUE)
cluster <- P$cluster
X <- environment

for (i in 1:nc) {
  pc<-princomp(X[,which(cluster==i)],cor=TRUE)$sdev^2
  print(pc)
}
P$cluster
clusterID <- P$var
clusterID
write.table(P$scores,"PCAScores.csv")

env2 <- read.table("PCAScores.csv",row.names=1,header=T,sep=" ")
dim(env2)
View(env2)
vif(env2) #Be sure that no variable presents VIF>10.

#In case of no colinearities, please run the command below:
env3=env2

#In case of collinearities, please run the commands below (or re-define the number of groups below):
#(v1<-vifcor(env2,th=0.8)) 
#env3 <- exclude(env2, v1)

names(env3)
write.table(env3,"env_without_collinearities.csv")
environment = env3
View(environment)

#################################################################
#### END OF HIERARCHICAL CLUSTERING OF VARIABLES ################
#################################################################

#-----------------------------------------------------------------

#################################################################
### FORWARD SELECTION OF PREDICTORS
#################################################################

# Forward selection of the environmental variables
env.rda <- rda(spp.h,environment) #If you prefer to work with the whole response matrix, please change 'scores.pcoa' by 'spp.h'
env.rda
anova(env.rda, permutations = how(nperm=999))
### According to Blanchet et al. (2008): "If, and only if, the global 
### test is significant, one can proceed with forward selection"
(env.R2a <- RsquareAdj(env.rda)$adj.r.squared)
env.fwd <- forward.sel(spp.h, environment, R2thresh=env.R2a)
env.fwd #List of selected variables
env.sign <- sort(env.fwd$order)
env.red <- environment[,c(env.sign)]
head(env.red)
#vif(env.red)
#env.rda.selected <- rda(spp.h, env.red)
#env.R2a <- RsquareAdj(env.rda.selected)$adj.r.squared)
# save.image()

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
                    p.adjust = TRUE, MEM.all = FALSE, nperm = 999)) 
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
space.rda.selected <- rda(spp.h, spatial.red)
(env.R2a <- RsquareAdj(space.rda.selected)$adj.r.squared)

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
#mem.positive<-mem(candidates.sel,MEM.autocor = "positive")
#mem.positive
#class(mem.positive)
#dim(mem.positive)

save.image()



#############################################
########## VARIATION PARTITIONING ###########
#############################################
vprda <- varipart(spp.h, env.red, spatial.red) #classic variation partitioning
vprda

vprdaMSR <- msr.varipart(vprda, mem.all, nrepet = 999) #new variation partitioning (Clappe et al. 2018)
vprdaMSR #observe the result for fraction [a]. To obtain the significance of fraction [c], consider the result obtained by vprda.



####################################################
## Testing fraction [b] of variation partitioning ###
####################################################

#Testing fraction [b] (Bauman et al. 2019)
modsel.env <- listw.select(env.red, candidates, MEM.autocor = "positive", method = "global")
(global.env <- anova(rda(spp.h, env.red))$Pr[1])#apenas checar se foi significativo. N?o precisa anotar. Se n?o for significativo, n?o prossiga.
VP <- varpart(spp.h, env.red, spatial.red)
# SSEF:
(SSEF <- VP$part$indfract$Adj.R.square[2])
SSEF.test <- envspace.test(spe = spp.h, env = env.red, coord = ll, MEM.spe = spatial.red,
		listw.env = candidates[[modsel.env$best.id]], MEM.autocor = "positive", regular = FALSE, nperm = 999)
# Significance of the permutation test performed with MSR:
SSEF.test$pval


###########################################################################
####### Assessing the environmental model, after considering #########  
################# the effect of selected MEMs #############################
###########################################################################
env.rda<-rda(spp.h,env.red,spatial.red)
head(summary(env.rda)) ### Please observe the explanation of each axis.
#plot(env.rda)
spenvcor(env.rda) # species-environment correlation
intersetcor(env.rda) #"interset" correlation


###########################################################################
######## Assessing the spatial model, after considering ##############  
################# the effect of selected env ##############################
###########################################################################
spatial.rda<-rda(spp.h,spatial.red,env.red)
head(summary(spatial.rda)) ### Please observe the explanation of each axis.
#plot(spatial.rda)
spenvcor(spatial.rda) # species-space correlation
intersetcor(spatial.rda) #"interset" correlation


###########################################################################
#################### Evaluating the whole model ###########################
###########################################################################
all.predictors <- cbind(env.red,spatial.red)
all<-rda(spp.h,all.predictors)
anova(all)
head(summary(all)) ### Please observe the explanation of each axis.
#plot(all)
spenvcor(all) # species-predictors correlation
intersetcor(all) #"interset" correlation
# To test each axis individually:
rda.formula <- rda(spp.h~., data=all.predictors)
anova(rda.formula, by="axis")



# RDA's graphs
#all<-rda(spp.h,all.predictors)
FACTOR.df<-read.table(file.choose(),row.names=1,header=T,sep=",")
all.f<-cbind(FACTOR.df, all.predictors) ### all.predictors + FACTOR
dim(all.f)
View(all.f)
write.table(all.f,"All_butFACTOR.csv")
rda.graph<-all.f[,1]
View(rda.graph) #Factor
cca.results <- all
cca.summary <- summary(cca.results)
(imp.axis.1 <- cca.summary$cont$importance[2,1]) # explica??o do eixo 1
(imp.axis.2 <- cca.summary$cont$importance[2,2]) # explica??o do eixo 2
(rowScores <- as.data.frame(cca.results$CCA$u))  # scores das parcelas
write.table(rowScores,"scores_parcelas.csv")
(colScores <- as.data.frame(cca.results$CCA$v))  # scores das esp?cies 
write.table(colScores,"scores_spp.csv")
(bi.var <- as.data.frame(cca.results$CCA$biplot)) # scores das vari?veis ambientais
write.table(bi.var,"scores_amb.csv")
rowScores$FACTOR<-rda.graph
str(rowScores)

# tiff(filename="figura.tiff", res=600, height=600/72*600, width=600/72*600, compression= "lzw")



library(ggplot2)
ggrda1 <- ggplot() +
    scale_shape_manual(values=c(16,17,18,19,20,21,17,19,20)) + # os numeros sao os tipos de pontos
    geom_point(data=rowScores, aes(x=RDA1, y=RDA2, shape=FACTOR, colour=FACTOR)) +
  geom_vline(xintercept = 0, alpha=0.9, linetype = "dashed") +
  geom_hline(yintercept = 0, alpha=0.9, linetype = "dashed") +
  geom_point(data = rowScores, aes(x = RDA1, y = RDA2, shape=FACTOR, colour=FACTOR),
             size = 3,
             fill="grey50") +
  theme(legend.justification=c(0,0), legend.position=c(0.05,0.65),legend.key.size = unit(0.6, "cm"))+
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_line(colour = NULL),
        panel.grid.major = element_line(colour = NULL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.text.x = element_text(color ="black", size = 12, angle = 0),
        axis.text.y = element_text(color ="black", size = 12, angle = 0))
ggrda1

#ggrda1 <- ggplot() +
#     geom_point(data=rowScores, aes(x=RDA1, y=RDA2, shape=FACTOR, colour=FACTOR),
#                size=2,
#                fill = "white"
#                  ) +
#    geom_vline(xintercept = 0, alpha=0.9, linetype = "dashed") +
#    geom_hline(yintercept = 0, alpha=0.9, linetype = "dashed") +
#    geom_point(data = rowScores, aes(x = RDA1, y = RDA2, shape=FACTOR, colour=FACTOR),
#               size = 3,
#                fill="grey50") +
#    theme(legend.justification=c(0,0), legend.position=c(0.05,0.65),legend.key.size = unit(0.6, "cm"))+
#    theme(panel.background = element_rect(fill = "white", colour = NA),
#          panel.grid.minor = element_line(colour = NULL),
#          panel.grid.major = element_line(colour = NULL),
#          plot.background = element_rect(fill = "white", colour = NA),
#          panel.border = element_rect(fill = NA, colour = "black"),
#          axis.text.x = element_text(color ="black", size = 12, angle = 0),
#          axis.text.y = element_text(color ="black", size = 12, angle = 0))
#ggrda1


ggrda2 <- ggrda1 +
    xlab(paste("RDA 1 ", "(", round(imp.axis.1*100, digits = 1), "%)")) +
    theme(axis.title.x=element_text(angle = 0, size = 15)) + #, face = "bold"
    ylab(paste("RDA 2 ", "(", round(imp.axis.2*100, digits = 1), "%)")) +
    theme(axis.title.y=element_text(angle = 90, size = 15)) #, face = "bold"

# Visualizar
ggrda2




ggrda3 <- ggrda2 +
    geom_segment(data = bi.var, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                 size = 1,
                 arrow = arrow(length = unit(0.2, "cm"),
                               type = "closed", angle = 15),
                 color = "black",
                 alpha = 0.5) +
    geom_text(data = bi.var, aes(RDA1*1.05, RDA2*1.05,
                                 label = rownames(bi.var)),
              color = "black",
              size = 2,
              fontface = "bold")
#geom_text_repel(aes(x = colScores$PC1, y = colScores$PC2, label = rownames(colScores)))

# Visualizar
ggrda3



ggrda4 <- ggrda3 + ggsave("RDA_MOD.tiff",width = 15, height = 15, units = "cm")

# Visualizar
ggrda4


save.image()

#END
