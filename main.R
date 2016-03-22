### main function of pathTiMEx
rm(list=ls())
load("~/Dropbox/PathTiMEx/RPathTiMEx/datas/colonMat.RData")

source("~/Dropbox/PathTiMEx/RPathTiMEx/pathTiMEx/R/libraries.R")
source("~/Dropbox/PathTiMEx/RPathTiMEx/pathTiMEx/R/funcsSearch.R")
source("~/Dropbox/PathTiMEx/RPathTiMEx/pathTiMEx/R/funcsCBN.R")
source("~/Dropbox/PathTiMEx/RPathTiMEx/pathTiMEx/R/funcsPoset.R")
source("~/Dropbox/PathTiMEx/RPathTiMEx/pathTiMEx/R/funcsTiMEx.R")
source("~/Dropbox/PathTiMEx/RPathTiMEx/pathTiMEx/R/funcsSims.R")
source("~/Dropbox/PathTiMEx/RPathTiMEx/pathTiMEx/R/funcsOth.R")
source("~/Dropbox/PathTiMEx/RPathTiMEx/pathTiMEx/R/funcsPathTiMEx.R")

## run TiMEx iteratively on a binary matrix (in this example on a colon cancer dataset: colonMat)
colonMatNew<-doMetagene(colonMat)$newMat
colonMetagroups<-doMetagene(colonMat)$groups
initialGroupsStruct<-iterTiMExBinaryMats(colonMatNew,pairMu = 0.5,pairPvalue = 0.01,groupPvalue = 0.1)

## arguments to the pathTiMEx function
Datamat<-initialGroupsStruct$mat # binary alteration matrix of genes
path<-"~/Dropbox/PathTiMEx/CBNs/" # path where the CBN folder will be created
name<-"Test" # name of the CBN folder

noReps<-100 # number of times to run the joint optimization
optionsSA<-"-s -T 1 -N 200" # options for the simulated annealing
noThreads<-4 # number of threads to use for running CBNs
numsave<-2000 # how many scores to save for plotting
skipsteps<-100 # once in how many steps to save the score
gamma<-0.1
noNoOpen<-1 # number of iterations in the beginning in which to allow the opening of a new mutually exclusive group
additionalGenes<-setdiff(c(1:dim(Datamat)[2]),unlist(initialGroupsStruct$geneIndices))
out<-pathTiMEx(initialGroupsStruct,Datamat,path,name,noReps,optionsSA,noThreads,numsave,skipsteps,gamma,noNoOpen,additionalGenes)
