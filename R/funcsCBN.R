###############################################################################
#' @title Produces \code{pat} files for H-CBN
#' 
#' @description \code{doPatForCBN} produces a binary matrix with dimensions 
#' patients x pathways, to be used as input by H-CBN. It can also write the 
#' matrix to a file.
#' 
#' @param groupsList list with as many elements as mutually exclusive groups,
#' each element consists of the genes assigned to the respective pathway
#' @param genesMat binary alteration matrix with rows representing patients 
#' and columns representing genes, where an 1 indicates that a gene has been
#' altered in one patient
#' @param write boolean, indicating whether to also write the output matrix to 
#' a file. Default to FALSE.
#' @param path path where the output matrix should be written. Default to
#' "~/Dropbox/PathTiMEx/RPathTiMEx/txts".
#' @param name name of the file to be written. Default to "outCBN".
#' 
#' @details The matrix returned by this function is the binary alteration matrix
#' on the basis of which order constraints between pathways are inferred. The
#' input to this function are the assignment of genes to pathways (as a list)
#' and the binary alteration matrix on the level of genes. A pathway is 
#' considered to be altered as soon as at least one of its gene members is 
#' altered. If \code{write=TRUE} and if \code{path} and \code{name} are 
#' provided, the matrix is written to the file \code{path/name.pat}. 
#' 
#' @return binary alteration matrix representing the alteration of pathways
#' in patients.
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases doPatForCBN
#' 
#' @export
doPatForCBN<-function(groupsList,genesMat,write,path,name)
{
    if (missing(write))
        write<-FALSE
    if (missing(path))
        path<-"~/Dropbox/PathTiMEx/RPathTiMEx/txts"
    if (missing(name))
        name<-"outCBN"
    
    noPaths<-length(groupsList)
    idxGenesUnused<-setdiff(c(1:dim(genesMat)[2]),unlist(groupsList))
    namesGenesUnused<-colnames(genesMat)[idxGenesUnused]
    
    pathsMat<-c()
    for (i in 1:noPaths)
    {
        genesNow<-groupsList[[i]]
        matGenesNow<-genesMat[,genesNow]
        if (length(genesNow)>1) # if the current pathway has more than one gene
        {
           mutInPath<-(apply(matGenesNow,1,sum)>0)+0
        }
        else
        {
            mutInPath<-matGenesNow
        }
        pathsMat<-cbind(pathsMat,mutInPath)
    }
    colnames(pathsMat)<-paste("path",c(1:noPaths),sep="")
    
    if (write)
    {
        write(c(dim(genesMat)[1],(noPaths+1)),file=paste(path,"/",name,".pat",sep=""))
        write.table(cbind(rep(1,dim(genesMat)[1]),pathsMat),file=paste(path,"/",name,".pat",sep=""),row.names=FALSE,col.names=FALSE,append=TRUE)
    }  
    
    return(pathsMat)
}



###############################################################################
#' @title Writes a poset to a file
#' 
#' @description \code{writePosetForCBN} writes the input poset to a file, in 
#' a form accepted as input by H-CBN.
#' 
#' @param poset the binary input matrix representing the order constraints
#' between events
#' @param path path where the output poset should be written. Default to
#' "~/Dropbox/PathTiMEx/RPathTiMEx/txts".
#' @param name name of the file to be written. Default to "outCBN".
#' 
#' @details The poset encodes order constraints between events and it is either
#' used as starting solution for H-CBN, or the final result of the H-CBN 
#' structure otpimization. The poset is written to file in a form compatible to 
#' H-CBN, namely on the first line the number of events, on each of the 
#' following lines a directed edge representing an order constraint, and on 
#' the last line the number 0. If \code{path} and \code{name} are provided, 
#' the poset is written to \code{path/name.poset}.
#' 
#' @return nothing
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases writePosetForCBN
#' 
#' @export
writePosetForCBN<-function(poset,path,name)
{
    if (missing(path))
        path<-"~/Dropbox/PathTiMEx/RPathTiMEx/txts"
    if (missing(name))
        name<-"outCBN"
        
    write(dim(poset)[1],file=paste(path,"/",name,".poset",sep=""))
    edges<-which(poset==1,arr.ind = TRUE)
    write.table(edges,file=paste(path,"/",name,".poset",sep=""),row.names = FALSE, col.names = FALSE, append = TRUE)
    write(0,file=paste(path,"/",name,".poset",sep=""),append = TRUE)
}



###############################################################################
#' @title Runs H-CBN
#' 
#' @description \code{runCBN} runs the H-CBN program from the the command line, 
#' for already existing .pat (binary alteration matrix) and .poset (initial
#' solution for structure) files.
#' 
#' @param path path where the folder containing the initial results is created, 
#' as well as where the final results are written. In the H-CBN help file, 
#' this is the path where the \code{filestem} directory is created.
#' @param name name of the directory where results are initially written. 
#' @param optionsSA string representing the command-line options of the
#' simmulated annealing optimization routine.
#' @param noThreads number of threads on which to run H-CBN. 
#' 
#' @details The folder were results of the H-CBN run are initially saved is 
#' \code{completeFilestem=path/name} (which, if not existent, is created). 
#' After running this function, five files are written: the optimized 
#' structure (poset), written to  \code{path/name/00000.poset} and also to 
#' \code{path/name.poset}, replacing the initial poset and becoming the initial
#' solution for the next iteration of the simmulated annealing optimization; 
#' the optimized lambda parameters, written to \code{path/name.lambda} and also
#' to \code{path/name/paramsEps.lambda}, where it also includes the estimate of
#' the error-rate epsilon; and the log of the simmulated annealing 
#' repetitions, which are written to \code{path/name.log}. The returned optimal
#' structure contains an additional dummy row and column identically set to 
#' 0, in preparation for the format required by the local group optimization.
#' 
#' @return list consisting of 
#' \itemize{
#' \item{\code{optStruct}} {a square matrix of dimension input number of 
#' pathways plus 1, representing the inferred optimal structure. The last row
#' and column correspond to a dummy pathway and are identically zero.}
#' \item{\code{eps}} {H-CBN estimate for the gene-wise error rate (positive real 
#' number)}
#' \item{\code{alpha}} {H-CBN estimate for alpha (positive real number).}
#' \item{\code{loglik}} {H-CBN estimate for the log likelihood.}
#' \item{\code{lamobs}} {H-CBN estimate for the waiting time rate of the 
#' observation time (fixed to 1).}
#' \item{\code{lams}} {vector, H-CBN estimate for the waiting time rates of the 
#' events (pathways).}
#' }
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases runCBN
#' 
#' @export
runCBN<-function(path,name,optionsSA,noThreads)
{
    # create directory where to write the poset
    completeFilestem<-paste(path,"/",name,sep="")
    dir.create(completeFilestem, showWarnings = FALSE)
    
    # no threads
    OMPthreads<-paste("OMP_NUM_THREADS=",noThreads,sep="")
    
    # 1. run CBN for the given .pat and .poset files with the specified simulated annealing options
    # write optimized structure to completeFilestem/00000.poset
    # writes optimized lambda parameters to path/name.lambda
    system(paste(OMPthreads, "h-cbn -f", completeFilestem, optionsSA, "-w"))
    
    # 2. read in and store the optimized structure form completeFilestem/00000.poset
    # add an extra NULL node, in preparation for the local group optimization
    optEdges<-scan(paste(completeFilestem,"/","00000.poset",sep=""))
    optEdges<-optEdges[-length(optEdges)]
    noPaths<-optEdges[1]
    optEdges<-optEdges[-1]
    optEdges<-matrix(optEdges,ncol=2,byrow = TRUE)
    optStruct<-matrix(0,nrow=(noPaths+1),ncol=(noPaths+1))
    optStruct[cbind(optEdges[,1],optEdges[,2])]<-1
    
    # 3. replace the initial solution for poset with the estimated one,
    # to be used as initial solution in the next iteration
    system(paste("mv ", completeFilestem, "/", "00000.poset", " ", completeFilestem, ".poset",sep=""))
    
    # 4. re-run again parameter estimation for the given poset for epsilon, as epsilon was not saved
    # save the parameters (eps, loglik, lambda_obs, lambdas) in the file completeFilestem/paramsEps.lambda
    system(paste(OMPthreads, " h-cbn -f ", completeFilestem, " -w -> ", completeFilestem, "/paramsEps.lambda",sep=""))
    
    # 5. store the estimated parameters
    params<-read.table(paste(completeFilestem,"/paramsEps.lambda",sep=""),sep="\t",header = TRUE)
    eps<-params[2]
    alpha<-params[3]
    loglik<-params[4]
    lamobs<-params[5]
    lams<-params[c(6:length(params))]
    
    result=list("optStruct"=optStruct,"eps"=eps,"alpha"=alpha,"loglik"=loglik,"lamobs"=lamobs,"lams"=lams)
    
    return(result)
}



###############################################################################
#' @title Transforms poset
#' 
#' @description \code{transformAndWritePoset} transforms the existing optimal
#' poset from the previous optimization step to a new compatible poset, to be
#' used as starting solution by the next run of H-CBN.
#' 
#' @param path path where the new compatible poset is written
#' @param name name of the poset file to be written
#' @param oldStruct the old poset to be transformed, encoding order relations
#' between pathways at the previous optimization step
#' @param newPathMat current binary alteration matrix at the level of pathways, 
#' with rows representing samples and columns pathways, where an 1 represents 
#' that a pathway is altered in a sample
#' 
#' @details The role of this function is to create a compatible starting 
#' poset for the new H-CBN optimization step, on the basis of the current 
#' optimal poset. This is done by simply removing any edge from between events 
#' which do not longer exist as pathways. The numbering of the evnets is kept 
#' as in the current poset, without checking for consistency between pathways 
#' between the old and current groupings. The new poset is written to the file 
#' path/name.poset.
#' 
#' @return nothing
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases transformAndWritePoset
#' 
#' @export
transformAndWritePoset<-function(path,name,oldStruct,newPathMat)
{
    # number of events for the new poset
    noEventsNow<-dim(newPathMat)[2]
    
    # if needed, remove edges between events which are no longer there
    #edgesOld<-which(oldStruct==1,arr.ind = TRUE)
    #allOldEvents<-unique(as.vector(edgesOld))
    #if (max(allOldEvents)>noEventsNow)
    #{
    #    oldStruct[c((noEventsNow+1):dim(oldStruct)[1]),]<-0
    #    oldStruct[,c((noEventsNow+1):dim(oldStruct)[1])]<-0        
    #}
    
    oldStruct<-make_linear_poset(p = noEventsNow)
    
    edgesOldToWrite<-which(oldStruct==1,arr.ind = TRUE)
    
    write(noEventsNow,file=paste(path,"/",name,".poset",sep=""))
    write.table(edgesOldToWrite,file=paste(path,"/",name,".poset",sep=""), row.names = FALSE, col.names = FALSE, append = TRUE)
    write(0,file=paste(path,"/",name,".poset",sep=""),append = TRUE)
}



# every time, a linear poset is used as the initial solution
transformAndWritePosetLinear<-function(path,name,oldStruct,newPathMat)
{
    noEventsNow<-dim(newPathMat)[2]
    
    oldStruct<-make_linear_poset(p = noEventsNow)
    
    edgesOldToWrite<-which(oldStruct==1,arr.ind = TRUE)
    
    write(noEventsNow,file=paste(path,"/",name,".poset",sep=""))
    write.table(edgesOldToWrite,file=paste(path,"/",name,".poset",sep=""), row.names = FALSE, col.names = FALSE, append = TRUE)
    write(0,file=paste(path,"/",name,".poset",sep=""),append = TRUE)
}



# the initially used recycling of the initial solution
transformAndWritePosetInitial<-function(path,name,oldStruct,newPathMat)
{
    # number of events for the new poset
    noEventsNow<-dim(newPathMat)[2]
    
    # if needed, remove edges between events which are no longer there
    edgesOld<-which(oldStruct==1,arr.ind = TRUE)
    allOldEvents<-unique(as.vector(edgesOld))
    if (max(allOldEvents)>noEventsNow)
    {
        oldStruct[c((noEventsNow+1):dim(oldStruct)[1]),]<-0
        oldStruct[,c((noEventsNow+1):dim(oldStruct)[1])]<-0        
    }
    
    edgesOldToWrite<-which(oldStruct==1,arr.ind = TRUE)
    
    write(noEventsNow,file=paste(path,"/",name,".poset",sep=""))
    write.table(edgesOldToWrite,file=paste(path,"/",name,".poset",sep=""), row.names = FALSE, col.names = FALSE, append = TRUE)
    write(0,file=paste(path,"/",name,".poset",sep=""),append = TRUE)
}


transformAndWritePoset2<-function(path,name,oldStruct,oldGroupys,newGroupys)
{
    oldGroupys<-lapply(oldGroupys,function(x){sort(x)})
    newGroupys<-lapply(newGroupys,function(x){sort(x)})
    #oldGroupys%in%newGroupys
    
    # number of events for the new poset
    noEventsNow<-dim(newPathMat)[2]
    
    # if needed, remove edges between events which are no longer there
    #edgesOld<-which(oldStruct==1,arr.ind = TRUE)
    #allOldEvents<-unique(as.vector(edgesOld))
    #if (max(allOldEvents)>noEventsNow)
    #{
    #    oldStruct[c((noEventsNow+1):dim(oldStruct)[1]),]<-0
    #    oldStruct[,c((noEventsNow+1):dim(oldStruct)[1])]<-0        
    #}
    
    oldStruct<-make_linear_poset(p = noEventsNow)
    
    edgesOldToWrite<-which(oldStruct==1,arr.ind = TRUE)
    
    write(noEventsNow,file=paste(path,"/",name,".poset",sep=""))
    write.table(edgesOldToWrite,file=paste(path,"/",name,".poset",sep=""), row.names = FALSE, col.names = FALSE, append = TRUE)
    write(0,file=paste(path,"/",name,".poset",sep=""),append = TRUE)
}



###############################################################################
#' @title Compares two posets
#' 
#' @description \code{compareStructs} compares two posets of same size, 
#' encoded as binary matrices. It returns various metrics comparing the two
#' matrices. 
#' 
#' @param oldStruct the first of the two structures to compare, corresponding to
#' the previous optimization step.
#' @param newStruct the second of the two structures to compare, corresponding
#' to the current optimization step.
#' 
#' @details For computing the false posive and false negative rates, it is 
#' assumed that the first argument (oldStruct) is the real structure. In this
#' help page, fhe first structure is referred to as the old structure, while 
#' the new structure is referred to as the new structure. The number of all 
#' possible directed edges, often used to compute rates, equals twice the 
#' number of events choose 2. In case the two structures are equal, the 
#' \code{accuracy} measure is 1.
#' 
#' @return list consisting of 
#' \itemize{
#' \item{\code{falseNeg}} {number of edges which are in the old structure,
#' but not in the new structure, divided by the number of all possible directed
#' edges.}
#' \item{\code{falsePos}} {number of edges which are in the new structure,
#' but not in the old structure, divided by the number of all possible directed
#' edges.}
#' \item{\code{accuracy}} {measure of similarity of the two structures, equal
#' to 1 - number of mismatches divided by the number of all possible directed 
#' edges. In case the two structures are equal, \code{accuracy = 1}.}
#' \item{\code{noMismatches}} {the number of edges which are found in only one
#' of either of the two structures}
#' \item{\code{densOld}} {number of edges in the old structure.}
#' \item{\code{densNew}} {number of edges in the new structure.}
#' }
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases compareStructs
#' 
#' @export
compareStructs<-function(oldStruct,newStruct)
{ 
    sizeMat<-dim(oldStruct)[1]
    
    noMismatches<-sum(abs(oldStruct-newStruct))
    noTotalEdges<-sizeMat^2-sizeMat
    accuracy<-1-noMismatches/noTotalEdges
    
    falseNeg<-length(which(oldStruct-newStruct==1))/noTotalEdges
    falsePos<-length(which(newStruct-oldStruct==1))/noTotalEdges
    
    densOld<-sum(oldStruct)
    densNew<-sum(newStruct)
    
    measures<-list("falseNeg" = falseNeg,"falsePos" = falsePos,
                   "accuracy" = accuracy, "noMismatches" = noMismatches,
                   "densOld" = densOld, "densNew" = densNew)
    return(measures)
}



# function to map two CBN structures of the same size, for comparison
# code is just copied from the main text and wrapped as a function
mapStructures<-function(vecAssignPrev,vecAssignNow,optStructHere)
{
    uniqueClust<-unique(vecAssignPrev)
    clustMapped<-c()
    for (idxC in uniqueClust)
    {
        membersOld<-which(vecAssignPrev==idxC)
        clustMapped[idxC]<-unique(vecAssignNow[membersOld]) # the new indices of clusters
    }
    
    # the adjacency matrix of the old structure, mapped to be compared to the new structure
    oldStruct<-optStructHere[match(c(1:length(uniqueClust)),clustMapped),match(c(1:length(uniqueClust)),clustMapped)]
    
    # add a coat of zeroes to both rows and columns, for comparison
    oldStruct<-rbind(oldStruct,rep(0,length(uniqueClust)))
    oldStruct<-cbind(oldStruct,rep(0,(length(uniqueClust)+1)))
    
    return(oldStruct)
}
