# Author: Simona Constantinescu; simona.constantinescu@bsse.ethz.ch

###############################################################################
#' @title Summarizes the inferred groups
#' 
#' @description \code{addPathBinary} summarizes the assignment of genes to 
#' pathways as a binary matrix.
#' 
#' @param resIterBinTiMEx structure with the results of running TiMEx 
#' iteratively on any input binary matrix, as returned by the function
#' \code{iterTiMExBinaryMats}.
#' 
#' @details This function appends to the input structure the field 
#' \code{genesBinary}, representing the assignment of the identified genes to 
#' pathway, in the form of a binary matrix.
#' 
#' @return resIterBinTiMEx, the input structure, with the added field 
#' \code{genesBinary}. \code{genesBinary} is a matrix with rows the 
#' identified mutually exclusive pathways and columns the genes. Each entry 
#' [i,j] = 1 if gene j was identified in pathway i and 0 otherwise.
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases addPathBinary
#' 
#' @export
addPathBinary<-function(resIterBinTiMEx)
{
    genesBinary<-list()
    noGenes<-dim(resIterBinTiMEx$mat)[2]
    noPaths<-length(resIterBinTiMEx$geneIndices)
    if (noPaths>0)
    {
        for (i in 1:noPaths)
        {
            genesBinary[[i]]<-rep(0,noGenes)
            genesBinary[[i]][resIterBinTiMEx$geneIndices[[i]]]<-1
        }
        genesBinary<-matrix(unlist(genesBinary),ncol=noGenes,byrow = TRUE)
    }
    
    resIterBinTiMEx$genesBinary<-genesBinary
    return(resIterBinTiMEx)   
}



###############################################################################
#' @title Counts duplicates in a data frame
#' 
#' @description \code{count.duplicates} counts how many times each row appears
#' in the input data frame
#' 
#' @param DF the data frame whose rows will be counted
#' 
#' @details This function returns a new data frame, with the unique rows of the
#' input one, and an additional column, counting how many times each row 
#' appears. 
#' 
#' @return a data frame with unique rows, followed by an additional column
#' which counts the number of times each row was found in the input data
#' frame.
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases count.duplicates
#' 
#' @export
count.duplicates <- function(DF){
    x <- do.call('paste', c(DF, sep = '\r'))
    ox <- order(x)
    rl <- rle(x[ox])
    p<-cbind(DF[ox[cumsum(rl$lengths)],,drop=FALSE],count = rl$lengths)
    p<-p[order(p$count,decreasing = TRUE),]
    return(p)
}



###############################################################################
#' @title Assesses how well each group was reconstructed
#' 
#' @description \code{countPerfSims} assesses the reconstruction performance
#' of each group, inferred from any binary matrix. 
#' 
#' @param pathwaysTable a binary matrix of all groups identified aross many
#' repetitions, as returned by multiple runs of \code{addPathBinary}. rows 
#' represent indetified groups and columns represent genes. the element [i,j] 
#' is 1 if gene j was identified as part of pathway i and 0 otherwise.
#' @param groups a list with the real assignment of genes to pathways. each 
#' element of the list is a vector with the genes belonging to that pathway.
#' 
#' @details This function summarizes multiple repetitions of the inference
#' of pathways from a simulated dataset. Ideally, the groups with largest 
#' counts are the real groups.
#' 
#' @return a binary matrix in which the first columns represent the unique 
#' pathways identified in the input matrix (columns named with \code{Gr}, 
#' followed by the group index, followed by \code{ge}, followed by the gene
#' index). next column, \code{counts}, counts how many times each group was 
#' identified as mutually exclusive. next columns correspond to each of the 
#' unique pathways identified, and show either the total number of genes form 
#' each group which were part of the group (\code{Tot}), or how much this 
#' number represents, as a percentage, from the real number of genes in each 
#' group (\code{Per}).
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases countPerfSims
#' 
#' @export
countPerfSims<-function(pathwaysTable,groups)
{
    pathwaysTable<-as.data.frame(pathwaysTable)
    pathwaysTable<-count.duplicates(pathwaysTable)
    pathwaysSimple<-pathwaysTable[,-dim(pathwaysTable)[2]]
    noGroups<-length(groups)
    # order the identified pathways by the number of elements
    #pathwaysTable<-pathwaysTable[do.call(order, as.data.frame(pathwaysSimple)),]
    
    allNoEls<-apply(pathwaysTable,1,function(x){
        numberEls<-c() #number of elements identified in each group
        for (i in 1:noGroups)
        {
            numberEls[[i]]<-length(which(x[groups[[i]]]==1))
        }
        return(numberEls)
    })
    
    percentNoEls<-apply(pathwaysTable,1,function(x){
        percentEls<-c() #number of elements identified in each group
        for (i in 1:noGroups)
        {
            percentEls[[i]]<-length(which(x[groups[[i]]]==1))/length(groups[[i]])
        }
        return(percentEls)
    })
    
    pathwaysTable<-cbind(pathwaysTable,t(allNoEls),t(percentNoEls))
    colnames(pathwaysTable)<-c(unlist(lapply(c(1:length(groups)),function(x){
        paste("Gr",x,"ge",groups[[x]],sep="")})),"counts",
        unlist(lapply(c(1:length(groups)),function(x){paste("gr",x,"Tot",sep="")})),
        unlist(lapply(c(1:length(groups)),function(x){paste("gr",x,"Per",sep="")})))
    
    pathwaysTable<-pathwaysTable[order(pathwaysTable$counts,decreasing = TRUE),]
    
    return(pathwaysTable)
}


# function to summarize how well each group was reconstructed
# inputs:
#   pathwaysTable: matrix with rows the indentified groups, together with the absolute and relative values of how well each group was reconstructed
#   groups: list with the original groups
#   noReps: number of repetitions
# outputs:
#   summaryPerf: list with as many elements as identified groups, and for each group a matrix with percentage of reconstruction and counts
###############################################################################
#' @title Assesses how well each group was reconstructed
#' 
#' @description \code{summarizeSims} assesses the reconstruction performance
#' of each group, inferred from any binary matrix. 
#' 
#' @param pathwaysTableSum a binary matrix of all unique groups identified 
#' across many repetitions, together with summary information about each of the 
#' identified groups, as returned by the function \code{countPerfSims}.
#' @param groups a list with the real assignment of genes to pathways. each 
#' element of the list is a vector with the genes belonging to that pathway.
#' @param noReps the number of repetitions of the simulation experiments
#' 
#' @details This function further summarizes multiple repetitions of the 
#' inference of pathways from a simulated dataset, having as input a summary
#' table with the unique groups inferred across many runs, as returned by the
#' function \code{countPerfSims}. The summary returned by this function doesn't
#' provide an assessment of the false positive rate.
#' 
#' @return a list with as many elements as real pathways. each element is
#' a matrix in which the first column is the reconstruction percentage of
#' each pathway, and the second element is the number of times in which that
#' reconstruction was attained.
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases summarizeSims
#' 
#' @export
summarizeSims<-function(pathwaysTableSum,groups,noReps)
{
    noGroups<-length(groups)
    summaryPerf<-list()
    # for each of the identified groups
    for (i in 1:noGroups) 
    {
        # the position of the current group
        pos<-which(colnames(pathwaysTableSum)==paste("gr",i,"Per",sep="")) 
        summaryPerf[[i]]<-matrix(0,nrow=(length(unique(pathwaysTableSum[,pos]))),ncol=2)
        c<-1
        for (j in setdiff(unique(pathwaysTableSum[[pos]]),0))
        {
            smallMat<-pathwaysTableSum[which(pathwaysTableSum[,pos]==j),]
            summaryPerf[[i]][c,]<-c(j,sum(smallMat$counts))
            c<-c+1
        }
        summaryPerf[[i]][dim(summaryPerf[[i]])[1],]<-c(0,(noReps-sum(summaryPerf[[i]][,2])))
        colnames(summaryPerf[[i]])<-c("PerBack","SumCounts")
        summaryPerf[[i]]<-as.data.frame(summaryPerf[[i]])
    }
    return(summaryPerf)
}



# function to join multiple experiments and prepare matrix for ggplot2
prepareMatGGplot<-function(listPlotPerf)
{
    noExps<-length(listPlotPerf) # number of different scenarios
    newMat<-c()
    k<-0
    for (i in 1:noExps) # for each scenario
    {
        for (j in 1:length(listPlotPerf[[i]])) # for each group
        {
            k<-k+1
            newMat<-rbind(newMat,cbind(listPlotPerf[[i]][[j]],rep(k,dim(listPlotPerf[[i]][[j]])[1])))
        }
    }
    colnames(newMat)[3]<-"Group"
    return(newMat)
}
