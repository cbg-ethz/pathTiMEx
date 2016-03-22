# Author: Simona Constantinescu; simona.constantinescu@bsse.ethz.ch

###############################################################################
#' @title Interative TiMEx on any binary dataset
#' 
#' @description \code{iterTiMExBinaryMats} performs TiMEx iteratively on any 
#' binary dataset, by each time removing the most significant group of largest 
#' size and runnning TiMEx again.
#' 
#' @param mat binary input matrix, with rows representing samples and columns
#' representing genes
#' @param pairMu threshold on the pair-level mu (input to the function 
#' \code{TiMEx}). Default to 0.5.
#' @param pairPvalue threshold on the pair-level pvalue (input to the 
#' function \code{TiMEx}). Default to 0.01.
#' @param groupPvalue threshold on the group-level pvalue (input to the 
#' function \code{TiMEx}). Default to 0.1.
#' @param corr type of correction employed by TiMEx. It has to either be "bonf"
#' or "fdr". Default to "bonf".
#' @param noRuns maximum number of times to run the iterative procedure. 
#' Default to 10.
#' 
#' @details This function reports the most significant mutually exclusive 
#' groups of largest size inferred iteratively, on the basis of any binary input
#' matrix. Once a group is reported, it is removed from the dataset and TiMEx is 
#' ran again. The groups are reported decreasingly by size.
#' 
#' @return list consisting of 
#' \itemize{
#' \item{\code{geneNames}} {list consisting of the names of the genes part
#' of the identified mutually exclusive groups. each element is a vector of 
#' gene names corresponding to one group.}
#' \item{\code{geneIndices}} {list consisting of the indices in the input
#' binary matrix of the genes part of the identified mutually exclusive 
#' groups. each element is a vector of gene indices corresponding to one 
#' group.}
#' \item{\code{mat}} {input binary matrix}
#' \item{\code{samplesModel}} {input structure containing the simulated 
#' pathTiMEx model, including the binary matrix used as input to the function
#' \code{TiMEx}.}
#' \item{\code{muEst}} {vector with the estimated mu for the identified 
#' mutually exclusive groups, of length the number of groups identified.}
#' \item{\code{corr}} {the input method for multiple testing correction used
#' by TiMEx.}
#' }
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases iterTiMExBinaryMats
#' 
#' @export
iterTiMExBinaryMats<-function(mat,pairMu,pairPvalue,groupPvalue,corr,noRuns)
{   
    if (missing(pairMu))
        pairMu<-0.5
    if (missing(pairPvalue))
        pairPvalue<-0.01
    if (missing(groupPvalue))
        groupPvalue<-0.1
    if (missing(noRuns))
        noRuns<-10
    if (missing(corr))
        corr<-"bonf"
    
    # record the initial matrix and gene names for reconstructing the indices
    # afterwards
    initialMat<-mat
    initialNames<-colnames(mat)
    
    geneNames<-list()
    muEst<-c()
    i<-1
    maxSize<-1
    
    while (i<=noRuns & !(is.null(dim(mat)) || ncol(mat)==0) & maxSize>0)
    {
        print(paste("TiMEx iteration number: ", i, sep=""))
        resTiMEx<-TiMEx(mat, pairMu, pairPvalue, groupPvalue)
        maxSize<-length(resTiMEx$idxSignif)
        
        if (maxSize>0)
        {
            idxName<-which(names(resTiMEx$idxSignif[[maxSize]])==corr)
            l<-length(resTiMEx$idxSignif[[maxSize]][[idxName]])
            
            while (l==0 && (maxSize>0))
            {
                resTiMEx$idxSignif[[maxSize]]<-NULL
                maxSize<-length(resTiMEx$idxSignif)
                if (maxSize>0)
                {
                    idxName<-which(names(resTiMEx$idxSignif[[maxSize]])==corr)
                    if (length(idxName)>0)
                        l<-length(resTiMEx$idxSignif[[maxSize]][[idxName]])
                    else
                        l<-0
                }
                else
                    l<-0
            }
            if (maxSize>0)
            {
                d<-dim(resTiMEx$idxSignif[[maxSize]][[idxName]])[1]
                while (maxSize>0 &&  !is.null(d) && d==0)
                {
                    maxSize<-maxSize-1
                    idxName<-which(names(resTiMEx$idxSignif[[maxSize]])==corr)
                    d<-dim(resTiMEx$idxSignif[[maxSize]][[idxName]])[1]
                }
                
                # if more than one group of maximal size were found
                if (!is.null(dim(resTiMEx$idxSignif[[maxSize]][[idxName]])))
                {
                    # keep the most significant group of largest size
                    geneNames[[i]]<-unlist(resTiMEx$genesSignif[[maxSize]][[idxName]][1,]) 
                    muEst[i]<-resTiMEx$MusGroup[[maxSize]][[idxName]][1]
                    toRemove<-unlist(resTiMEx$idxSignif[[maxSize]][[idxName]][1,])
                } else # only one group of maximal size was found
                {
                    geneNames[[i]]<-unlist(resTiMEx$genesSignif[[maxSize]][[idxName]])
                    muEst[i]<-resTiMEx$MusGroup[[maxSize]][[idxName]]
                    toRemove<-unlist(resTiMEx$idxSignif[[maxSize]][[idxName]])
                }
                # remove group from the data before runnign TiMEx again
                mat<-mat[,-toRemove]
                i<-i+1
            }
        }
    }
    
    # arrange the groups by size decreasingly
    geneNames<-geneNames[order(sapply(geneNames,length),decreasing = TRUE)]
    
    # retrieve the initial indices
    geneIndices<-lapply(geneNames,function(x){
        match(x,initialNames)
    })
    
    # create return structure
    l<-list("geneNames" = geneNames, "geneIndices" = geneIndices, 
            "mat" = initialMat, "muEst" = muEst,"corr" = corr)
    
    return(l)  
}



###############################################################################
#' @title Add noise to a binary matrix
#' 
#' @description \code{addNoise} adds binomially sampled noise to a binary 
#' matrix, with success probability epsilon.
#' 
#' @param mat binary input matrix, with rows representing samples and columns
#' representing genes
#' @param epsilon the error rate at the level of genes, real number 
#' between 0 and 1
#' 
#' @details Any binary matrix can be used as input.
#' 
#' @return The input binary matrix with added noise.
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases addNoise
#' 
#' @export
addNoise<-function(mat,epsilon)
{
    # add noise to the genes with probability epsilon
    noPats<-dim(mat)[1]
    for (i in 1:noPats) # for each patient
    {
        # the genes to which noise is added, sampled binomially with success 
        # probability epsilon (1 turns into 0 and 0 turns into 1)
        genesNoise<-which(rbinom(prob=epsilon,n=dim(mat)[2],size=1)==1)
        # if there is at least one gene for which noise is added
        if (length(genesNoise)>0) 
        {
            mat[i,genesNoise]<-(!mat[i,genesNoise])+0
        }
    }
    return(mat)
}



###############################################################################
# function to include additional genes in the input group structure of 
# pathTiMEx
# inputs: initialGroupsStruct (input structure resulting after iteratively
# running pathTiMEx on an input dataset, as resulting from the function ...)
# additionalGenes: vector of indices (in the Datamat matrix) of additional 
# genes to include
# Datamat: binary input matrix patients x genes, which will be used throughout
# the optimization routine. Default to the matrix used to run TiMEx
createInitialGroups<-function(initialGroupsStruct,additionalGenes,Datamat)
{
    if (missing(Datamat))
        Datamat<-initialGroupsStruct$mat
    
    groupReturned<-initialGroupsStruct$geneIndices
    genesToAdd<-intersect(additionalGenes,c(1:dim(Datamat)[2]))
    genesToAdd<-setdiff(genesToAdd,unlist(groupReturned))
    
    if (length(genesToAdd)>20)
    {
        print("Too many genes to add; please update the list!")
        return(groupReturned)
    }
    
    if (length(genesToAdd)==0)
    {
        print("No other genes added")
        return(groupReturned)
    }
    
    groupReturned<-append(groupReturned,as.list(genesToAdd))
    return(groupReturned)
    
}





# plots a mutually exclusive group, customized
plotMEGroup<-function(Datamat,groupysInput,colsBroman,pathPlot,namesOfGenes)
{
    colGroups<-colsBroman[c(1:length(groupysInput))]
    pdf("~/Desktop/ee.pdf",width=8,height=4)
    par(mfcol=c(2,3))
    for (idx in 1:length(groupysInput))
    {
        group<-groupysInput[[idx]]
        l<-length(group)
        submatrix<-Datamat[,group]
        if (l>1)
        {
            freqs<-round(apply(submatrix,2,sum)*100/dim(submatrix)[1],2)
            ord<-order(freqs,decreasing=TRUE)
            submatrix<-submatrix[,ord]
            mat2<-submatrix
            for (i in l:1)
                mat2<-mat2[order(mat2[,i]),]
            mat2<-mat2[nrow(mat2):1,]
            mat2<-mat2[,ncol(mat2):1]
        }
        if (l==1)
        {
            mat2<-submatrix
            freqs<-sum(submatrix)/length(submatrix)*100
            mat2<-cbind(submatrix,rep(0,length(mat2)))
            for (i in 2:1)
                mat2<-mat2[order(mat2[,i]),]
            mat2<-mat2[nrow(mat2):1,]
            mat2<-mat2[,ncol(mat2):1]
            ord<-c(1)
            
        }
                
        image(mat2,col=c("white",colGroups[idx]),yaxt="n",xaxt='n',fg=colGroups[idx])
        axis(2,at=seq(from=0,to=1,length.out=l),las=2,tck=0,cex.axis=0.7,labels=colnames(Datamat)[groupysInput[[idx]]][rev(ord)])
    }
    dev.off()
    
}
