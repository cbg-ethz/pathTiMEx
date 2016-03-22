# Author: Simona Constantinescu; simona.constantinescu@bsse.ethz.ch

###############################################################################
#' @title Randomly generates a pathTiMEx model
#' 
#' @description \code{generateRandomModel} generates randomly the structure 
#' and the parameters of a pathTiMEx model, including the exponential rates 
#' of the waiting times to mutations, the intensities of mutual exclusivity
#' of the pathways, the error rate, the poset and the pathway assignments. 
#' 
#' @param n number of genes
#' @param k number of pathways
#' @param lamobs exponential rate of the observation time (usually fixed
#' to 1)
#' @param minLam lower limit for simulating the exponential rates of the
#' waiting times to mutations of the genes
#' @param maxLam upper limit for simulating the exponential rates of the
#' waiting times to mutations of the genes
#' @param muInterval vector of values from which the intensities of mutual
#' exclusivity of the pathways are uniformly sampled
#' @param ldensesInterval vector of values from which the density of edges in 
#' the poset is uniformly sampled
#' @param epsilon the error rate at the level of genes, real number between 0 
#' and 1
#' 
#' @details The exponential rate of the observation time is kept fixed (usually
#' to 1) when generating parameters for the model.
#' 
#' @return lists of parameters and structure for a pathTiMEx model, as follows:
#' \itemize{
#' \item{\code{lams}} {vector of exponential rates (lambda) of the waiting
#' times to mutations of the genes}
#' \item{\code{lamobs}} {the exponential rate of the observation time, usually
#' fixed to 1}
#' \item{\code{mus}} {vector of intensities of mutual exclusivity (mu) for 
#' each pathway}
#' \item{\code{ldense}} {the density of edges in the poset}
#' \item{\code{poset}} {the poset of the model}
#' \item{\code{assignment}} {the assignment of genes to pathways}
#' \item{\code{epsilon}} {the error rate at the level of genes, real number 
#' between 0 and 1}
#' \item{\code{type}} {way in which the model was simulated, in this case 
#' \code{random}}
#' }
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases generateRandomModel
#' 
#' @export
generateRandomModel<-function(n,k,lamobs,minLam,maxLam,muInterval,ldensesInterval,epsilon)
{
    # exponential rates for the waiting times of the genes
    lams<-runif(n,min=minLam,max=maxLam)
    # mutual exclusivity intensities of the pathways
    mus<-sample(muInterval,k,replace=TRUE)
    
    # density of edges
    ldense<-sample(ldensesInterval,1,replace=TRUE)
    poset<-genRandomPoset(k, ldense)
    
    # assignments of the n genes to the k pathways; repeat until each pathway has at least one member
    assignment<-sample(c(1:k),n,replace=TRUE)
    while (length(unique(assignment))!=k)
        assignment<-sample(c(1:k),n,replace=TRUE)
    names(assignment)<-c(1:n)
    
    l<-list("lams" = lams, "lamobs" = lamobs, "mus" = mus, "ldense" = ldense, 
            "poset" = poset, "assignment" = assignment, "epsilon"=epsilon, "type" = "random")
    
    return(l)
}



###############################################################################
#' @title Generates a pathTiMEx model for given parameters
#' 
#' @description \code{generateStructureModel} randomly generates the
#' structure for a pathTiMEx model under given parameters, i.e. the 
#' poset and the pathway assignments are generated, while the 
#' exponential rates of the waiting times to mutations 
#' and of the observation time, the intensities of mutual exclusivity and the 
#' error rate are known.
#' 
#' @param lams vector of exponential rates of the waiting times to mutations, 
#' equal in size to the number of genes
#' @param lamobs the exponential rate of the observation time (usually fixed
#' to 1)
#' @param mus vector of intensities of mutual exclusivity, equal in size to
#' the number of pathways
#' @param ldensesInterval vector of values from which the density of edges in 
#' the poset is uniformly sampled
#' @param epsilon the error rate at the level of genes, real number 
#' between 0 and 1
#' 
#' @details The exponential rate of the observation time is kept fixed (usually
#' to 1) when generating parameters for the model.
#' 
#' @return lists of parameters and structure for a pathTiMEx model, as follows:
#' \itemize{
#' \item{\code{lams}} {input vector of exponential rates of the waiting times 
#' to mutation of each gene}
#' \item{\code{lamobs}} {input exponential rate of the observation time}
#' \item{\code{mus}} {input vector of mutual exclusivity intensities of each 
#' pathway}
#' \item{\code{ldense}} {the density of edges in the poset}
#' \item{\code{poset}} {the poset of the model}
#' \item{\code{assignment}} {the assignment of genes to pathways}
#' \item{\code{epsilon}} {the error rate at the level of genes, real number 
#' between 0 and 1}
#' \item{\code{type}} {way in which the model was simulated, in this case 
#' \code{structure}}
#' }
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases generateStructureModel
#' 
#' @export
generateStructureModel<-function(lams,lamobs,mus,epsilon)
{    
    # density of edges
    ldense <- sample(ldensesInterval,1,replace=TRUE)
    poset <- genRandomPoset(length(mus), ldense)
    
    # assignments of the n genes to the k pathways; repeat until each pathway has at least one member
    assignment <- sample(c(1:length(mus)),length(lams),replace=TRUE)
    while (length(unique(assignment))!=length(mus))
        assignment <- sample(c(1:length(mus)),length(lams),replace=TRUE)
    names(assignment) <- c(1:length(lams))
    
    l<-list("lams" = lams, "lamobs" = lamobs, "mus" = mus, "ldense" = ldense, 
            "poset" = poset, "assignment" = assignment, "epsilon"= epsilon, 
            "type" = "structure")
    
    return(l)
}



###############################################################################
#' @title Generates a pathTiMEx model for given structure and parameters
#' 
#' @description \code{generateFixedModel} generates a pathTiMEx model 
#' for given structure (poset and pathway assignment) and parameters
#' 
#' @param lams the exponential rates of the waiting times to mutations
#' @param lamobs the exponential rate of the observation time (usually fixed
#' to 1)
#' @param mus the intensities of mutual exclusivity
#' @param listEdges list of pairs [i,j] such that the element [i,j] in the 
#' poset equals 1
#' @param assignment vector of assignments of genes to pathways
#' @param epsilon {the error rate at the level of genes, real number 
#' between 0 and 1}
#' 
#' @details The exponential rate of the observation time is kept fixed (usually
#' to 1) when generating parameters for the model
#' 
#' @return lists of parameters and structure for a pathTiMEx model, as follows:
#' \itemize{
#' \item{\code{lams}} {input vector of exponential rates of the waiting times 
#' to mutation of each gene}
#' \item{\code{mus}} {input vector of mutual exclusivity intensities of each 
#' pathway}
#' \item{\code{poset}} {the poset of the model}
#' \item{\code{assignment}} {the assignment of genes to pathways}
#' \item{\code{epsilon}} {the error rate at the level of genes, real number 
#' between 0 and 1}
#' \item{\code{type}} {way in which the model was simulated, in this case 
#' \code{fixed}}
#' }
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases generateFixedModel
#' 
#' @export
generateFixedModel<-function(lams,lamobs,mus,listEdges,assignment,epsilon)
{    
    # all the elements in listEdges are an 1 in the poset input matrix
    poset<-matrix(0,nrow=length(mus),ncol=length(mus))
    mats<-lapply(listEdges,function(x){poset[cbind(x[1],x[2])]<-1; return(poset)})
    poset<-(Reduce("+",mats)>0)+0 
    
    l<-list("lams" = lams, "lamobs" = lamobs, "mus" = mus, "poset" = poset, 
            "assignment" = assignment, "epsilon" = epsilon, "type" = "fixed")
    
    return(l)
}



###############################################################################
#' @title Draws samples from a pathTiMEx model
#' 
#' @description \code{drawSamples} generates samples form a pathTiMEx model
#' 
#' @param m number of samples to be simulated
#' @param model the pathTiMEx model, as a list, as returned by the functions 
#' which generate data: \code{generateRandomModel}, \code{generateFixedModel} 
#' and \code{generateStructureModel}
#' 
#' @details This function generates a binary matrix on the level of genes 
#' from the pathTiMEx model, as well as various other features of the model, 
#' such as the waiting times to alteration, the observation times etc.
#' 
#' @return list consisting of
#' \itemize{
#' \item{\code{XGenes}} {binary matrix with the simulated observations on the
#' level of genes, including noise. Rows represent samples and columns 
#' represent genes. }
#' \item{\code{XGenesClean}} {binary matrix with the simulated observations on
#' the level of genes, before noise was added. Rows represent samples and 
#' columns represent genes.}
#' \item{\code{XPathwaysClean}} {binary matrix with the simulated observations
#' on the level of pathways, before noise was added. Rows represent samples  
#' and columns represent pathways.}
#' \item{\code{TGenes}} {matrix with the waiting times of the genes. Rows 
#' represent samples and columns represent genes.}
#' \item{\code{TPathways}} {matrix with the waiting times of the pathways. The 
#' waiting time of each pathway equals the minimum waiting time of its gene 
#' members. Rows represent samples and columns represent pathways.}
#' \item{\code{Tobs}} {vector with the observation times. For each sample,
#' an independent observation time is drawn exponentially.}
#' \item{\code{ZGenes}} {matrix with the independent exponential components 
#' contributing to the waiting times of the genes. Rows represent samples 
#' and columns represent genes.}
#' \item{\code{m}} {number of samples}
#' \item{\code{model}} {structure with details about the model, s.a. waiting 
#' time rates, poset etc. This type of structure is returned by the functions
#' whcih generate pathTiMEx models: \code{generateRandomModel}, 
#' \code{generateFixedModel} and \code{generateStructureModel}.}
#' }
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases drawSamples
#' 
#' @export
drawSamples<-function(m,model)
{    
    poset <- model$poset
    lams <- model$lams
    lamobs <- model$lamobs
    assignment <- model$assignment
    mus <- model$mus
    epsilon <- model$epsilon
    
    n <- length(lams) # number of genes
    k <- ncol(poset) # number of pathways
    
    # Z is the independent exponential variable component specific  to each gene
    ZGenes <- matrix(NA, m, n)
    # T is Z plus the waiting times of the parents
    TGenes <- matrix(NA, m, n)
    
    # T is the actual value of the waiting time of the pathway, as the minimum 
    # of the waiting times of its gene members
    TPathways <- matrix(NA, m, k)
    
    # sample observation times
    Tobs <- rexp(m,rate=lamobs)
    
    
    # generate Z of the genes (exponentially independent)
    ZGenes<-t(replicate(m,rexp(n,rate=lams)))
        
    # go over the poset of pathways iteratively, by its topological sort
    topoPath <- topologicalSort(poset)
    
    for (e in topoPath)
    {
        members <- which(assignment == e) # members of pathway e
        parents <- which(poset[,e] == 1) # parents of pathway e
        
        # if the pathway has no parents, then the waiting times of its gene members
        # equal their exponential components
        if (length(parents) == 0)
        {
            TGenes[,members] <- ZGenes[,members]
            if (length(members)>1)
                TPathways[,e] <- apply(TGenes[,members],1,min)
            else
                TPathways[,e] <- TGenes[,members]
        }

        # if the pathway has a single parent, then the waiting times of its gene members
        # equal their exponential component plus the waiting time of its pathway parent
        else if (length(parents) == 1)
        {
            TGenes[,members] <- ZGenes[,members] + TPathways[,parents]
            if (length(members)>1)
                TPathways[,e] <- apply(TGenes[,members],1,min)
            else 
                TPathways[,e] <- TGenes[,members]
        }
        
        # if the pathway has more than one parent, then the waiting times of its gene members
        # equal their exponential component plus the maximum waiting time of its pathway parents
        else
        {
            TGenes[,members] <- ZGenes[,members] + apply(TPathways[,parents],1,max)
            if (length(members)>1)
                TPathways[,e] <- apply(TGenes[,members],1,min)
            else 
                TPathways[,e] <- TGenes[,members]
        }
    }
    
    XGenes <- matrix(NA,nrow=m,ncol=n)
    
    # binary matrices of observation X for genes and pathways
    XPathways <- (TPathways<Tobs)+0
    
    for (e in topoPath)
    {
        # X for genes
        members <- which(assignment == e) # members of pathway e
        # simulate a mutually exclusive pathway without noise
        XGenes[,members] <- sampleMEPathways(m,TGenes[,members],Tobs,mus[e],0)
    }
    
    colnames(XGenes) <- paste("gene",c(1:n),sep="")
    colnames(XPathways) <-  paste("pathway",c(1:k),sep="")
    
    # the binary matrix before adding noise
    XGenesClean<-XGenes
    
    # add noise to the genes with probability epsilon
    for (i in 1:m) # for each patient
    {
        # the genes to which noise is added, sampled binomially with success 
        # probability epsilon (1 turns into 0 and 0 turns into 1)
        genesNoise<-which(rbinom(prob=epsilon,n=dim(XGenes)[2],size=1)==1)
        if (length(genesNoise)>0) # if there is at least one gene for which noise is added
        {
            XGenes[i,genesNoise]<-(!XGenes[i,genesNoise])+0
        }
    }
    
    # the pathway matrix is returned before adding noise
    XPathwaysClean<-XPathways
    
    l <- list("XGenes" = XGenes, "XGenesClean" = XGenesClean, 
              "XPathwaysClean" = XPathwaysClean, "TGenes" = TGenes, 
              "TPathways" = TPathways, "Tobs" = Tobs, "ZGenes" = ZGenes, 
              "m" = m, "model" = model)
    
    return(l)   
}



###############################################################################
#' @title Draws samples from a mutually exclusive pathway
#' 
#' @description \code{sampleMEPathways} simulates a group of mutually exclusive
#' genes, for given waiting times of genes, observation time, intensity of
#' mutual exclusivity and noise level.
#' 
#' @param m number of samples to simulate, as a positive integer
#' @param TGenesTrunc matrix with waiting times of the genes; columns represent
#' the genes, rows represent the samples
#' @param Tobs vector of observation times, with length the number of samples
#' @param muTrunc intensity of mutual exclusivity of the pathway to be 
#' simulated, real number between 0 and 1
#' @param eps error rate at the level of genes, real number between 0 and 1
#' 
#' @details In the case of perfect mutual exclusivity (\code{mu=1}), for each
#' sample, only the gene with the waiting time can be mutated, and only if its
#' waiting time is lower than the observation time.
#' 
#' @return a binary matrix representing the simulated mutually exclusive 
#' group, with columns as genes and rows as samples
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases sampleMEPathways
#' 
#' @export
sampleMEPathways <- function(m,TGenesTrunc,Tobs,muTrunc,eps)
{
    # the number of genes in the pathway
    nTrunc <- dim(TGenesTrunc)[2]
    if (is.null(nTrunc))
    {
        nTrunc <- 1
        TGenesTrunc <- matrix(TGenesTrunc)
    }
    
    # the matrix of observation corresponding to the pathway
    xTrunc <- matrix(0,nrow = m,ncol = nTrunc)
    
    for (i in 1:m)
    {
        smalleq <- which(TGenesTrunc[i,] <= Tobs[i])
        minwt <- smalleq[which.min(TGenesTrunc[i,smalleq])]
        
        # only necessary if not all waiting times are larger than the observation time
        if (length(minwt)>0)
        {
            if (runif(1)<muTrunc) {
                # with probability mu, only X_min = 1, and all the rest are 0
                xTrunc[i,minwt] <- 1
                xTrunc[i,setdiff(c(1:nTrunc),minwt)] <- 0
            } else {
                # with probability 1-mu, all X corresponding to times < obs are 1 and the rest are 0
                xTrunc[i,smalleq] <- 1
                xTrunc[i, setdiff(c(1:nTrunc), smalleq)] <- 0
            }  
        }
    }
    
    return(xTrunc)
}



###############################################################################
#' @title Performs TiMEx iteratively on a simulated dataset
#' 
#' @description \code{iterTiMExSims} performs TiMEx iteratively on a simulated 
#' dataset, by each time removing the highest ranking group of largest size and 
#' runnning TiMEx again.
#' 
#' @param samplesModel structure as returned by the function \code{drawSamples},
#' containing the simulated pathTiMEx model, including the simulated binary 
#' matrices, the waiting times etc
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
#' @details The groups reported are always the most significant ones of largest
#' size. Once a group is reported, it is removed from the dataset and TiMEx is 
#' ran again. The groups are reported decreasingly by size. This function is
#' designed to be applied on data simulated by the functions in this package.
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
#' \item{\code{geneAssignments}} {list consisting of the initial assignment to
#' pathways of the genes in the identified groups. Each element is a vector of 
#' pathway indices corresponding to an identified mutually exclusive group. 
#' The names of each vector represent the indices of the genes in the binary 
#' input matrix.}
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
#' @aliases iterTiMExSims
#' 
#' @export
iterTiMExSims <- function(samplesModel,pairMu,pairPvalue,groupPvalue,corr,noRuns)
{
    mat<-samplesModel$XGenes
    
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
    
    realAssignment<-samplesModel$model$assignment
    
    # record the initial matrix and gene names for reconstructing the indices afterwards
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
            } else
            {
                geneNames[[i]]<-unlist(resTiMEx$genesSignif[[maxSize]][[idxName]])
                muEst[i]<-resTiMEx$MusGroup[[maxSize]][[idxName]]
                toRemove<-unlist(resTiMEx$idxSignif[[maxSize]][[idxName]])
            }
            mat<-mat[,-toRemove] # remove the identified group from the data
        }
        i<-i+1
    }
    
    # arrange the groups by size decreasingly
    muEst<-muEst[order(sapply(geneNames,length),decreasing = TRUE)]
    geneNames<-geneNames[order(sapply(geneNames,length),decreasing = TRUE)]
    
    # retrieve the initial indices
    geneIndices<-lapply(geneNames,function(x){
        match(x,initialNames)
    })
    
    # retrieve the assignments to groups
    i<-1
    geneAssignments <- list()
    while (i <= length(geneIndices))
    {
        vecIdxs <- geneIndices[[i]]
        geneAssignments[[i]] <- realAssignment[vecIdxs]
        i <- i+1
    }
    
    # create return structure
    l<-list("geneNames" = geneNames, "geneIndices" = geneIndices, 
            "geneAssignments" = geneAssignments, "samplesModel" = samplesModel,
            "muEst" = muEst,"corr" = corr)
    
    return(l)
}



###############################################################################
#' @title Assesses the performance of the inferred groups after running 
#' TiMEx iteratively
#' 
#' @description \code{clustMetrIterTiMEx} computes per-edge and per-pathway 
#' metrics evaluating the performance of the pathway assignment inferred after
#' running TiMEx iteratively, based on data simulated from PathTiMEx
#' 
#' @param resIterTiMEx structure containing the results from iteratively
#' running TiMEx on a binary dataset, as returned by the function
#' \code{iterTiMExSims}.
#' 
#' @details This function assesses how well the real groups were recovered 
#' after performing iterative inference with TiMEx.
#' 
#' @return list consisting of 
#' \itemize{
#' \item{\code{reconstruction}} {structure with metrics of the reconstuction of
#' the initial groups, as follows
#' \itemize{
#' \item{\code{groupsSpread}}{list containing as many elements as real 
#' pathways. each element is a list with at most as many elements as the 
#' number of inferred pathways (including potential NULL elements). 
#' each element of this second list  is a vector of genes part of the current
#' real pathway and also part of the inferred pathway. In othder words, the 
#' first index represents the real group and the second index represents the 
#' inferred group.}
#' \item{\code{groupsSpreadPercent}} {list containing as many elements as real 
#' pathways. each element is a list with at most as many elements as the 
#' number of inferred pathways (including potential NULL elements). 
#' each element of this second list  is the percentage of 
#' reconstruction of the respective first pathway by the genes identified as 
#' part of the second pathway. In othder words, the first index represents 
#' the real group and the second index represents the inferred group.}
#' \item{\code{groupsSpreadOrd}} {list containing as many elements as real 
#' pathways, ordered by the real mu values of the pathways. each element is a 
#' list with at most as many elements as the number of inferred pathways 
#' (including potential NULL elements). each element of this second list is a 
#' vector of genes part of that particular real pathway 
#' and also inferred as part of the current pathway. In othder words, the 
#' first index represents the real group and the second index represents the 
#' inferred group.}
#' \item{\code{groupsSpreadPercentOrd}} {list containing as many elements as 
#' real pathways, orderde by the real mu values of the pathways. 
#' each element is a list with at most as many elements as the 
#' number of inferred pathways (including potential NULL elements). each 
#' element of this second list  is the percentage of reconstruction of the 
#' respective first pathway by the genes identified as part of the second 
#' pathway. In othder words, the first index represents the real group and the 
#' second index represents the inferred group.}
#' \item{\code{groupsCompact}} {list containing as many elements as real 
#' pathways. each element is a vector consisting of the genes members of
#' that particular pathway which were identified as part of a mutually 
#' exclusive group, regardless of the group in which they were identified.}
#' \item{\code{groupsCompactPercent}} {list containing as many elements as real 
#' pathways. each element is the percentage of reconstruction of that 
#' particular pathway, i.e. the percentage of genes which were identified as 
#' mutually exclusive, regardless of the group in which they were identified.}
#' \item{\code{groupsCompactOrd}} {list containing as many elements as real 
#' pathways, ordered by the real mu values of the pathways. each element is a 
#' vector consisting of the genes members of that particular pathway which 
#' were identified as part of a mutually exclusive group, regardless of the 
#' group in which they were identified.}
#' \item{\code{groupsCompactPercentOrd}} {list containing as many elements as 
#' real pathways, ordered by the real mu values of the pathwayss. each 
#' element is the percentage of reconstruction of that particular pathway, 
#' i.e. the percentage of genes which were identified as mutually exclusive, 
#' regardless of the group in which they were identified.}
#' \item{\code{groupsEmpty}} {list containing as many elements as 
#' real pathways. each element contains the genes from that particular pathway
#' which were not identified as part of mutually exclusive groups.}
#' \item{\code{groupsEmptyPercent}} {list containing as many elements as 
#' real pathways. each element is the percentage of the respective pathway 
#' which was not identified as part of mutually exclusive groups.}
#' \item{\code{groupsEmptyOrd}} {list containing as many elements as 
#' real pathways, ordered by the real mu values of the pathways. each element 
#' contains the genes from that particular pathway which were not identified 
#' as mutually exclusive.}
#' \item{\code{groupsEmptyPercentOrd}} {list containing as many elements as 
#' real pathways, ordered by the real mu values of the pathways. each element 
#' is the percentage of the respective pathway which was not identified as 
#' part of mutually exclusive groups.}
#' \item{\code{booleanOthers}} {list containing as many elements as 
#' real pathways. each element is a list with at most as many elements as the 
#' number of inferred pathways (including potential NULL elements). 
#' each element of this second list is either 0, if the 
#' inferred group contained elements from more than one real pathway, or 1, if
#' the inferred group only contained elements from one real pathway. In othder 
#' words, the first index represents the real group and the second index 
#' represents the inferred group.}
#' \item{\code{booleanOthersOrd}}{list containing as many elements as 
#' real pathways, orderde by the real mu values of the pathways. each element 
#' is a list with at most as many elements as the number of inferred pathways 
#' (including potential NULL elements). each element of this second list is 
#' either 0, if the inferred group contained elements from more than one real 
#' pathway, or 1, if the inferred group only contained elements from one real 
#' pathway. In othder words, the first index represents the real group and the 
#' second index represents the inferred group.}
#' }}
#' \item{\code{characteristics}} {structure with metrics per pathway 
#' characterizing the identified groups, as follows
#' \itemize{
#' \item{\code{freqGenes}} {list containing as many elements as real pathways.
#' each element is a vector with frequencies of the genes inside that pathway.}
#' \item{\code{freqGenesAvg}} {list containing as many elements as real 
#' pathways. each element is the average frequency of all genes inside that 
#' pathway (averaged over all the genes in each pathway).}
#' \item{\code{freqGenesOrd}} {list containing as many elements as real 
#' pathways. each element is a vector with the frequencies of the 
#' genes inside that pathway, ordered decreasingly by the intensities of 
#' mutual exclusivity of each pathway.}
#' \item{\code{realMus}} {vector containing the intensities of mutual 
#' exclusivity of the real pathways}
#' \item{\code{realMusOrd}} {vector containing the intensities of mutual 
#' exclusivity of the real pathways, ordered decreasingly}
#' }}
#' \item{\code{peredge}} {structure with metrics per edge characterizing the 
#' reconstruction of groups, as follows
#' \itemize{
#' \item{\code{tp}} {true positive rate: the number of true mutual exclusivity
#' pairwise connections which were also identified, divided by the total 
#' number of true mutual exclusivity pairwise connections.}
#' \item{\code{fp}} {false positive rate: the number of false mutual 
#' exclusivity edges identified, divided by the total number of possible
#' pairwise connections.}
#' }}
#' \item{\code{resIterTiMEx}} {input structure containing the identified 
#' mutually exclusive groups after running TiMEx itartively, together with 
#' the simulated data on the basis of which the inference was done.}
#' }
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases clustMetrIterTiMEx
#' 
#' @export
clustMetrIterTiMEx <- function(resIterTiMEx)
{
    # number of genes
    n <- length(resIterTiMEx$samplesModel$model$lams)
    # real assignment used for simulating the data
    realAssignment <- resIterTiMEx$samplesModel$model$assignment
    # indices of the inferred groups
    detectedGroups <- resIterTiMEx$geneIndices 
    # mus of the real groups
    realMus<-resIterTiMEx$samplesModel$model$mus
    # binary matrix used as input (including noise)
    Xmatrix<-resIterTiMEx$samplesModel$XGenes 
    # the real poset
    poset<-resIterTiMEx$samplesModel$model$poset
    
    R <- as.relation(poset)
    RT <- transitive_closure(R)    
    res <- relation_incidence(RT)
    # number of parents of all upstream parents of each pathway in the poset
    noParents <- apply(res,2,sum)
    
    ## metrics per edge
    
    # the real mutual exclusivity connections, as edges in a square binary 
    # matrix
    edgesReal<-matrix(0,nrow=n,ncol=n)
    for (i in 1:n)
    {
        currentAssignment<-realAssignment[i]
        coMembers<-which(realAssignment==currentAssignment)
        edgesReal[i,coMembers]<-1
    }
    # the entries in the diagonal are spurious, so remove them
    diag(edgesReal)<-0
    # desymmetrize the matrix
    edgesReal[lower.tri(edgesReal)]<-0
    
    # count how many times each edge between any two genes was detected as 
    # part of any inferred mutually exclusive group
    edgesFound<-matrix(0, nrow = n, ncol = n)
    
    # for each detected group
    for (x in detectedGroups) 
    {
        # if there was at least one group of the current size detected
        if (!is.null(x)) 
        {
            # record each edge from each group
            if (length(x)>1)
            {
                edgesToAdd<-combinations(n=length(x),r=2,v=x)
                edgesFound[cbind(edgesToAdd[,1],edgesToAdd[,2])]<-edgesFound[cbind(edgesToAdd[,1],edgesToAdd[,2])] + 1
            }
            
        }
    }
    
    # compute the true positive and the false negative rates per identified edges
    
    # the number of true mutual exclusivity edges which were also identified
    tp<-length(intersect(which(edgesReal==1),which(edgesFound==1)))/length(which(edgesReal==1))
    # the number of falsely identified mutually exclusive edges
    fp<-length(setdiff(which(edgesFound==1),which(edgesReal==1)))/length(which(edgesReal==1))
    
    peredge<-list()
    peredge$tp<-tp
    peredge$fp<-fp
    
    ## metrics per pathway
    
    # the first index represents the real group, and the second index is 
    # the inferred group. the elements are the genes
    groupsSpread<-list() # as part of which groups each pathway was inferred
    
    freqsGenes<-list() # average frequencies of genes inside each real pathway
    
    for (i in c(1:length(unique(realAssignment)))) # for each real pathway
    {
        # all the members part of the current real pathway
        currentMembers<-which(realAssignment==i) 
        
        # average frequency of genes inside each pathway
        if (length(currentMembers)>1)
            freqsGenes[[i]]<-apply(Xmatrix[,currentMembers],2,mean)
        else
            freqsGenes[[i]]<-mean(Xmatrix[,currentMembers])
        
        # the reconstruction of each pathway
        groupsSpread[[i]]<-list()
        if (length(detectedGroups)>0)
        {
            for (x in c(1:length(detectedGroups))) # for each detected group
            {
                # genes part of pathway i which were found in group with index x
                newEls<-intersect(detectedGroups[[x]],currentMembers) 
                if (length(newEls)>0)
                    groupsSpread[[i]][[x]]<-newEls
            }  
        }
    }
    
    # average the frequencies of all genes in each group
    freqsGenesAvg<-lapply(freqsGenes,function(x){mean(x)})
    
    # the elements identified for each group, regardless of the group in 
    # which they were identified
    groupsCompact<-lapply(groupsSpread,function(x){unlist(x)})
    
    # the elements of the groups which were not at all idnentified
    groupsEmpty<-lapply(c(1:length(groupsCompact)),function(x){
        return(setdiff(which(realAssignment==x),groupsCompact[[x]]))
    })
    
    # whether any of the groups in groupsSpread also contained other elements
    booleanOthers<-lapply(groupsSpread,function(x){
        if (length(x)>0)
        {
            lapply(c(1:length(x)),function(y){
                if (!is.null(x[[y]]))
                {
                    # all members in the identified group
                    allMemsInGroup<-resIterTiMEx$geneIndices[[y]]
                    if (length(setdiff(allMemsInGroup,x[[y]])>0))
                        return (1)
                    else
                        return (0)                    
                }
            })  
        }
    })
    
    # groupsSpread as percentage
    groupsSpreadPercent<-lapply(c(1:length(groupsSpread)),function(x){
        sizeOfGroup<-length(which(realAssignment==x))
        lapply(groupsSpread[[x]],function(y){
            return(length(y)/sizeOfGroup)
        })  
    })
    
    # groupsCompact as percentage
    groupsCompactPercent<-lapply(c(1:length(groupsCompact)),function(x){
        sizeOfGroup<-length(which(realAssignment==x))
        return(length(groupsCompact[[x]])/sizeOfGroup)
    })
    
    # groupsEmpty as percentage
    groupsEmptyPercent<-lapply(c(1:length(groupsEmpty)),function(x){
        return(1-groupsCompactPercent[[x]])
    })
    
    orderedMus<-order(realMus,decreasing = TRUE)
    realMusOrd<-realMus[orderedMus]
    
    # order the return structures decreasingly by the real mus of the real pathways    
    groupsSpreadOrd<-groupsSpread[orderedMus]
    groupsSpreadPercentOrd<-groupsSpreadPercent[orderedMus]
    groupsCompactOrd<-groupsCompact[orderedMus]
    groupsCompactPercentOrd<-groupsCompactPercent[orderedMus]
    groupsEmptyOrd<-groupsEmpty[orderedMus]
    groupsEmptyPercentOrd<-groupsEmptyPercent[orderedMus]
    booleanOthersOrd<-booleanOthers[orderedMus]
    freqsGenesOrd<-freqsGenes[orderedMus]
    freqsGenesAvgOrd<-freqsGenesAvg[orderedMus]
    noParentsOrd<-noParents[orderedMus]
    
    # use lists as return structures
    reconstruction<-list()
    reconstruction$groupsSpread<-groupsSpread
    reconstruction$groupsSpreadPercent<-groupsSpreadPercent
    reconstruction$groupsSpreadOrd<-groupsSpreadOrd
    reconstruction$groupsSpreadPercentOrd<-groupsSpreadPercentOrd 
    reconstruction$groupsCompact<-groupsCompact
    reconstruction$groupsCompactPercent<-groupsCompactPercent
    reconstruction$groupsCompactOrd<-groupsCompactOrd
    reconstruction$groupsCompactPercentOrd<-groupsCompactPercentOrd
    reconstruction$groupsEmpty<-groupsEmpty
    reconstruction$groupsEmptyPercent<-groupsEmptyPercent
    reconstruction$groupsEmptyOrd<-groupsEmptyOrd
    reconstruction$groupsEmptyPercentOrd<-groupsEmptyPercentOrd
    reconstruction$booleanOthers<-booleanOthers
    reconstruction$booleanOthersOrd<-booleanOthersOrd
    
    characteristics<-list()
    characteristics$freqsGenes<-freqsGenes
    characteristics$freqsGenesAvg<-freqsGenesAvg
    characteristics$freqsGenesOrd<-freqsGenesOrd
    characteristics$realMus<-realMus
    characteristics$realMusOrd<-realMusOrd
    
    l<-list("reconstruction"=reconstruction,"characteristics"=characteristics,"peredge"=peredge,
            "resIterTiMEx"=resIterTiMEx)
    return(l)
}



###############################################################################
#' @title Evaluates performance of inference of groups
#' 
#' @description \code{returnMetricsClust} returns metrics evaluating the 
#' performance in group recovery of running TiMEx iteratively
#' 
#' @param metricsIterTiMEx structure with metrics of inference of the mutually
#' exclusive groups after running TiMEx iteratively, as returned by the 
#' function \code{clustMetrIterTiMEx}
#' 
#' @details 
#' 
#' @return list consisting of 
#' \itemize{
#' \item{\code{percent}} {vector with length the number of real mutually 
#' exclusive pathways. each element represents the maximal percentage of the 
#' current group which was reconstructed by another group (real number between
#' 0 and 1).} 
#' \item{\code{alone}} {vector with length the number of real mutually 
#' exclusive pathways. each element is either 0, if the initial group has not
#' been recovered by multiple groups, NA if it was and 2 if it hasn't been 
#' identified at all.}
#' }
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases returnMetricsClust
#' 
#' @export
returnMetricsClust<-function(metricsIterTiMEx)
{
    
    percent<-sapply(metricsIterTiMEx$reconstruction$groupsSpreadPercent,function(x){
        if (length(x)>0)
        {
            v<-unlist(x)
            return(max(v))  
        }
        else return(0)
    })
    
    alone<-sapply(metricsIterTiMEx$reconstruction$groupsSpreadPercent,function(x){
        v<-unlist(x)
        if (length(v)>0)
        {
            if (v[which.max(v)]==1)
            {
                v1<-v[-(which.max(v))]
                if (length(v1)>0)
                    p<-max(v1)
                else p<-0
            }
            else return(NA)  
        }
        else return(2)
    })
    
    l<-list("percent" = percent, "alone" = alone)
    return(l)
}



###############################################################################
#' @title Appends the genes as single clusters
#' 
#' @description \code{addAllGenesNotInTiMEx} appends the genes not indentified
#' in any mutually exclusive groups as clusters of single genes.
#' 
#' @param resIterTiMEx structure containing the results from iteratively
#' running TiMEx on a binary dataset, as returned by the function
#' \code{iterTiMExSims}.
#' 
#' @details The function returns the updated clustering which also includes
#' the single genes.
#' 
#' @return groupysL list with the identified clusters. each element is a cluster
#' with gene indices.
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases addAllGenesNotInTiMEx
#' 
#' @export
addAllGenesNotInTiMEx<-function(resIterTiMEx)
{
    Datamat<-resIterTiMEx$samplesModel$XGenes
    groupysL<-resIterTiMEx$geneIndices
    # also add the genes which were not part of the groupings, as separate events
    groupysL<-append(groupysL,as.list(setdiff(c(1:dim(Datamat)[2]),unlist(groupysL))))
    
    return(groupysL)
}