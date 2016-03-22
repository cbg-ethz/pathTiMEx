###############################################################################
#' @title Locally optimizes pathway assignment
#' 
#' @description \code{optGroups} optimizes the assignment of genes to pathways,
#' for given fixed structure. 
#' 
#' @param Datamat binary alteration matrix where rows are samples and columns
#' are genes and an 1 on position [i,j] means that gene j is altered in sample
#' i.
#' @param groupys list with as many elements as mutually exclusive groups,
#' each element consists of the genes assigned to the respective pathway.
#' @param gamma parameter used in optimization to tune whether the MCMC 
#' sampler: a larger value encourages the chain to explore each local minimum 
#' more thoroughly, while a lower value allows the chain to escape each local 
#' area more easily and explore more of the global space. Usually set to 0.5.
#' @param numsave number of minimal scores attained by the chains, to bs used 
#' for plotting.
#' @param skipsteps number representing once in how many iterations to save
#' the minimal score for plotting
#' @param adjmat binary square matrix, representing the optimal structure which 
#' is kept fixed. The matrix contains one dummy column and row at the end, 
#' identically 0, hence its dimension is number of groups +1.
#' 
#' @details The optimization is done by via an MCMC chain, by minimizing the 
#' number of contradictions due to both mutual exclusivity and progression, 
#' as the number of ones that would need to be changed to zeroes to ensure 
#' consistency of the data with the clusters and structure. The number of
#' iterations in each chain equals \code{numsave*skipsteps}. For each run, 
#' opening a new mutually exclusive group is allowed. 
#' 
#' @return list consisting of 
#' \itemize{
#' \item{\code{mingroupys}} {optimized clustering solution, which minimizes the
#' number of contradictions due to both mutual exclusivity and progression: 
#' list with as many elements as mutually exclusive groups, each element 
#' consists of the genes assigned to the respective pathway. }
#' \item{\code{mingroupysscore}} {minimal score, namely the number of 
#' contradictions due to both mutual exclusivity and progression.}
#' \item{\code{accepted}} {number of times the proposed move was accepted.}
#' \item{\code{allscores}} {the scores attained at each \code{skipsteps} 
#' iterations.}
#' }
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases optGroups
#' 
#' @export
optGroups<-function(Datamat,groupys,gamma,numsave,skipsteps,adjmat)
{
    ng<-length(groupys)
    accepted<-0
    np<-ng+1
    groupys<-append(groupys,list(NULL))
    
    NN<-dim(Datamat)[1] # number of observations
    nn<-dim(Datamat)[2] # number of genes
    
    numits<-numsave*skipsteps # how many iterations in the chain
    allscores<-rep(0,numsave+1) # stores the scores in the chain
    
    posey<-rep(0,nn) # to store the group of each gene
    for(kk in 1:ng){
        posey[groupys[[kk]]]<-kk
    }
    
    
    ### Summarise the values of the nodes
    
    nodesums<-matrix(0,nrow=NN,ncol=np)
    for(kk in 1:np){ # add all the observations in each group
        nodesums[,kk]<-rowSums(as.matrix(Datamat[,groupys[[kk]]]))
    }
    nodevalues<-sign(nodesums) # value of the mutual exclusive pathways
    
    ### Within each group count anything with more than two 1s as errors
    
    nodeerrors<-nodesums-1 # count anything above 1 as errors
    nodeerrors[which(nodeerrors<0)]<-0 # give no errors to 0s
    
    withinnodeerrors<-colSums(nodeerrors) # store the total number of errors within each grouping
    
    ### Between groups, count 0 parents leading to a 1 as errors
    ### note that having move than a single 1 will already be counted as mutually exclusive errors
    
    betweennodeerrors<-rep(0,np)
    for(kk in 1:np){ # find parent values
        parentnodes<-which(adjmat[,kk]>0)
        lp<-length(parentnodes)
        if(lp>0){ # find where any parent nodes are 0
            anyparentsoff<-which(rowSums(as.matrix(nodevalues[,parentnodes]))<lp)
            betweennodeerrors[kk]<-sum(nodevalues[anyparentsoff,kk])
        }
    }
    
    groupysscore<-sum(withinnodeerrors+betweennodeerrors)
    
    allscores[1]<-groupysscore # store the first score
    
    ### store the minimum number of errors found
    mingroupysscore<-groupysscore
    mingroupys<-groupys
    
    for(iterations in 1:numits){
        genetomove<-sample.int(nn,1) # sample which gene to reassign
        nodetomovefrom<-posey[genetomove] # which node it was in
        nodetomoveto<-sample.int(np,1) # sample which node to place it in
        if(nodetomoveto!=nodetomovefrom){ # check we don't try to place it in the same node
            affectednodes<-c(posey[genetomove],nodetomoveto)
            # create the proposed grouping
            proposedgroupys<-groupys
            proposedgroupys[[nodetomoveto]]<-c(proposedgroupys[[nodetomoveto]],genetomove) # put it in its new home and take it away from the old one
            proposedgroupys[[nodetomovefrom]]<-proposedgroupys[[nodetomovefrom]][-which(proposedgroupys[[nodetomovefrom]]==genetomove)]
            
            # update the list of where each node is
            proposedposey<-posey
            proposedposey[genetomove]<-nodetomoveto
            
            # update the node sums and values
            proposednodesums<-nodesums
            proposednodesums[,affectednodes[1]]<-nodesums[,affectednodes[1]]-Datamat[,genetomove] # take away the 1s from the moved gene
            proposednodesums[,affectednodes[2]]<-nodesums[,affectednodes[2]]+Datamat[,genetomove] # and add them here
            
            propaffnodes<-proposednodesums[,affectednodes] # take the sums of the affected nodes
            proposednodevalues<-nodevalues
            proposednodevalues[,affectednodes]<-sign(propaffnodes) # value of the mutual exclusive pathways
            
            propaffnodeerrors<-propaffnodes-1 # count anything above 1 as errors
            propaffnodeerrors[which(propaffnodeerrors<0)]<-0 # give no errors to 0
            
            proposedwithinnodeerrors<-withinnodeerrors
            proposedwithinnodeerrors[affectednodes]<-colSums(propaffnodeerrors) # store the total number of errors within each grouping
            
            ### Between groups, count 0 parents leading to a 1 as errors
            ### note that we only need to check parents and children of the affected nodes
            
            proposedbetweennodeerrors<-betweennodeerrors
            # the nodes we need to rescore are those that gain or lose an element and their children  
            allaffectednodes<-union(union(which(adjmat[affectednodes[1],]>0),which(adjmat[affectednodes[2],]>0)),union(affectednodes[1],affectednodes[2]))
            for(kk in allaffectednodes){ # find parent values
                parentnodes<-which(adjmat[,kk]>0)
                lp<-length(parentnodes)
                if(lp>0){ # find where any parent nodes are 0
                    anyparentsoff<-which(rowSums(as.matrix(proposednodevalues[,parentnodes]))<lp)
                    proposedbetweennodeerrors[kk]<-sum(proposednodevalues[anyparentsoff,kk])
                }
            }
            
            proposedgroupysscore<-sum(proposedwithinnodeerrors+proposedbetweennodeerrors)
            
            ### Now do an MCMC type move
            if(runif(1)< exp(gamma*(groupysscore-proposedgroupysscore)) ){# accept the move
                accepted<-accepted+1
                groupys<-proposedgroupys
                posey<-proposedposey
                groupysscore<-proposedgroupysscore
                nodesums<-proposednodesums
                nodevalues<-proposednodevalues
                withinnodeerrors<-proposedwithinnodeerrors
                betweennodeerrors<-proposedbetweennodeerrors
                if(groupysscore<mingroupysscore){
                    mingroupysscore<-groupysscore
                    mingroupys<-groupys
                }
            }
            
            
        }
        if((iterations%%skipsteps)==0)
            allscores[iterations/skipsteps+1]<-groupysscore
    }
    
    return(list("mingroupys" = mingroupys, "mingroupysscore" = mingroupysscore, 
                "accepted" = accepted, "allscores" = allscores))
}



###############################################################################
#' @title Removes empty groups
#' 
#' @description \code{cleanGroupys} removes the empty elements in a list, 
#' representing the assignment of genes to pathways. 
#' 
#' @param groupysInput list with as many elements as mutually exclusive groups,
#' each element consists of the genes assigned to the respective pathway. It
#' possibily contains empty elements.
#' 
#' @details The empty elements are a consequence of the optimization step, in
#' which a grouping with less pathways than initially attains the minimal 
#' score.
#' 
#' @return list with as many elements as mutually exclusive groups,
#' each element consists of the genes assigned to the respective pathway. It
#' doesn't contain any empty elements.
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases cleanGroupys
#' 
#' @export
cleanGroupys<-function(groupysInput)
{
    return(groupysInput[which(lapply(groupysInput,length)>0)])
}



###############################################################################
#' @title Returns assignment to groups as vector
#' 
#' @description \code{returnAssignment} returns the assignment to pathways of 
#' a list of gene groups, as a vector of groups.
#' 
#' @param groupysInput list with as many elements as mutually exclusive groups,
#' each element consists of the genes assigned to the respective pathway
#' 
#' @details The vector representation of a clustering is used for comparison
#' between clusterings.
#' 
#' @return a vector where each position represents a gene and each entry 
#' represents the group to which it is assigned
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases returnAssignment
#' 
#' @export
returnAssignment<-function(groupysInput)
{
    assgn<-rep(0,max(unlist(groupysInput)))
    for (l in 1:length(groupysInput))
        assgn[groupysInput[[l]]]<-l
    return(assgn)
}



###############################################################################
#' @title Computes clustering score
#' 
#' @description \code{computeScore} computes the score of a given clustering, 
#' for a given structure.
#' 
#' @param Datamat binary alteration matrix where rows are samples and columns
#' are genes and an 1 on position [i,j] means that gene j is altered in sample
#' i.
#' @param groupys list with as many elements as mutually exclusive groups,
#' each element consists of the genes assigned to the respective pathway.
#' @param adjmat binary matrix, representing the optimal structure which is
#' kept fixed. The matrix contains one dummy column and row at the end, 
#' identically 0, hence its dimension is number of groups +1.
#' 
#' @details The score equals the number of contradictions due to both mutual 
#' exclusivity and progression, or the number of ones that would need to be 
#' changed to zeroes to ensure consistency of the data with the clusters and 
#' structure.
#' 
#' @return clustering score (integer)
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases computeScore
#' 
#' @export
computeScore<-function(Datamat,groupys,adjmat)
{
    NN<-dim(Datamat)[1] # number of observations
    nn<-dim(Datamat)[2] # number of genes
    ng<-length(groupys) # number of mutually exclusive groups
    np<-ng # no extra group opened
    
    posey<-rep(0,nn) # to store the group of each gene
    for(kk in 1:ng){
        posey[groupys[[kk]]]<-kk
    }
    
    
    ### Summarise the values of the nodes
    nodesums<-matrix(0,nrow=NN,ncol=np)
    for(kk in 1:np){ # add all the observations in each group
        nodesums[,kk]<-rowSums(as.matrix(Datamat[,groupys[[kk]]]))
    }
    nodevalues<-sign(nodesums) # value of the mutual exclusive pathways
    
    ### Within each group count anything with more than two 1s as errors
    nodeerrors<-nodesums-1 # count anything above 1 as errors
    nodeerrors[which(nodeerrors<0)]<-0 # give no errors to 0s
    withinnodeerrors<-colSums(nodeerrors) # store the total number of errors within each grouping
    
    ### Between groups, count 0 parents leading to a 1 as errors
    # note that having move than a single 1 will already be counted as mutually exclusive errors
    betweennodeerrors<-rep(0,np)
    for(kk in 1:np){ # find parent values
        parentnodes<-which(adjmat[,kk]>0)
        lp<-length(parentnodes)
        if(lp>0){ # find where any parent nodes are 0
            anyparentsoff<-which(rowSums(as.matrix(nodevalues[,parentnodes]))<lp)
            betweennodeerrors[kk]<-sum(nodevalues[anyparentsoff,kk])
        }
    }
    
    groupysscore<-sum(withinnodeerrors+betweennodeerrors)
    
    return(groupysscore)
}



optGroups2<-function(Datamat,groupys,gamma,numsave,skipsteps,adjmat)
{
    ng<-length(groupys)
    
    accepted<-0
    np<-ng
    
    dimadj<-dim(adjmat)[1]
    adjmat<-adjmat[-dimadj,]
    adjmat<-adjmat[,-dimadj]
    
    
    NN<-dim(Datamat)[1] # number of observations
    nn<-dim(Datamat)[2] # number of genes
    
    numits<-numsave*skipsteps # how many iterations in the chain
    allscores<-rep(0,numsave+1) # stores the scores in the chain
    
    posey<-rep(0,nn) # to store the group of each gene
    for(kk in 1:ng){
        posey[groupys[[kk]]]<-kk
    }
    
    
    ### Summarise the values of the nodes
    
    nodesums<-matrix(0,nrow=NN,ncol=np)
    for(kk in 1:np){ # add all the observations in each group
        nodesums[,kk]<-rowSums(as.matrix(Datamat[,groupys[[kk]]]))
    }
    nodevalues<-sign(nodesums) # value of the mutual exclusive pathways
    
    ### Within each group count anything with more than two 1s as errors
    
    nodeerrors<-nodesums-1 # count anything above 1 as errors
    nodeerrors[which(nodeerrors<0)]<-0 # give no errors to 0s
    
    withinnodeerrors<-colSums(nodeerrors) # store the total number of errors within each grouping
    
    ### Between groups, count 0 parents leading to a 1 as errors
    ### note that having move than a single 1 will already be counted as mutually exclusive errors
    
    betweennodeerrors<-rep(0,np)
    for(kk in 1:np){ # find parent values
        parentnodes<-which(adjmat[,kk]>0)
        lp<-length(parentnodes)
        if(lp>0){ # find where any parent nodes are 0
            anyparentsoff<-which(rowSums(as.matrix(nodevalues[,parentnodes]))<lp)
            betweennodeerrors[kk]<-sum(nodevalues[anyparentsoff,kk])
        }
    }
    
    groupysscore<-sum(withinnodeerrors+betweennodeerrors)
    
    allscores[1]<-groupysscore # store the first score
    
    ### store the minimum number of errors found
    mingroupysscore<-groupysscore
    mingroupys<-groupys
    
    for(iterations in 1:numits){
        genetomove<-sample.int(nn,1) # sample which gene to reassign
        nodetomovefrom<-posey[genetomove] # which node it was in
        nodetomoveto<-sample.int(np,1) # sample which node to place it in
        if(nodetomoveto!=nodetomovefrom){ # check we don't try to place it in the same node
            affectednodes<-c(posey[genetomove],nodetomoveto)
            # create the proposed grouping
            proposedgroupys<-groupys
            proposedgroupys[[nodetomoveto]]<-c(proposedgroupys[[nodetomoveto]],genetomove) # put it in its new home and take it away from the old one
            proposedgroupys[[nodetomovefrom]]<-proposedgroupys[[nodetomovefrom]][-which(proposedgroupys[[nodetomovefrom]]==genetomove)]
            
            # update the list of where each node is
            proposedposey<-posey
            proposedposey[genetomove]<-nodetomoveto
            
            # update the node sums and values
            proposednodesums<-nodesums
            proposednodesums[,affectednodes[1]]<-nodesums[,affectednodes[1]]-Datamat[,genetomove] # take away the 1s from the moved gene
            proposednodesums[,affectednodes[2]]<-nodesums[,affectednodes[2]]+Datamat[,genetomove] # and add them here
            
            propaffnodes<-proposednodesums[,affectednodes] # take the sums of the affected nodes
            proposednodevalues<-nodevalues
            proposednodevalues[,affectednodes]<-sign(propaffnodes) # value of the mutual exclusive pathways
            
            propaffnodeerrors<-propaffnodes-1 # count anything above 1 as errors
            propaffnodeerrors[which(propaffnodeerrors<0)]<-0 # give no errors to 0
            
            proposedwithinnodeerrors<-withinnodeerrors
            proposedwithinnodeerrors[affectednodes]<-colSums(propaffnodeerrors) # store the total number of errors within each grouping
            
            ### Between groups, count 0 parents leading to a 1 as errors
            ### note that we only need to check parents and children of the affected nodes
            
            proposedbetweennodeerrors<-betweennodeerrors
            # the nodes we need to rescore are those that gain or lose an element and their children  
            allaffectednodes<-union(union(which(adjmat[affectednodes[1],]>0),which(adjmat[affectednodes[2],]>0)),union(affectednodes[1],affectednodes[2]))
            for(kk in allaffectednodes){ # find parent values
                parentnodes<-which(adjmat[,kk]>0)
                lp<-length(parentnodes)
                if(lp>0){ # find where any parent nodes are 0
                    anyparentsoff<-which(rowSums(as.matrix(proposednodevalues[,parentnodes]))<lp)
                    proposedbetweennodeerrors[kk]<-sum(proposednodevalues[anyparentsoff,kk])
                }
            }
            
            proposedgroupysscore<-sum(proposedwithinnodeerrors+proposedbetweennodeerrors)
            
            ### Now do an MCMC type move
            if(runif(1)< exp(gamma*(groupysscore-proposedgroupysscore)) ){# accept the move
                accepted<-accepted+1
                groupys<-proposedgroupys
                posey<-proposedposey
                groupysscore<-proposedgroupysscore
                nodesums<-proposednodesums
                nodevalues<-proposednodevalues
                withinnodeerrors<-proposedwithinnodeerrors
                betweennodeerrors<-proposedbetweennodeerrors
                if(groupysscore<mingroupysscore){
                    mingroupysscore<-groupysscore
                    mingroupys<-groupys
                }
            }
            
            
        }
        if((iterations%%skipsteps)==0)
            allscores[iterations/skipsteps+1]<-groupysscore
    }
    
    #print(mingroupysscore)
    return(list("mingroupys" = mingroupys, "mingroupysscore" = mingroupysscore, 
                "accepted" = accepted, "allscore" = allscores))
}


# merge groups with the same parents
mergeGroups<-function(Datamat,posetCurrent,groupysCurrent)
{
    # the final merged groups to be returned
    groupysMerged<-groupysCurrent
    
    # boolean to indicate whether another solution has been found
    found<-0
    
    # remove the last row and column of the poset (identical 0)
    dimpos<-dim(posetCurrent)[1]
    posetCurrent<-posetCurrent[-dimpos,]
    posetCurrent<-posetCurrent[,-dimpos]
    
    # the final merged poset to be returned
    posetMerged<-posetCurrent
    
    toChange<-c()
    
    # search duplicated columns, which correspond to genes having the exact same parents
    dups<-which(duplicated(posetCurrent, MARGIN = 2))
    if (length(dups)>0) # in case there exist duplicates
    {
        for (d in dups) # for each column which has at least one duplicate
        {
            # all the events identical to d
            equalNow<-which(apply(posetCurrent,2,function(x){identical(x,posetCurrent[,d])})==1)
            
            # all pairs of duplicated events
            pairs<-t(combn(x=equalNow,m=2))
            
            # the current score, before the tentative merge
            scoreCurrent<-computeScore(Datamat,groupysCurrent,posetCurrent)
            
            for (p in 1:dim(pairs)[1]) # for each pair, which can tentatively be merged
            {
                # the pair to work with now
                pairNow<-pairs[p,]
                
                # the new tentative grouping after the merge
                groupysTentative<-groupysCurrent
                groupysTentative[[pairNow[1]]]<-unlist(groupysTentative[pairNow])
                groupysTentative[[pairNow[2]]]<-NULL
                groupysTentative<-cleanGroupys(groupysTentative)
                
                # the new tentative poset
                posetTentative<-posetCurrent
                posetTentative<-posetTentative[-pairNow[2],]
                posetTentative<-posetTentative[,-pairNow[2]]
                posetTentative<-cbind(posetTentative,0)
                posetTentative<-rbind(posetTentative,0)
                
                # the new tentative score
                scoreTentative<-computeScore(Datamat,groupysTentative,posetTentative)
                
                # change if the new score is the same as the current one
                if (scoreTentative==scoreCurrent)
                {
                    toChange<-rbind(toChange,pairNow)
                    rownames(toChange)<-NULL
                }
            }
        }
    }
    
    # update the groups with the pairs
    if (length(toChange)>0)
    {
        changedVec<-rep(0,length(groupysCurrent))
        found<-1
        dups2<-which(duplicated(toChange, MARGIN = 1))
        toChange<-toChange[-dups2,]
        if (is.null(dim(toChange)))
            toChange<-t(as.matrix(toChange))
        unq<-unique(as.vector(toChange))
        for (kk in unq)
        {
            if (changedVec[kk]==0)
            {
                sameGroup<-unique(as.vector(toChange[which(toChange==kk,arr.ind = TRUE)[,1],]))
                changedVec[sameGroup]=1
                groupysMerged[[sameGroup[1]]]<-unlist(groupysCurrent[sameGroup])
                groupysMerged[setdiff(sameGroup,sameGroup[1])]<-NULL
                
                posetMerged<-posetMerged[-setdiff(sameGroup,sameGroup[1]),]
                posetMerged<-posetMerged[,-setdiff(sameGroup,sameGroup[1])]
            }
        }
    }
    
    posetMerged<-cbind(posetMerged,0)
    posetMerged<-rbind(posetMerged,0)
    
    return(list("groupysMerged" = groupysMerged, "posetMerged" = posetMerged,
                "found" = found))
}

