pathTiMEx<-function(initialGroupsStruct,Datamat,path,name,noReps,optionsSA,
                    noThreads,numsave,skipsteps,gamma,noNoOpen,additionalGenes,
                    limitChanged)
{
    if (missing(noReps))
        noReps<-100
    
    if (missing(optionsSA))
        optionsSA<-"-s -T 1 -N 200"
    
    if (missing(noThreads))
        noThreads<-4
    
    if (missing(numsave))
        numsave<-2000
    
    if (missing(skipsteps))
        skipsteps<-100
    
    if (missing(gamma))
        gamma<-0.1
    
    if (missing(noNoOpen))
        noNoOpen<-100
    
    if (missing(additionalGenes))
        additionalGenes<-NULL
    
    # indicates the number of times the grouping and progression need to not change to finish the optimization
    if (missing(limitChanged))
        limitChanged<-2 
    
    groupysL<-list() # list to store all groupings at every step
    PathmatsL<-list() # list to store all matrices patients x pathways at every step
    unidetifEvsL<-list() # list to store whether the PathmatsL elements had only unique events at every step (via doMetagene)
    cbnResultL<-list() # list to store all cbnResult structures at every step
    optStructL<-list() # list to store all optimal structure at every step
    epsL<-list() # list to store the estimated values of epsilon for the ML poset at every step
    alphaL<-list() # list to store the values of alpha of the ML poset at every step
    loglikeL<-list() # list to store the optimal log likelihod of the ML poset at every step
    lamobsL<-list() # list to store the estimated values of lambda obs for the ML poset at every step
    lambdasL<-list() # list to store the estimated values of lambda for the ML poset at every step
    minScore<-c() # vector to store the minimum scores from the local optimization of groups at every step
    clustSimConsecJacc<-c() # vector to store the Jaccard similarity index between two consecutive clusterings
    clustSimConsecRand<-c() # vector to store the Rand similarity index between two consecutive clusterings
    
    # use the groups inferred by TiMEx as an initial solution
    groupysL[[1]]<-createInitialGroups(initialGroupsStruct,additionalGenes,Datamat)
    #groupysL[[1]]<-append(groupysL[[1]],as.list(setdiff(c(1:dim(Datamat)[2]),unlist(groupysL[[1]])))) # also add the genes which were not part of the groupings, as separate events
    
    PathmatsL[[1]]<-doPatForCBN(groupysL[[1]], Datamat, write=TRUE, path=path,name=name) # the initial binary alteration matrix of pathways
    unidetifEvsL[[1]]<-doMetagene(PathmatsL[[1]])$groups
    startPos<-make_linear_poset(dim(PathmatsL[[1]])[2]) # do a linear poset as initial solution for SA and the structure used in the first optimization run
    writePosetForCBN(poset=startPos,path=path,name=name) # writes the poset to a file, to be used as input by H-CBN
    cbnResultL[[1]]<-runCBN(path=path,name=name,optionsSA=optionsSA,noThreads=noThreads)
    optStructL[[1]]<-cbnResultL[[1]]$optStruct
    
    epsL[[1]]<-cbnResultL[[1]]$eps
    alphaL[[1]]<-cbnResultL[[1]]$alpha
    loglikeL[[1]]<-cbnResultL[[1]]$loglik
    lamobsL[[1]]<-cbnResultL[[1]]$lamobs
    lambdasL[[1]]<-cbnResultL[[1]]$lams
    
    i<-2
    cont<-TRUE # indicates whether to still continue the optimization
    notChanged<-0 # indicates whether the inferred grouping and progression haven't changed since the previous iteration
    timeElapsed<-c() # time elapsed for one iteration run of the optimizer
    timeElapsed[1]<-NA
    accuracyStructConsec<-rep(NA,noReps)
    accuracyStructConsec[1]<-NA
    
    while (i <= (noReps+1) && cont==TRUE)
    {
        print(paste("Repetition =", i, sep=" "))
        ptm <- proc.time()
        
        # keep the structure from repetition i-1 fixed and optimize for groups
        groupys<-groupysL[[i-1]]
        adjmat<-optStructL[[i-1]]
        suppressWarnings(rm(optimGroups))
        
        print("optimizing groupings")
        if (i<=noNoOpen)
            optimGroups<-optGroups2(Datamat = Datamat, groupys = groupys, 
                                    gamma = gamma, numsave = numsave, 
                                    skipsteps = skipsteps, adjmat = adjmat)
        else
            optimGroups<-optGroups(Datamat = Datamat, groupys = groupys, 
                                   gamma = gamma, numsave = numsave, 
                                   skipsteps = skipsteps, adjmat = adjmat)
        
        # store the new groups
        groupysL[[i]]<-cleanGroupys(optimGroups$mingroupys)
        minScore[i-1]<-optimGroups$mingroupysscore
        
        # store and write the new pathway matrix
        PathmatsL[[i]]<-doPatForCBN(groupysL[[i]], Datamat, write=TRUE, path=path,name=name)
        unidetifEvsL[[i]]<-doMetagene(PathmatsL[[i]])$groups
        
        # transform the existing optimal poset from step i-1 to a new compatible one
        transformAndWritePosetLinear(path,name,optStructL[[i-1]],PathmatsL[[i]])
        
        # keep the groups from repetition i fixed and optimize for structure, 
        # using as initial solution the transformed previously optimized structure
        cbnResultL[[i]]<-runCBN(path=path,name=name,optionsSA=optionsSA,noThreads=noThreads)
        optStructL[[i]]<-cbnResultL[[i]]$optStruct
        
        # try to merge groups with same parents and which leave the score unchanged
        mergedG<-mergeGroups(Datamat = Datamat, posetCurrent = optStructL[[i]], groupysCurrent = groupysL[[i]])
        if (mergedG$found==1) 
        {
            print(paste("merging of groups done"))
            groupysL[[i]]<-mergedG$groupysMerged
            optStructL[[i]]<-mergedG$posetMerged
        }
        
        epsL[[i]]<-cbnResultL[[i]]$eps
        alphaL[[i]]<-cbnResultL[[i]]$alpha
        loglikeL[[i]]<-cbnResultL[[i]]$loglik
        lamobsL[[i]]<-cbnResultL[[i]]$lamobs
        lambdasL[[i]]<-cbnResultL[[i]]$lams
        
        # check if the clustering has changed between this run and the previous one
        vecAssignPrev<-returnAssignment(groupysL[[i-1]])
        vecAssignNow<-returnAssignment(groupysL[[i]])
        newStruct<-optStructL[[i]]
        clustSimConsecJacc[i-1]<-cluster_similarity(vecAssignNow,vecAssignPrev,similarity = "jaccard")
        clustSimConsecRand[i-1]<-cluster_similarity(vecAssignNow,vecAssignPrev,similarity = "rand")
        print(paste("Consecutive Jaccard cluster similarity = ",clustSimConsecJacc[i-1]),sep="")
        print(paste("Consecutive Rand cluster similarity = ",clustSimConsecRand[i-1]),sep="")
        if (clustSimConsecRand[i-1]==1)
        {   
            # map structure i to structure i-1, for comparison
            oldStruct<-mapStructures(vecAssignPrev,vecAssignNow,optStructL[[i-1]])
            accuracyStructConsec[i]<-compareStructs(oldStruct,newStruct)$accuracy
            
            # if clustering didn't change, check if structure stayed the same by 
            # comparing the two structures
            if (compareStructs(oldStruct,newStruct)$accuracy==1)
                notChanged<-notChanged+1
            else
                notChanged<-0
            print(paste("progression accuracy = ",compareStructs(oldStruct,newStruct)$accuracy),sep="")
        } else
            notChanged<-0
        
        if (notChanged==limitChanged)
            cont<-FALSE
        
        print(paste("notChanged =", notChanged))
        timeElapsed[i]<-(proc.time() - ptm)[3]
        i<-i+1
    }
    
    result<-list()
    result$epsL<-epsL
    result$Datamat<-Datamat
    result$groupysL<-groupysL
    result$PathmatsL<-PathmatsL
    result$unidetifEvsL<-unidetifEvsL
    result$cbnResultL<-cbnResultL
    result$optStructL<-optStructL
    result$epsL<-epsL
    result$alphaL<-alphaL
    result$loglikeL<-loglikeL
    result$lamobsL<-lamobsL
    result$lambdasL<-lambdasL
    result$minScore<-minScore
    result$clustSimConsecJacc<-clustSimConsecJacc
    result$clustSimConsecRand<-clustSimConsecRand
    result$accuracyStructConsec<-accuracyStructConsec
    result$limitChanged<-limitChanged
    result$timeElapsed<-timeElapsed
    result$noIters<-i-1
    result$notChanged<-notChanged
    result$limitChanged<-limitChanged
    
    return(result)
}