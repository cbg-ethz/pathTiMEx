# Author: Simona Constantinescu; simona.constantinescu@bsse.ethz.ch

###############################################################################
#' @title Generates a poset
#' 
#' @description \code{genRandomPoset} generates a random poset (partial
#' order at the level of pathways) with a given number of events (pathways) 
#' and a given edge density.
#' 
#' @param k number of events (pathways)
#' @param ldense ratio between 0 and 1 indicating how dense (in terms of 
#' edges) the poset is
#' @param trans_reduced boolean indicating whether the poset should be a
#' transitive reduction (default to \code{TRUE})
#' 
#' @details The poset is a binary upper triangular matrix in which if 
#' element [i,j]=1, then there is an order relation between i and j, i.e. event
#' i happens before or at the same time as event j. The 
#' parameter \code{ldenses} represents the fraction (between 0 and 1) of 
#' edges from the maximal number of edges (which contains nr_muts choose 2 
#' edges) to be included in the poset. 
#' 
#' @return a poset (a binary upper triangular matrix)
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases genRandomPoset
#' 
#' @export
genRandomPoset <- function(k , ldenses, trans_reduced = TRUE)
{    
    # number of edges for the poset
    nr_edges <- ceiling(choose(k, 2)*ldenses)
    
    # randomly add the determined number of edges in the upper triangular 
    # part of the poset
    poset <- matrix(0, k, k)
    poset[upper.tri(poset)][sample(choose(k, 2), nr_edges)] <- 1
    
    if(trans_reduced)
        poset = trans_reduction(poset)
    
    colnames(poset) <- paste("path", c(1:k), sep='')
    rownames(poset) <- paste("path", c(1:k), sep='')
    
    return(poset)
}



###############################################################################
#' @title Transforms a poset into its transitive reduction
#' 
#' @description \code{trans_reduction} transforms a poset, i.e. a binary 
#' traiangular matrix corresponding to a directed acyclic graph, into another
#' binary matrix representing its transitive reduction.
#' 
#' @param poset the input binary matrix representing the partial order between 
#' pathways
#' 
#' @details the transitive reduction of a directed acyclic graph is the 
#' directed acyclic graph with as few edges as possible that has the same 
#' reachability relation as the initial graph.
#' 
#' @return a binary matrix, representing the transitive reduction of the input
#' poset
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases trans_reduction
#' 
#' @export
trans_reduction <- function(poset)
{   
    old_names = dimnames(poset)
    
    # rename the events as "n", followed by index
    colnames(poset) =  rownames(poset) = paste("path", 1:nrow(poset), sep='')
    
    # transforms the matrix of the poset into a binary relation and extract its
    # transitive reduction, and then transforms it back into a matrix
    R <- as.relation(poset)
    RT = transitive_reduction(R)    
    res = relation_incidence(RT)
    
    # rearrange the rows and columns of the poset
    res = res[rownames(poset),colnames(poset)]
    
    dimnames(res) = old_names
    
    return(res)
}



###############################################################################
#' @title Plots a poset
#' 
#' @description \code{plotPoset} plots a poset given as input, by representing 
#' the causal relations between events. 
#' 
#' @param poset the poset to be plotted, given as a binary matrix
#' @param size the size of the names of events, default to 12
#' 
#' @details The column names of the poset should be the names of the events, 
#' and they will be used for plotting.
#' 
#' @return the plot of the input poset
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases plot_poset
#' 
#' @export
plotPoset <- function(poset, size=12)
{    
    require("graph")
    Names = colnames(poset)
    rownames(poset) <- Names
    am.graph <- new("graphAM", adjMat = poset, edgemode="directed")
    plot(am.graph, attrs = list( node = list(color = "transparent", fontsize = size, fontcolor="dodgerblue4"), 
                                 edge = list(arrowsize=0.5, color="antiquewhite4")))
}



###############################################################################
#' @title Plots a poset on pathways
#' 
#' @description \code{plotPosetPaths} plots a poset on the level of pathways 
#' given as input, by representing  the causal relations between pathways,
#' together with gene membership. 
#' 
#' @param poset the poset to be plotted, given as a binary matrix with 
#' relationships between pathways
#' @param groupsInput the assignment of genes to pathways
#' @param Datamat the alterations of genes, represented as a binary matrix
#' @param size the size of the names of events, default to 30
#' 
#' @details The column names of the poset should be the names of the events, 
#' and they will be used for plotting. This function only works for groups 
#' generated after the stochastic search step.
#' 
#' @return the plot of the input poset, together with assignment of genes to 
#' pathways
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases plot_poset
#' 
#' @export
plotPosetPaths<-function(poset,groupsInput,Datamat,size=30)
{
    require("graph")
    namesEvs<-colnames(Datamat)
    
    # remove the empty event at the end
    poset<-poset[-dim(poset)[1],-dim(poset)[1]]
    namesPos<-rep(0,dim(poset)[1])
    for (l in 1:length(groupsInput))
    {
        namesPos[l]<-paste(namesEvs[groupsInput[[l]]],collapse=", ")
    }
    colnames(poset)<-namesPos
    rownames(poset)<-namesPos
    am.graph <- new("graphAM", adjMat = poset, edgemode="directed")
    plot(am.graph, attrs = list( node = list(color = "transparent", fontsize = size, fontcolor="dodgerblue4"), 
                                 edge = list(arrowsize=0.5, color="antiquewhite4")))   
}



###############################################################################
#' @title Computes the topological sort of a poset
#' 
#' @description \code{topologicalSort} outputs the events in a poset according 
#' to their topological order.
#' 
#' @param poset the poset to be sorted
#' 
#' @details the topological sorting of a poset is useful in iteratively 
#' going over all vertices in a poset and making sure every vertex is only
#' visited once
#' 
#' @return a vector of events (vertices) topologically sorted
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases topologicalSort
#' 
#' @export
topologicalSort <- function(poset)
{
    p <- ncol(poset)
    sorted.list <- rep(0, p)
    for(i in 1:p)
    {
        sorted.list[i] = min.degree = which.min(apply(poset, 2, sum))  
        poset[min.degree, ] <- 0
        poset[, min.degree] <- Inf
    }
    
    return(sorted.list)    
}



###############################################################################
#' @title Simulates a linear poset
#' 
#' @description \code{make_linear_poset} simulates a linear poset, for a given
#' number of events.
#' 
#' @param p the number of events
#' 
#' @details A linear poset is used as starting solution for simulated annealing
#' for CBN structure optimization.
#' 
#' @return a binary matrix with relationships between events representing
#' order constraints as a chain (linear poset).
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @aliases make_linear_poset
#' 
#' @export
make_linear_poset <- function(p) {
    pos <- matrix(0, p, p)
    if(p > 1){
        for(i in 1:(p-1))
            pos[i, i+1] = 1
    }
    return(pos)
}