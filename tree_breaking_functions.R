#' Function to decide whether a branch can be broken
#'
#' @param b A vector of length 1 or 2 describing. Each element
#' represents a branch, the branches are siblings, and the number in
#' the vector gives the fraction of samples that contain descendants
#' of that branch.
#'
#' @return A Boolean vector of the same length as b, describing
#' whether the branch should be broken.
branch_breakings = function(b) {
    if(length(b) == 1) {
        ## if the branch doesn't have a sibling, we don't break
        return(rep(FALSE, length(b)))
    }
    if(b[1] != 1 && b[2] != 1) {
        ## If both branches don't have descendants in all of the
        ## samples, we don't break either of them
        return(c(FALSE, FALSE))
    } else if(b[1] != 1 && b[2] == 1) {
        ## If b2 has descendants in all the samples and b1 doesn't, we
        ## break at b1
        return(c(TRUE, FALSE))
    } else if(b[1] == 1 && b[2] != 1) {
        ## If b1 has descendants in all the samples and b2 doesn't, we
        ## break at b2
        return(c(FALSE, TRUE))
    } else if(b[1] == 1 && b[2] == 1) {
        ## If both branches have descendants in all the samples, we
        ## can break at either, here we arbitrarily choose the first
        ## one
        return(c(TRUE, FALSE))
    } else {
        ## branches with more than one sibling aren't supported, but
        ## we'll just print a warning and not break any of the
        ## branches
        warning("Trees with more than two descendants from a node are not supported.")
        return(rep(FALSE, length(b)))
    }
}


#' Finds all the edges in the tree that should be broken.
#'
#' @param tr The phylogenetic tree.
#' @param X A n x p matrix containing abundances of the species in the
#' samples. n is the number of samples, p is the number of species,
#' which should be the same as the number of leaves in the
#' phylogenetic tree.
#'
#' @return A boolean vector of length equal to the number of edges in
#' the tree, describing which edges can be broken off to form new
#' subtrees.
#' @importFrom dplyr group_by, mutate
#' @importFrom magrittr %>%
#' @importFrom igraph graph_from_edgelist
break_tree <- function(tr, X) {
    ## we have an old function that computes for each branch what
    ## fraction of the samples a descendant of that branch is present
    ## in.
    A = treeDA:::makeDescendantMatrix(tr)
    ## X %*% A gives node abundances. Node abundances are the same as
    ## branch abundances for the branch terminating in that node,
    ## hence taking the columns tr$edge[,2]. Finally, we query which
    ## branch abundances are positive, and take the average
    p = Matrix::colMeans((X %*% A[,tr$edge[,2]]) > 0)
    output_df = data.frame(tr$edge, frac = p)
    colnames(output_df)[1:2] = c("parent", "child")
    annotated_edges = output_df %>%
        dplyr::group_by(parent) %>%
        dplyr::mutate(break_branch = branch_breakings(frac))
    ## each broken branch will require a new node in the tree
    n_new_nodes = sum(annotated_edges$break_branch)
    broken_tree = as(annotated_edges[,1:2], "matrix")
    if(n_new_nodes >= 1) {
        broken_tree[annotated_edges$break_branch,1] = (max(annotated_edges) + 1):(max(annotated_edges) + n_new_nodes)
    }
    return(igraph::graph_from_edgelist(broken_tree))
    
}


# copied from the phyloseq function tree_layout, tree must be in postorder
get_node_positions <- function(tr) {
    Nedge = nrow(tr$edge)[1]
    Nnode = tr$Nnode
    Ntip = length(tr$tip.label)
    ROOT = Ntip + 1
    TIPS = tr$edge[(tr$edge[,2] <= Ntip), 2]
    NODES = (ROOT):(Ntip + Nnode)
    nodelabels = tr$node.label
    xx = phyloseq:::ape_node_depth_edge_length(Ntip, Nnode, tr$edge, Nedge, tr$edge.length)
    yy = numeric(Ntip + Nnode)
    yy[TIPS] = 1:Ntip
    ape_node_height <- function(Ntip, Nnode, edge, Nedge, yy) {
        .C(ape:::node_height, PACKAGE = "ape", as.integer(Ntip), 
           as.integer(Nnode), as.integer(edge[, 1]),
           as.integer(edge[, 2]), as.integer(Nedge), as.double(yy))[[6]]
    }
    yy = ape_node_height(Ntip, Nnode, tr$edge, Nedge, yy)
    return(cbind(x = yy, y = -xx))
}

# tree must be in postorder
get_uf_contribs <- function(tr, v1, v2, include.lengths = FALSE) {
    tipComp = apply(cbind(v1, v2), 1, function(x) {
        if(x[1] != 0 & x[2] == 0)
            return("a")
        else if(x[1] == 0 & x[2] != 0)
            return("b")
        else if(x[1] != 0 & x[2] != 0)
            return("m")
        return("none")
    })
    overallComp = tipComp[tr$edge[,2]]
    for(i in 1:nrow(tr$edge)) {
        if(!is.na(overallComp[i]))
            next
        node = tr$edge[i,2]
        children = overallComp[which(tr$edge[,1] == node)]
        if(children[1] == children[2])
            overallComp[i] = children[1]
        else if(any(children == "m"))
            overallComp[i] = "m"
        else if(all(children == c("a", "b")) | all(children == c("b", "a")))
            overallComp[i] = "m"
        else if(any(children == "none"))
            overallComp[i] = children[children != "none"]
    }
    pure = sapply(overallComp, function(x) {
        if(x == "a" | x == "b")
            return(1)
        return(0)
    })
    if(include.lengths) {
        contributions = pure * tr$edge.length
    } else {
        contributions = pure
    }
    
    return(list(tree = tr, edgeContributions = contributions / sum(contributions),
                types = overallComp))
}


# tree must be in postorder
get_wuf_contribs <- function(tr, v1, v2, include.lengths = FALSE) {
    tipFrac1 = v1 / sum(v1)
    tipFrac2 = v2 / sum(v2)
    overallFrac1 = tipFrac1[tr$edge[,2]]
    overallFrac2 = tipFrac2[tr$edge[,2]]
    for(i in 1:nrow(tr$edge)) {
        if(!is.na(overallFrac1[i]))
            next
        node = tr$edge[i,2]
        childrenIdx = which(tr$edge[,1] == node)
        overallFrac1[i] = sum(overallFrac1[childrenIdx])
        overallFrac2[i] = sum(overallFrac2[childrenIdx])
    }
    allWeights = abs(overallFrac1 - overallFrac2)
    if(include.lengths) {
        contribution = allWeights * tr$edge.length
    } else {
        contribution = allWeights
    }
    contribution = contribution / sum(contribution)
    return(list(tree = tr, edgeContributions = contribution))
}

get_guf_contribs <- function(tr, v1, v2, alpha, include.lengths = FALSE) {
    tipFrac1 = v1 / sum(v1)
    tipFrac2 = v2 / sum(v2)
    overallFrac1 = tipFrac1[tr$edge[,2]]
    overallFrac2 = tipFrac2[tr$edge[,2]]
    for(i in 1:nrow(tr$edge)) {
        if(!is.na(overallFrac1[i]))
            next
        node = tr$edge[i,2]
        childrenIdx = which(tr$edge[,1] == node)
        overallFrac1[i] = sum(overallFrac1[childrenIdx])
        overallFrac2[i] = sum(overallFrac2[childrenIdx])
    }
    numerator = function(p, alpha) {
        if(length(p) != 2) {
            stop("p needs to be a length 2 vector")
        }
        if(sum(p) == 0) {
            return(0)
        }
        return((p[1] + p[2])^alpha * abs((p[1] - p[2]) / (p[1] + p[2])))
    }
    allWeights = apply(cbind(overallFrac1, overallFrac2), 1, numerator, alpha)
    if(include.lengths) {
        contribution = allWeights * tr$edge.length
    } else {
        contribution = allWeights
    }
    contribution = contribution / sum(contribution)
    return(list(tree = tr, edgeContributions = contribution))
}


# tree must be in postorder
get_ndescendants <- function(tr) {
    ndescendants = rep(1, length(tr$tip.label))[tr$edge[,2]]
    for(i in 1:nrow(tr$edge)) {
        if(!is.na(ndescendants[i]))
            next
        node = tr$edge[i,2]
        childrenIdx = which(tr$edge[,1] == node)
        ndescendants[i] = sum(ndescendants[childrenIdx])
    }
    return(ndescendants)
}

plotEdgeContribs <- function(tr, c, ...) {
    p = plot_tree(tr)
    pdata = as.data.frame(p$data)
    contribDF = data.frame(tr$edge)
    contribDF$contrib = c
    contribDF$edgeID = apply(contribDF, 1, function(x)
        paste(as.numeric(x[1]), as.numeric(x[2]), sep = "_"))
    pdata$edgeID = apply(pdata, 1, function(x)
        paste(as.numeric(x[1]), as.numeric(x[2]), sep = "_"))
    merged = merge(pdata, contribDF, by = "edgeID")
    plot = p + 
        geom_segment(aes(x = xleft, y = y, xend = xright, yend = y,
                         color = contrib), data = merged) +
                             scale_color_viridis(...) +
                                 coord_flip() + scale_x_reverse()
    return(plot)
}

rv_coef <- function(X1, X2) {
    X1 = scale(X1, scale = FALSE)
    X2 = scale(X2, scale = FALSE)
    s12 = sum(diag(t(X1) %*% X2 %*% t(X2) %*% X1))
    s11 = sum(diag(X1 %*% t(X1)))
    s22 = sum(diag(X2 %*% t(X2)))
    return(s12 / (s11 * s22))
}


dpcoa <- function(otutab, tree, k) {
    Q = ape::vcv(tree)
    Qeig = eigen(Q)
    X = otutab
    D = rowSums(X)
    for(i in 1:nrow(X)) {
        X[i,] = X[i,] / sum(X[i,])
    }
    X = scale(X, scale = FALSE)
    out.gpca = adaptiveGPCA:::gpcaEvecs(X, evecs = Qeig$vectors, evals= Qeig$values, D = D, k = 2)
    return(out.gpca)
}

dpcoaDist <- function(otutab, tree, r = 1) {
    Q = ape::vcv(tree)
    Qeig = eigen(Q)
    X = otutab
    D = rowSums(X)
    for(i in 1:nrow(X)) {
        X[i,] = X[i,] / sum(X[i,])
    }
    X = scale(X, scale = FALSE)
    if(r < 1 && r > 0) {
        out.gpca = adaptiveGPCA:::gpcaEvecs(X,
            evecs = Qeig$vectors,
            evals = (rep(r^(-1), ncol(X)) + (1-r)^(-1) * Qeig$values^(-1))^(-1),
            D = D, k = min(nrow(X), ncol(X)))
    } else if(r == 1) {
        out.gpca = adaptiveGPCA:::gpcaEvecs(X,
            evecs = Qeig$vectors,
            evals = Qeig$values,
            D = D, k = min(nrow(X), ncol(X)))
    } else if(r == 0){
        out.gpca = adaptiveGPCA:::gpcaEvecs(X,
            evecs = Qeig$vectors,
            evals = rep(1, ncol(X)),
            D = D, k = min(nrow(X), ncol(X)))
    }
    Uscaled = sweep(out.gpca$U, 2, STATS = out.gpca$lambda, FUN = "*")
    return(dist(Uscaled, method = "euclidean"))
}


tip_glom_k <- function (physeq, k = ntaxa(physeq), hcfun = agnes, ...)  {
    dd = as.dist(cophenetic.phylo(phy_tree(physeq)))
    psclust = cutree(as.hclust(hcfun(dd, ...)), k = k)
    cliques = levels(factor(psclust))[tapply(psclust, factor(psclust), 
        function(x) {
            length(x) > 1
        })]
    for (i in cliques) {
        physeq = merge_taxa(physeq, eqtaxa = names(psclust)[psclust == i])
    }
    return(physeq)
}

makeGlommedList <- function(phy, kvec) {
    glommedList = lapply(kvec, function(k) tip_glom_k(phy, k = k))
}

makeGlommingScores <- function(glommedList, kvec, method, d) {
    ordList = lapply(glommedList, function(p) ordinate(p, method = method, distance = d))
    for(i in 2:length(kvec)) {
        r = adaptiveGPCA:::findReflection(ordList[[i-1]]$vectors[,1:2], ordList[[i]]$vectors[,1:2])
        ordList[[i]]$vectors = sweep(ordList[[i]]$vectors, 2, STATS = r, FUN = "*")
    }
    scoresList = lapply(ordList, function(o) data.frame(o$vectors[,1:2], sample_data(glommedList[[1]])))
    scores = Reduce(rbind, scoresList)
    scores$ntaxa = rep(kvec, each = nsamples(glommedList[[1]]))
    return(scores)
}

makeGlommingDistances <- function(glommedList, kvec, d) {
    distList = lapply(glommedList, function(p) distance(p, method = d))
    n = nsamples(glommedList[[1]])
    distanceArray = array(dim = c(n, n, length(distList)))
    for(i in 1:length(distList)) {
        distanceArray[,,i] = as.matrix(distList[[i]])
    }
    out.distatis = distatis(distanceArray)
    colnames(out.distatis$res4Cmat$C) = rownames(out.distatis$res4Cmat$C) = kvec
    return(out.distatis)
}

makeGlommingDistancesDPCoA <- function(glommedList, kvec, r = 1) {
    n = nsamples(glommedList[[1]])
    distList = lapply(glommedList, function(g) {
        if(taxa_are_rows(g)) {
            X = t(as(otu_table(g), "matrix"))
        } else {
            X = as(otu_table(g), "matrix")
        }
        dpcoaDist(X, phy_tree(g), r = r)
    })
    distanceArray = array(dim = c(n, n, length(distList)))
    for(i in 1:length(distList)) {
        distanceArray[,,i] = as.matrix(distList[[i]])
    }
    out.distatis = distatis(distanceArray)
    colnames(out.distatis$res4Cmat$C) = rownames(out.distatis$res4Cmat$C) = kvec
    return(out.distatis)
}

makeGlommingDistancesGUF <- function(glommedList, kvec, alpha) {
    distList = lapply(glommedList, function(p) fastUniFrac(p, alpha = alpha, weighted = TRUE))
    n = nsamples(glommedList[[1]])
    distanceArray = array(dim = c(n, n, length(distList)))
    for(i in 1:length(distList)) {
        distanceArray[,,i] = as.matrix(distList[[i]])
    }
    out.distatis = distatis(distanceArray)
    colnames(out.distatis$res4Cmat$C) = rownames(out.distatis$res4Cmat$C) = kvec
    return(out.distatis)
   
}

makeGlommingScoresDPCoA <- function(glommedList, kvec) {
    dpcoaList = lapply(glommedList, function(g) {
        if(taxa_are_rows(g)) {
            X = t(as(otu_table(g), "matrix"))
        } else {
            X = as(otu_table(g), "matrix")
        }
        dpcoaDist(X, phy_tree(g))
    })
    ## find the sign changes for the axes that make the points
    ## align the best from one plot to another
    for(i in 2:length(kvec)) {
        r = adaptiveGPCA:::findReflection(dpcoaList[[i-1]][,1:2], dpcoaList[[i]][,1:2])
        dpcoaList[[i]] = sweep(dpcoaList[[i]], 2, STATS = r, FUN = "*")
    }
    scoresList = lapply(dpcoaList, function(o)
        data.frame(o[,1:2], sample_data(glommedList[[1]])))
    scores = Reduce(rbind, scoresList)
    scores$ntaxa = rep(kvec, each = nsamples(glommedList[[1]]))
    colnames(scores)[1:2] = c("Axis.1", "Axis.2")
    return(scores)
}

# copied from ape but with the collapse.singles line deleted
my.drop.tip <- 
function (phy, tip, trim.internal = TRUE, subtree = FALSE, root.edge = 0, 
    rooted = is.rooted(phy), interactive = FALSE) 
{
    if (!inherits(phy, "phylo")) 
        stop("object \"phy\" is not of class \"phylo\"")
    Ntip <- length(phy$tip.label)
    if (interactive) {
        cat("Left-click close to the tips you want to drop; right-click when finished...\n")
        xy <- locator()
        nToDrop <- length(xy$x)
        tip <- integer(nToDrop)
        lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        for (i in 1:nToDrop) {
            d <- sqrt((xy$x[i] - lastPP$xx)^2 + (xy$y[i] - lastPP$yy)^2)
            tip[i] <- which.min(d)
        }
    }
    else {
        if (is.character(tip)) 
            tip <- which(phy$tip.label %in% tip)
    }
    out.of.range <- tip > Ntip
    if (any(out.of.range)) {
        warning("some tip numbers were larger than the number of tips: they were ignored")
        tip <- tip[!out.of.range]
    }
    if (!length(tip)) 
        return(phy)
    if (length(tip) == Ntip) {
        warning("drop all tips of the tree: returning NULL")
        return(NULL)
    }
    wbl <- !is.null(phy$edge.length)
    if (!rooted && subtree) {
        phy <- root(phy, (1:Ntip)[-tip][1])
        root.edge <- 0
    }
    phy <- reorder(phy)
    NEWROOT <- ROOT <- Ntip + 1
    Nnode <- phy$Nnode
    Nedge <- dim(phy$edge)[1]
    if (subtree) {
        trim.internal <- TRUE
        tr <- reorder(phy, "postorder")
        N <- .C(node_depth, as.integer(Ntip), as.integer(Nnode), 
            as.integer(tr$edge[, 1]), as.integer(tr$edge[, 2]), 
            as.integer(Nedge), double(Ntip + Nnode), 1L)[[6]]
    }
    edge1 <- phy$edge[, 1]
    edge2 <- phy$edge[, 2]
    keep <- !logical(Nedge)
    keep[match(tip, edge2)] <- FALSE
    if (trim.internal) {
        ints <- edge2 > Ntip
        repeat {
            sel <- !(edge2 %in% edge1[keep]) & ints & keep
            if (!sum(sel)) 
                break
            keep[sel] <- FALSE
        }
        if (subtree) {
            subt <- edge1 %in% edge1[keep] & edge1 %in% edge1[!keep]
            keep[subt] <- TRUE
        }
        if (root.edge && wbl) {
            degree <- tabulate(edge1[keep])
            if (degree[ROOT] == 1) {
                j <- integer(0)
                repeat {
                  i <- which(edge1 == NEWROOT & keep)
                  j <- c(i, j)
                  NEWROOT <- edge2[i]
                  degree <- tabulate(edge1[keep])
                  if (degree[NEWROOT] > 1) 
                    break
                }
                keep[j] <- FALSE
                if (length(j) > root.edge) 
                  j <- 1:root.edge
                NewRootEdge <- sum(phy$edge.length[j])
                if (length(j) < root.edge && !is.null(phy$root.edge)) 
                  NewRootEdge <- NewRootEdge + phy$root.edge
                phy$root.edge <- NewRootEdge
            }
        }
    }
    if (!root.edge) 
        phy$root.edge <- NULL
    phy$edge <- phy$edge[keep, ]
    if (wbl) 
        phy$edge.length <- phy$edge.length[keep]
    TERMS <- !(phy$edge[, 2] %in% phy$edge[, 1])
    oldNo.ofNewTips <- phy$edge[TERMS, 2]
    if (subtree) {
        i <- which(tip %in% oldNo.ofNewTips)
        if (length(i)) {
            phy$tip.label[tip[i]] <- "[1_tip]"
            tip <- tip[-i]
        }
    }
    n <- length(oldNo.ofNewTips)
    phy$edge[TERMS, 2] <- rank(phy$edge[TERMS, 2])
    phy$tip.label <- phy$tip.label[-tip]
    if (subtree || !trim.internal) {
        node2tip <- oldNo.ofNewTips[oldNo.ofNewTips > Ntip]
        new.tip.label <- if (subtree) {
            paste("[", N[node2tip], "_tips]", sep = "")
        }
        else {
            if (is.null(phy$node.label)) 
                rep("NA", length(node2tip))
            else phy$node.label[node2tip - Ntip]
        }
        phy$tip.label <- c(phy$tip.label, new.tip.label)
    }
    phy$Nnode <- dim(phy$edge)[1] - n + 1L
    newNb <- integer(Ntip + Nnode)
    newNb[NEWROOT] <- n + 1L
    sndcol <- phy$edge[, 2] > n
    newNb[sort(phy$edge[sndcol, 2])] <- (n + 2):(n + phy$Nnode)
    phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]]
    phy$edge[, 1] <- newNb[phy$edge[, 1]]
    storage.mode(phy$edge) <- "integer"
    if (!is.null(phy$node.label)) 
        phy$node.label <- phy$node.label[which(newNb > 0) - Ntip]
    return(phy)
}



################################################################################
# Fast Generalized UniFrac for R.
# Very slightly adapted from Fast Unifrac in phyloseq (github.com/joey711/phyloseq)
################################################################################
#' @param physeq A phyloseq object.
#' @param weighted If TRUE, perform generalized weighted Unifrac with parameter alpha. If FALSE, perform unweighted Unifrac.
#' @param alpha The parameter for generalized weighted Unifrac.
#' @param normalized If TRUE, computed normalized weighted Unifrac
#' distances. Unweighted Unifrac distances are always normalized, as
#' are generalized Unifrac distances.
#' @param parallel If TRUE, perform some computations in parallel.
#' @importFrom ape prop.part
#' @importFrom ape reorder.phylo
#' @importFrom ape node.depth
#' @importFrom ape node.depth.edgelength
#' @keywords internal
#' @import foreach
fastUniFrac <- function(physeq, alpha=1, weighted=TRUE, normalized=TRUE, parallel=FALSE){
	# Access the needed components. Note, will error if missing in physeq.
	OTU  <- otu_table(physeq)
	tree <- phy_tree(physeq)
	# Some important checks.
	if( is.null(tree$edge.length) ) {
	  stop("Tree has no branch lengths, cannot compute UniFrac")
	}
	if( !is.rooted(tree) ) {
	  stop("Rooted phylogeny required for UniFrac calculation")
	}
	### Some parallel-foreach housekeeping.
	# If user specifies not-parallel run (the default), register the sequential "back-end"
	if( !parallel ){ registerDoSEQ() }
	# create N x 2 matrix of all pairwise combinations of samples.
	spn <- combn(sample_names(physeq), 2, simplify=FALSE)
	# Make sure OTU is in species-are-rows orientation
	if( !taxa_are_rows(physeq) ){OTU <- t(OTU)}
  # Convert to standard matrix
	OTU <- as(OTU, "matrix")  
	# Enforce that tree and otu_table indices are the same order, 
	# by re-ordering OTU, if needed
	if( !all(rownames(OTU) == taxa_names(tree)) ){
	  OTU <- OTU[taxa_names(tree), ]
	}
	########################################
	# Build the requisite matrices as defined 
	# in the Fast UniFrac article.
	########################################
	## This only needs to happen once in a call to UniFrac.
	## Notice that A and B do not appear in this section.
	# Begin by building the edge descendants matrix (edge-by-sample)
  # `edge_array`
  #
	# Create a list of descendants, starting from the first internal node (root)
	ntip <- length(tree$tip.label)
	if(ntip != ntaxa(physeq)) stop("Incompatible tree and OTU table!")
	# Create a matrix that maps each internal node to its 2 descendants
	# This matrix doesn't include the tips, so must use node#-ntip to index into it
	node.desc <- matrix(tree$edge[order(tree$edge[,1]),][,2],byrow=TRUE,ncol=2)
	# Define the edge_array object
	# Right now this is a node_array object, each row is a node (including tips)
	# It will be subset and ordered to match tree$edge later
	edge_array <- matrix(0, nrow=ntip+tree$Nnode, ncol=nsamples(physeq), 
	                     dimnames=list(NULL, sample_names=sample_names(physeq)))
	# Load the tip counts in directly
	edge_array[1:ntip,] <- OTU
	# Get a list of internal nodes ordered by increasing depth
	ord.node <- order(node.depth(tree))[(ntip+1):(ntip+tree$Nnode)]
	# Loop over internal nodes, summing their descendants to get that nodes count
	for(i in ord.node){
	  edge_array[i,] <- colSums(edge_array[node.desc[i-ntip,], , drop=FALSE], na.rm = TRUE)
	}
	# Keep only those with a parental edge (drops root) and order to match tree$edge
	edge_array <- edge_array[tree$edge[,2],]
	# Remove unneeded variables.
	rm(node.desc)
	# If unweighted-UniFrac, coerce to a presence-absence contingency, occ
	if(!weighted){
		# For unweighted UniFrac, convert the edge_array to an occurrence (presence/absence binary) array
		edge_occ <- (edge_array > 0) - 0
	}
	if( weighted & normalized ){
		# This is only relevant to weighted-UniFrac.
		# For denominator in the normalized distance, we need the age of each tip.
	  # 'z' is the tree in postorder order used in calls to .C
	  # Descending order of left-hand side of edge (the ancestor to the node)
	  z = reorder.phylo(tree, order="postorder")
	  # Call phyloseq-internal function that in-turn calls ape's internal
	  # horizontal position function, in C, using the re-ordered phylo object, `z`
	  tipAges = node.depth.edgelength(tree)
	  # Keep only the tips, and add the tip labels in case `z` order differs from `tree`
	  tipAges <- tipAges[1:length(tree$tip.label)]
	  names(tipAges) <- z$tip.label	
    # Explicitly re-order tipAges to match OTU
	  tipAges <- tipAges[rownames(OTU)]
	}
########################################	
        ## optionally-parallel implementation with foreach
########################################
	samplesums = sample_sums(physeq)
	distlist <- foreach( i = spn, .packages="phyloseq") %dopar% {
            A  <- i[1]
            B  <- i[2]
            AT <- samplesums[A]
            BT <- samplesums[B]
            if( weighted ){
                ## weighted UniFrac
                wUF_branchweight <- abs(edge_array[, A]/AT - edge_array[, B]/BT) *
                    (edge_array[, A]/AT + edge_array[, B]/BT)^(alpha-1)
                ## calculate the w-UF numerator
                numerator <- sum({tree$edge.length * wUF_branchweight}, na.rm = TRUE)
                ## if not-normalized weighted UniFrac, just return "numerator";
                ## the u-value in the w-UniFrac description
                if(!normalized){
                    return(numerator)
                } else {
                    ## denominator, changed from the phyloseq implementation
                    denominator <-
                        sum({tree$edge.length * (edge_array[,A] / AT + edge_array[,B] / BT)^alpha}, na.rm = TRUE)
                    ## return the normalized weighted UniFrac values
                    return(numerator / denominator)
                }
            } else {
                ## Unweighted UniFrac
                ## Subset matrix to just columns A and B
                edge_occ_AB <- edge_occ[, c(A, B)]
                ## Keep only the unique branches. Sum the lengths
                edge_uni_AB_sum <- sum((tree$edge.length * edge_occ_AB)[rowSums(edge_occ_AB, na.rm=TRUE) < 2, ], na.rm=TRUE)
                ## Normalize this sum to the total branches among these two samples, A and B
                uwUFpairdist <- edge_uni_AB_sum / sum(tree$edge.length[rowSums(edge_occ_AB, na.rm=TRUE) > 0])
                return(uwUFpairdist)
            }
	}
	## Initialize UniFracMat with NAs
	UniFracMat <- matrix(NA_real_, nsamples(physeq), nsamples(physeq))
	rownames(UniFracMat) <- colnames(UniFracMat) <- sample_names(physeq)
  # Matrix-assign lower-triangle of UniFracMat. Then coerce to dist and return.
  	matIndices <- do.call(rbind, spn)[, 2:1]
  	# Take care of edge case where there are two samples -> 1 pair of indices -> rbind doesn't return a matrix
  	if(!is.matrix(matIndices)) matIndices <- matrix(matIndices, ncol=2)
	UniFracMat[matIndices] <- unlist(distlist)
	return(as.dist(UniFracMat))	
}
################################################################################
