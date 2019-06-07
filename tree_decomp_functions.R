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
#' @return A graph describing the forest resulting from breaking the
#' trees.
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
    new_roots = (max(annotated_edges) + 1):(max(annotated_edges) + n_new_nodes)
    if(n_new_nodes >= 1) {
        broken_tree[annotated_edges$break_branch,1] = new_roots
    }
    return(igraph::graph_from_edgelist(broken_tree))
}

#' Computes the number of leaves per tree
#'
#' Given an igraph-formatted forest, computes the number of leaves for
#' each tree in the forest.
#'
#' @param forest An igraph object describing the forest. The function
#' assumes that the graph is directed, so that the leaves are the
#' nodes in the graph that have no outgoing edges.
#'
#' @param return A vector, length equal to the number of trees in the
#' forest, elements equal to the number of leaves in the corresponding
#' tree.
get_leaf_sizes <- function(forest) {
    cc = components(forest)
    adj_mat = as_adjacency_matrix(forest)
    tree_sizes = adply(1:cc$no, 1, function(comp_idx) {
        nodes = which(cc$membership == comp_idx)
        adj_submat = adj_mat[nodes,nodes]
        n_leaves = sum(rowSums(adj_submat) == 0)
        return(c("n_leaves" = n_leaves))
    }, .id = "tree_id")
    return(tree_sizes)
}


#' Partitions the leaves by tree.
#'
#' Given an igraph-formatted forest,
#'
#' @param forest An igraph object describing the forest. The function
#' assumes that the graph is directed, so that the leaves are the
#' nodes in the graph that have nooutgoing edges.
#'
#' @param return A list, length equal to the number of trees in the
#' forest, each element a vector of nodes giving the leaf nodes in the
#' corresponding tree.
get_leaf_list <- function(forest) {
    cc = components(forest)
    adj_mat = as_adjacency_matrix(forest)
    leaves = which(rowSums(adj_mat) == 0)
    node_and_component = data.frame(node = 1:length(cc$membership), component = cc$membership)
    leaf_components = node_and_component[leaves,]
    leaf_list = split(leaf_components$node, leaf_components$component)
    return(leaf_list)
}

#' Compute branch probabilities.
#'
#' @param X An otu abundance matrix, n x p, rows are samples.
#' @param tr A phylogenetic tree with p leaves.
#'
#' @return A n x n_branches matrix containing branch probabilities for
#' each sample/branch combination.
get_branch_probs <- function(X, tr) {
    A = treeDA:::makeDescendantMatrix(tr)
    P = plyr::aaply(X, 1, function(x) return(x / sum(x)))
    ## P %*% A are the node contributions. The contribution of a
    ## branch is the same as the contribution of the child node in
    ## that branch, hence the following:
    child_node = tr$edge[,2]
    return(P %*% A[,child_node])
}

#' Computes branch contributions for a set of unifrac distances
#'
#' @param X An otu abundance matrix, n x p.
#' @param tr A phylogenetic tree with p leaves.
#' @param alpha_list A list containing the alpha values for
#' guf. Numeric elements correspond to values of alpha, and a
#' character element "unweighted" corresponds to unweighted Unifrac.
#'
#' @return 
get_all_branch_contribs <- function(X, tr, alpha_list, per_unit_branch = TRUE) {
    P = get_branch_probs(X, tr)
    branch_lengths = tr$edge.length
    all_contribs = plyr::laply(alpha_list, function(a) {
        if(a == "unweighted") {
            return(uf_contribs_from_branch_probs(P, branch_lengths, per_unit_branch))
        } else if(is.numeric(a) && a >= 0 && a <= 1) {
            return(guf_contribs_from_branch_probs(P, branch_lengths, a, per_unit_branch))
        } else {
            warning("alpha not in [0,1] or equal to 'unweighted'")
            return(NULL)
        }        
    })
    
}


#' Generalized Unifrac branch contributions
#'
#' Computes generalized Unifrac branch contributions with pre-computed
#' branch probabilities.
#'
#' @param P A n x nbranches matrix with branch probabilities for each
#' branch and each sample.
#' @param branch_lengths An nbranches-length vector with branch lengths.
#' @param alpha Alpha parameter for generalized Unifrac.
#'
#' @return A vector of length nbranches giving branch contributions.
guf_contribs_from_branch_probs <- function(P, branch_lengths, alpha, per_unit_branch = TRUE) {
    n = nrow(P)
    pairs = combn(n,2)
    branch_contrib =
        apply(pairs, 2, function(idx) {
            p1 = P[idx[1],]
            p2 = P[idx[2],]
            guf_dist = sum(branch_lengths * (p1 + p2)^alpha *
                               abs((p1 - p2) / (p1 + p2)), na.rm = TRUE)
            contrib = (p1 + p2)^alpha * abs((p1 - p2) / (p1 + p2)) / guf_dist
            if(!per_unit_branch) {
                contrib = contrib * branch_lengths
            }
            contrib[(p1 + p2 == 0)] = 0
            if(guf_dist == 0) contrib = rep(0, length(contrib))
            return(contrib)
        })
    return(rowMeans(branch_contrib))
}

#' Unweighted Unifrac branch contributions.
#'
#' Computes unweighted Unifrac branch contributions from pre-computed
#' branch probabilities.
#'
#' @param P An n x nbranches matrix with branch probabilities for each
#' branch and each sample.
#' @param branch_lengths An nbranches-length vector with branch
#' lengths.
#'
#' @return A vector of length nbranches giving branch contributions.
uf_contribs_from_branch_probs <- function(P, branch_lengths, per_unit_branch = TRUE) {
    P_ind = P > 0
    ## for all pairs of indices
    n = nrow(P)
    pairs = combn(n, 2)    
    branch_contrib =
        apply(pairs, 2, function(idx) {
            p1 = P_ind[idx[1],]
            p2 = P_ind[idx[2],]
            uf_dist = sum(branch_lengths * abs(p1 - p2)) / sum(branch_lengths)
            if(is.na(uf_dist)) browser()
            if(uf_dist == 0) return(rep(0, length(p1)))
            if(per_unit_branch) {
                return(abs(p1 - p2) / uf_dist)
            } else {
                return(branch_lengths * abs(p1 - p2) / uf_dist)
            }
        })
    return(rowMeans(branch_contrib))
}

#' Returns the number of descendants for each branch in the tree.
#'
#' @param tr A tree.
#'
#' @return A vector, length equal to the number of branches in the
#' tree, giving the number of descendants of that branch.
get_ndescendants <-  function(tr) {
    A = treeDA:::makeDescendantMatrix(tr)
    ## The number of descendants of a branch is the number of
    ## descendants of the child node of that branch, which is the
    ## second column of tr$edge
    return(Matrix::colSums(A[,tr$edge[,2]]))
}


#' Makes accumulation plot for branch contributions.
#' 
#' @param contribs A matrix with columns corresponding to branches.
#' @param tr A phylogenetic tree, from ape
#' @param contrib_names A character vector, giving the names for the
#' types corresponding to the rows of contribs.
contrib_accumulation_plot_old <- function(contribs, tr, contrib_names) {
    desc = get_ndescendants(tr)
    breaks = unique(desc)
    out = sapply(breaks, function(b) rowSums(contribs[,desc <= b]) / rowSums(contribs))
    contribs_and_breaks = data.frame(t(out), breaks)
    plotting_df = melt(contribs_and_breaks, id.vars = "breaks")
    levels(plotting_df$variable) = contrib_names
    attributes(plotting_df$variable)$class = c("ordered", "factor")
    return(ggplot(plotting_df, aes(x = breaks, y = value, color = variable)))
}

#' Makes accumulation plot for branch contributions.
#'
#' @param contribs A matrix with columns corresponding to branches.
#' @param tr A phylogenetic tree, from ape
#' @param contrib_names A character vector, giving the names for the
#' types corresponding to the rows of contribs.
contrib_accumulation_plot = function(contribs, tr, contrib_names) {
    desc = get_ndescendants(tr)
    breaks = unique(desc)
    out = sapply(breaks, function(b) rowSums(contribs[,desc <= b]) / rowSums(contribs))
    contribs_and_breaks = data.frame(t(out),
        proportion_of_branches = sapply(breaks, function(b) mean(desc <= b)))
    plotting_df = melt(contribs_and_breaks, id.vars = "proportion_of_branches")
    levels(plotting_df$variable) = contrib_names
    attributes(plotting_df$variable)$class = c("ordered", "factor")
    return(ggplot(plotting_df, aes(x = proportion_of_branches, y = value, color = variable)))
}

#' Tree plotting
#'
#' Idea is that plotting goes recursively: If you know where a node
#' is, its descendants are offset from that position and one level
#' down, with the offset computed to leave room for all of the
#' descendants of the child nodes.
#'
#' We need this function because the ape plotting functions don't work on trees 
my_plot_tree <- function(tr) {
    ntips = length(tr$tip.label)
    nnodes = tr$Nnode
    v = numeric(ntips + nnodes)
    d = numeric(ntips + nnodes)
    children = makeChildMapping(tr)
    ## step 1: assign values v to all the leaves, set number of
    ## descendants for leaves equal to 1
    tr_collapsed = collapse.singles(tr)
    Q = ape::vcv(tr_collapsed)
    v[1:ntips] = svd(Q, nu = 0, nv = 1)$v[,1]
    d[1:ntips] = 1
    ## step 2: for every node, find the number of descendants and the
    ## average v for the leaves descending from that node
    getDescendants <- function(node) {
        c = children[[node]]
        if(is.null(c))
            return(1)
        else {
            ndescendants = sum(sapply(c, getDescendants))
            d[node] <<- ndescendants
        }
    }
    getValues <- function(node) {
        c = children[[node]]
        if(is.null(c)) {
            return(v[node])
        }
        else {
            values = sapply(c, getValues)
            v[node] <<- sum(values * d[c]) / sum(d[c])
        }
    }
    root = findRoot(tr)
    getDescendants(root)
    getValues(root)

    nsegments = 2 * nnodes + ntips - 1
    npoints = nnodes + ntips
    ## step 3: tree traversed starting from the root. We find the two
    ## children of the current node. Base case is that we're at a
    ## leaf, in which case we just put a point at the position
    ## given. If we're at an internal node, we need to find the
    ## children, find the x and y positions for the children, put a
    ## point at the node and draw three segments (one horizontal and
    ## two vertical). 
    plot_subtree <- function(node, position) {
        c = children[[node]]
        x = position[1]
        y = position[2]
        if(is.null(c)) {
            pointDF = data.frame(x = position[1], y = position[2], node = node)
            segmentDF = NULL
            return(list(pointDF = pointDF, segmentDF = segmentDF))
        }
        ## put the children in order of their values
        c = c[order(v[c])]
        ## compute positions of child nodes
        xc = position[1] + get_child_offsets(d[c])
        yc = rep(position[2] - 1, length(c))
        pointList = list()
        segmentList = list()
        ## subtree plotting frames for each of the child nodes
        for(i in seq_along(c)) {
            subtree_plot = plot_subtree(c[i], c(xc[i], yc[i]))
            pointList[[i]] = subtree_plot$pointDF
            segmentList[[i]] = subtree_plot$segmentDF
        }        
        ### start function to combine subtree data frames ###
        subtree_frames = combine_subtree_frames(pointList, segmentList,
            x, y, node,
            xc, yc, c)
        return(subtree_frames)
    }
    plottingDFs = plot_subtree(root, c(0,0))
    ggplot(plottingDFs$segmentDF) +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend))
    return(plottingDFs)
}

#' Combines frames for subtrees and adds connecting points
#'
#' @param pointList A list of data frames for plotting points in each
#' of the child subtrees.
#' @param segmentList A list of data frames for plotting segmenst in
#' each of the child subtrees.
#' @param xroot The x position of the node connecting the subtrees.
#' @param yroot The y position of the node connecting the subtrees.
#' @param root_name
#' @param xc The x positions of the roots of the child subtrees.
#' @param yc The y positions of the roots of the child subtrees.
#' @param child_names
#'
#' @return A list with pointDF and segmentDF, for plotting the tree.
combine_subtree_frames <- function(pointList, segmentList,
                                   xroot, yroot, root_name,
                                   xc, yc, child_names) {
    pointDF = Reduce(rbind, pointList)
    if(all(sapply(segmentList, is.null))) {
        segmentDF = data.frame(x = numeric(0),
            y = numeric(0),
            xend = numeric(0),
            yend = numeric(0),
            type = character(0),
            stringsAsFactors = FALSE)
    } else {
        segmentDF = Reduce(rbind, segmentList)
    }
    npoints = nrow(pointDF)
    nsegs = nrow(segmentDF)
    ## add the base point for this tree
    pointDF[npoints + 1,1:2] = c(xroot, yroot)
    pointDF$node[npoints + 1] = root_name
    ## horizontal segment
    segmentDF[nsegs + 1,1:4] = c(xc[1], yroot, xc[length(xc)], yroot)
    segmentDF$type[nsegs + 1] = "horizontal"
    ## vertical segments
    for(i in seq_along(child_names)) {
        segmentDF[nsegs + 1 + i, 1:4] = c(xc[i], yc[i], xc[i], yroot)
        segmentDF$type[nsegs + 1 + i] = "vertical"
        segmentDF$parentNode[nsegs + 1 + i] = root_name
        segmentDF$childNode[nsegs + 1 + i] = child_names[i]
    }
    return(list(pointDF = pointDF, segmentDF = segmentDF))
}

#' Compute offsets for children of a node
#'
#' @param descendant_vec A vector of length equal to the number of
#' children, giving the number of descendants each of the children
#' has.
#'
#' @return A vector of length equal to the number of children, giving
#' the offset for each of the children.
get_child_offsets <- function(descendant_vec) {
    -sum(descendant_vec) / 2 + cumsum(descendant_vec) - descendant_vec/ 2
}


makeChildMapping <- function(tr) {
    map = list()
    for(n in unique(tr$edge[,1])) {
        map[[n]] = tr$edge[tr$edge[,1] == n,2]
    }
    return(map)
}

findRoot <- function(tr) {
    setdiff(tr$edge[,1], tr$edge[,2])
}



tree_theme = theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()
    
)

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

## This function is a modified version of the drop.tip function in ape (https://CRAN.R-project.org/package=ape)
## A line containing drop.singles was deleted to allow for nodes with only one descendant.
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
