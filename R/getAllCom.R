#' Community Detection Based on Feature Correlations
#'
#' This function identifies communities (clusters) within features using correlation-based clustering.
#' Features that are highly correlated are grouped together to form communities.
#'
#' @param X A numeric matrix where rows represent samples and columns represent features.
#'
#' @return A list where each element is a vector of feature indices belonging to a detected community.
#'
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(100), 10, 10)
#' communities <- getAllCom(X)
#' str(communities)
#'
#'@importFrom stats cor
#'@export
getAllCom <- function(X) {
  N <- nrow(X)  # Number of samples
  p <- ncol(X)  # Number of features
  
  # Compute Pearson correlation coefficient matrix
  G <- cor(X)
  G <- (G + t(G)) / 2  # Ensure symmetry
  
  # Initialize community structure
  lastCOM <- list(seq_len(p))
  maxCOMSize <- min(round(N * 0.4), 50)
  Nochange.times <- 0
  
  # Iteratively refine communities
  repeat {
    com.num.old <- length(lastCOM)
    lastCOM <- ComtyDet(G, lastCOM, maxCOMSize)
    com.num.new <- length(lastCOM)
    
    cmty.sizemax <- max(lengths(lastCOM))
    
    if (com.num.new == com.num.old) {
      Nochange.times <- Nochange.times + 1
    } else {
      Nochange.times <- 0
    }
    
    if (Nochange.times == 3 || cmty.sizemax <= maxCOMSize) break
  }
  
  # Assign feature indices to community IDs
  NumComty <- length(lastCOM)
  ComtyID <- integer(p)
  for (i in seq_len(NumComty)) {
    ComtyID[lastCOM[[i]]] <- i
  }
  
  # Merge single-feature communities into nearest correlated community
  cmty.sizeis1 <- which(lengths(lastCOM) == 1)
  if (length(cmty.sizeis1) > 0) {
    for (i in cmty.sizeis1) {
      S1.metID <- lastCOM[[i]]
      allCorwithS1 <- G[, S1.metID]
      allCorwithS1[S1.metID] <- 0  # Ignore self-correlation
      Maxcor.metID <- which.max(allCorwithS1)
      
      mergeCmty <- ComtyID[Maxcor.metID]
      lastCOM[[mergeCmty]] <- c(lastCOM[[mergeCmty]], S1.metID)
      ComtyID[S1.metID] <- mergeCmty
    }
    lastCOM <- lastCOM[-cmty.sizeis1]
  }
  
  # Merge two-feature communities into nearest correlated community
  cmty.sizeis2 <- which(lengths(lastCOM) == 2)
  if (length(cmty.sizeis2) > 0) {
    for (i in cmty.sizeis2) {
      S2.metID <- lastCOM[[i]]
      allCorwithS2 <- G[, S2.metID]
      allCorwithS2[S2.metID, ] <- 0  # Ignore self-correlation
      Maxcor.metID <- which(allCorwithS2 == max(allCorwithS2), arr.ind = TRUE)[1]
      
      mergeCmty <- ComtyID[Maxcor.metID]
      lastCOM[[mergeCmty]] <- c(lastCOM[[mergeCmty]], S2.metID)
    }
    lastCOM <- lastCOM[-cmty.sizeis2]
  }
  
  return(lastCOM)
}

# ----------------------------------------------
# Community detection
# ----------------------------------------------
#' @importFrom igraph graph_from_adjacency_matrix cluster_louvain membership 
#' sizes
ComtyDet <- function(G, InputCOM, minCOMSize){

  nextCOM <- list()
  k <- 0
  
  InputCOM.size <- lengths(InputCOM)
  ii <- (InputCOM.size <= minCOMSize)
  if (sum(ii) != 0) {
    k <- sum(ii)
    nextCOM <- InputCOM[ii]
  }
  
  # get communities whose sizes are larger than minCOMSize
  InputCOM <- InputCOM[!ii]
  
  # Community detection
  for (i in c(1: length(InputCOM))) {
    InputG <- abs(G[InputCOM[[i]], InputCOM[[i]]])
    InputG <- InputG - diag(diag(InputG))
    Graph.InputG <- graph_from_adjacency_matrix(InputG, 
                                                weighted = TRUE, 
                                                mode = "undirected")
    COMTY <- cluster_louvain(Graph.InputG)
    cmty.id <- membership(COMTY)
    cmty.size <- as.vector(sizes(COMTY))
    NumComty <- length(cmty.size)
    
    for (n in c(1: NumComty)) {
      Node.Idx <- (cmty.id == n)
      Node.InitID <- InputCOM[[i]][Node.Idx]
      k <- k + 1
      nextCOM[[k]] <- Node.InitID
    }
  }
  
  return(nextCOM)
}

