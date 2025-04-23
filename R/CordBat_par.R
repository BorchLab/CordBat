CordBat_par <- function(X, batch, group = NULL, grouping = F, ref.batch, 
                            eps = 1e-5, print.detail = T, n_cores = 1, skip.impute = FALSE) {
  # X is the log-transformed data matrix of all batches
  # batch is the batch IDs for samples, a vector
  # group is the group IDs for samples, a vector
  # ref.batch is the ID of the reference batch
  if (is.null(group)) {
    containQC <- FALSE
    group <- rep(1, nrow(X))
  }
  else {
    if (sum(group == 'QC') != 0) {
      containQC <- T
      
      X.init <- X
      batch.init <- batch
      group.init <- group
      
      X_QC <- X[group == 'QC', ]
      X <- X[group != 'QC', ]
      
      QC_batch <- batch[group == 'QC']
      batch <- batch[group != 'QC']
      group <- group[group != 'QC']
    }
    else {
      containQC <- F
      X_QC <- NULL # AA
      batch.init <- NULL
      group.init <- NULL
      QC_batch <- NULL
      
    }
    
    if (!grouping) {
      group <- rep(1, nrow(X))
    }
  }
  
  X <- as.matrix(X)
  p <- ncol(X)
  
  # Force batch and group to be character vectors for comparisons
  batch <- as.character(batch)
  group <- as.character(group)
  
  batch.f <- factor(batch)
  batch.levels <- levels(batch.f)
  batch.num <- length(batch.levels)
  
  group.f <- factor(group)
  group.levels <- levels(group.f)
  group.num <- length(group.levels)
  
  delsampIdx <- integer()
  if(!skip.impute) {
    # Outlier Removal and Imputation
    X.delout <- X
    for (i in seq_len(batch.num)) {
      cur.batch <- batch.levels[i]
      bati.idx <- which(batch == cur.batch)
      if (length(bati.idx) == 1) next  # Skip batches with only 1 sample
      X.bati <- DelOutlier(X[bati.idx, , drop = FALSE])
      delsamp.bati <- X.bati$delsampIdx
      
      if (length(delsamp.bati) > 0) {
        dat.bati <- ImputeOutlier(X.bati$X.out)  # Perform imputation only if deletion happened
        bati.delinitIdx <- bati.idx[delsamp.bati]
        delsampIdx <- c(delsampIdx, bati.delinitIdx)
        X.delout[bati.idx[-delsamp.bati], ] <- dat.bati
      } else {
        X.delout[bati.idx, ] <- X.bati$X.out  # No deletion, use original data
      }
    }
    if (length(delsampIdx) > 0) {
      X.nodel <- X.delout
    } else {
      X.nodel <- X
    }
    # Refresh factors after outlier removal
    batch.f <- factor(batch)
    batch.levels <- levels(batch.f)
    batch.num <- length(batch.levels)
    group.f <- factor(group)
    group.levels <- levels(group.f)
    group.num <- length(group.levels)
  } else {
    X.nodel <- X # Skip outlier detection/imputation
    X.delout <- X
  }
  
  # Initialize correction parameters
  Theta.list <- vector("list", group.num)
  for (i in seq_len(group.num)) Theta.list[[i]] <- matrix(0, p, p)
  a <- rep(1, p)
  b <- rep(0, p)
  para <- list(Theta = Theta.list, coef.a = a, coef.b = b)
  
  Xcor.para <- vector("list", batch.num)
  for (i in seq_len(batch.num)) Xcor.para[[i]] <- para
  
  # Initialize corrected matrices 
  #n_del <- nrow(X.delout)
  n_del <- nrow(X)
  X.cor   <- matrix(0, n_del, p)
  X.cor.1 <- matrix(0, n_del, p)
  
  ref.batch_char <- as.character(ref.batch)
  ref.idx <- which(batch == ref.batch_char)
  X.cor[ref.idx, ]   <- X.delout[ref.idx, ]
  X.cor.1[ref.idx, ] <- X.delout[ref.idx, ]
  
  X.cor.withQC <- NULL
  if (containQC) {
    X.cor.withQC <- matrix(0, nrow(X.init), p)
    ref.idx.init <- which(as.character(batch.init) == as.character(ref.batch))
    X.cor.withQC[ref.idx.init, ] <- X.init[ref.idx.init, ]
  }
  
  #### community detection ------------------------
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  Xb0.mat <- X.delout[batch.new == ref.batch, ]
  COM <- getAllCom(Xb0.mat)
  cat("Community detection: ", length(COM), "communities", "\n",
      "Size: ", lengths(COM), "\n")
  
  for (i in seq_along(COM)) {
    metID <- COM[[i]]
    Xb0.COMi <- Xb0.mat[, metID]
    Nb0.COMi <- nrow(Xb0.COMi)
    
    Xb0.COMi.glist <- list()
    Nb0.COMi.gvec <- rep(0, group.num)
    for (g in 1:group.num) {
      Xb0.gi.COMi <- X.delout[batch.new == ref.batch & group.new == grp.level[g], metID]
      Xb0.COMi.glist[[g]] <- Xb0.gi.COMi
      Nb0.COMi.gvec[g] <- nrow(Xb0.gi.COMi)
    }
    
    # Select rho for each group
    rhos <- numeric(group.num)
    for (g in 1:group.num) {
      if (length(metID) > 5) {
        rho_g <- StARS(Xb0.COMi.glist[[g]], round(0.7 * Nb0.COMi.gvec[g]), 100, print.detail = print.detail)[1]
      } else {
        rho_g <- selrho.useCVBIC(Xb0.COMi.glist[[g]], print.detail = print.detail)[1]
      }
      rhos[g] <- rho_g
    }
    rho <- mean(rhos)
    if (print.detail) cat('Set rho = ', rho, '\n')
    
    # Parallel batch correction for this community
    batch_results <- foreach(k = 1:batch.num, .packages = c(), .combine = 'c', 
                             .export = c("X.delout", "X.nodel", "X_QC", "findBestPara", "BEgLasso", "ref.batch", 
                                         "rho", "eps", "group.num", "grp.level", "metID", "batch.new", 
                                         "batch.level", "group.new", "batch", "batch.init", "group.init", "QC_batch", 
                                         "print.detail", "containQC", "tr","CDfgL","soft","update.CorrectCoef")) %dopar% {
                                           if (batch.level[k] == ref.batch) return(NULL)
                                           
                                           Xb1.Batk.COMi.glist <- list()
                                           for (g in 1:group.num) {
                                             Xb1.gi.Batk <- X.delout[batch.new == batch.level[k] & group.new == grp.level[g], ]
                                             Xb1.Batk.COMi.glist[[g]] <- Xb1.gi.Batk[, metID]
                                           }
                                           
                                           penterm <- findBestPara(Xb0.COMi.glist, Xb1.Batk.COMi.glist, rho, eps)
                                           para.out <- BEgLasso(Xb0.COMi.glist, Xb1.Batk.COMi.glist, rho, penterm$penal.ksi, penterm$penal.gamma, eps)
                                           
                                           batch_data <- list(
                                             k = k,
                                             X1.cor.glist = para.out$X1.cor,
                                             Theta.list = para.out$Theta,
                                             coef.a = para.out$coef.a,
                                             coef.b = para.out$coef.b,
                                             metID = metID
                                           )
                                           return(list(batch_data))
                                         }
    
    # Integrate results
    for (res in batch_results) {
      if (is.null(res)) next
      k <- res$k
      for (g in 1:group.num) {
        X.cor.1[batch.new == batch.level[k] & group.new == grp.level[g], res$metID] <- res$X1.cor.glist[[g]]
        Xcor.para[[k]]$Theta[[g]][res$metID, res$metID] <- res$Theta.list[[g]]
      }
      Xcor.para[[k]]$coef.a[res$metID] <- res$coef.a
      Xcor.para[[k]]$coef.b[res$metID] <- res$coef.b
      
      Xb1.nodel <- X.nodel[batch == batch.level[k], ]
      N1 <- nrow(Xb1.nodel)
      coef.A <- diag(res$coef.a)
      coef.B <- matrix(1, N1, 1) %*% t(res$coef.b)
      X.cor[batch == batch.level[k], res$metID] <- Xb1.nodel[, res$metID] %*% coef.A + coef.B
      
      if (!is.na(containQC) && containQC) {
        X.cor.withQC[batch.init == batch.level[k] & group.init != 'QC', res$metID] <- Xb1.nodel[, res$metID] %*% coef.A + coef.B
        QC.batk <- X_QC[QC_batch == batch.level[k], ]
        Nqc.batk <- nrow(QC.batk)
        coef.B <- matrix(1, Nqc.batk, 1) %*% t(res$coef.b)
        X.cor.withQC[batch.init == batch.level[k] & group.init == 'QC', res$metID] <- QC.batk[, res$metID] %*% coef.A + coef.B
      }
    }
    cat('Finish correction of community ', i, '\n')
  }
  
  stopCluster(cl)
  
  return(list(
    batch.level = batch.level,
    delsampIdx = delsampIdx,
    batch.new = batch.new,
    group.new = group.new,
    X.delout = X.delout,
    X.cor = X.cor,
    X.cor.1 = X.cor.1,
    X.cor.withQC = X.cor.withQC,
    Xcor.para = Xcor.para
  ))
}