# -------------------------------------------------------------
# select a proper fold number for CV
# -------------------------------------------------------------
selfoldforCV <- function(N){
  foldstosel <- 2:9
  Num.quo <- N %/% foldstosel
  Num.rem <- N %% foldstosel
  
  # if no fold meets requirements, then change N
  selIdx <- (Num.rem == 0 & Num.quo >= 10)
  folds <- foldstosel[selIdx]
  folds.num <- length(folds)
  
  if (folds.num != 0) {
    d5 <- abs(rep(5, folds.num) - folds)
    fold <- folds[d5 == min(d5)]
    if (length(fold) > 1) {
      if (N > 300) {
        fold <- max(fold)
      } else {
        fold <- min(fold)
      }
    }
  } else {
    fold <- 1
  }
  
  return(fold)
}


# ------------------
# soft threshold
# ------------------
soft <- function(x, lambda){
  s <- sign(x) * pmax(abs(x) - lambda, 0)
  return(s)
}


# -------------------
# coordinate descent for graphical lasso sub-problem
# -------------------
CDfgL <- function(V, beta_i, u, rho, maxIter = 200, print.detail){
  p_1 <- ncol(V)
  
  # initialize
  beta.new <- rep(0, p_1)
  finished <- rep(FALSE, p_1)
  eps <- 1e-4
  times0 <- 0
  times1 <- 0
  
  iter <- 0
  
  while(TRUE){
    iter <- iter + 1
    beta.old <- beta_i
    for (j in seq_len(p_1)) {
      df <- V %*% beta_i - u
      x <- beta_i[j] - df[j] / V[j, j]
      beta_i[j] <- soft(x, rho / V[j, j])
    }
    beta.new <- beta_i
    
    zeroIdx <- (beta_i == 0)
    if (any(zeroIdx)) {
      if (all(beta.new[zeroIdx] == beta.old[zeroIdx])) {
        times0 <- times0 + 1
      } else {
        times0 <- 0
      }
      if (times0 >= 5) {
        finished[zeroIdx] <- TRUE
      }
    }
    
    if (any(!zeroIdx)) {
      rel_change <- abs((beta.new[!zeroIdx] - beta.old[!zeroIdx]) / beta.old[!zeroIdx])
      if (all(rel_change < eps)) {
        times1 <- times1 + 1
      } else {
        times1 <- 0
      }
      if (times1 >= 3) {
        finished[!zeroIdx] <- TRUE
      }
    }
    
    if (all(finished)) break
    if (iter >= maxIter) {
      if (print.detail) message("CDfgL reached max iteration; breaking.")
      break
    }
  }
  return(beta.new)
}


# -------------------------------------------------------------
# update correction coefficients a and b
# -------------------------------------------------------------
update.CorrectCoef <- function(X0.glist, X1.glist, Theta.list, 
                               a.i, b.i, penal.ksi, penal.gamma, print.detail) {
  G <- length(X0.glist)
  p <- ncol(X0.glist[[1]])
  
  N0_gvec <- rep(0, G)
  N1_gvec <- rep(0, G)
  N_gvec <- rep(0, G)
  for (g in 1:G) {
    N0_gvec[g] <- nrow(X0.glist[[g]])
    N1_gvec[g] <- nrow(X1.glist[[g]])
    N_gvec[g] <- N0_gvec[g] + N1_gvec[g]
  }
  
  a.o <- a.i
  b.o <- b.i
  
  for (j in seq_len(p)) {
    
    # update a[j]
    tmp1.gvec <- rep(0, G)
    tmp2.gvec <- rep(0, G)
    for (g in seq_len(G)) {
      A <- diag(a.o)
      B_gi <- matrix(rep(b.o, each = N1_gvec[g]), nrow = N1_gvec[g])
      
      X1.gi.cor <- X1.glist[[g]] %*% A + B_gi
      X.gi <- rbind(X0.glist[[g]], X1.gi.cor)
      
      X.gi.sca <- scale(X.gi, center = TRUE, scale = TRUE)
      X.gi.sca.attr <- attributes(X.gi.sca)
      Mu_g <- X.gi.sca.attr$`scaled:center`
      Sigma_g <- X.gi.sca.attr$`scaled:scale`
      
      # Ensure Sigma_g[j] is non-zero and not NA
      if (is.na(Sigma_g[j]) || abs(Sigma_g[j]) < 1e-6) {
        Sigma_g[j] <- 1e-6
        if (print.detail) warning(sprintf("Sigma_g[%d] was zero or NA for group
                                          %d in a-update; replaced with 1e-6", j, g))
      }
      
      # Use only the rows corresponding to group 1 in the scaled data
      if ((N0_gvec[g] + 1) > N_gvec[g]) {
        stop("Not enough rows in the combined data for group 1 extraction")
      }
      Y_g <- X.gi.sca[(N0_gvec[g] + 1):N_gvec[g], ]
      
      # Replace the j-th column: adjust b.o[j] by the scaling factors
      if(is.null(nrow(Y_g))) {
        Y_g[j] <- (rep(b.o[j], N1_gvec[g]) - Mu_g[j]) / Sigma_g[j]
      }else {
        Y_g[, j] <- (rep(b.o[j], N1_gvec[g]) - Mu_g[j]) / Sigma_g[j]
      }
      Z_g <- Y_g
      
      tmp <- as.vector(Z_g %*% Theta.list[[g]][, j])
      tmp1.gvec[g] <- 2 / (N_gvec[g] * Sigma_g[j]) * sum(X1.glist[[g]][, j] * tmp, na.rm = TRUE)
      tmp2.gvec[g] <- 2 * Theta.list[[g]][j, j] / (N_gvec[g] * Sigma_g[j]^2) * sum(X1.glist[[g]][, j]^2, na.rm = TRUE)
    }
    denominator_a <- sum(tmp2.gvec, na.rm = TRUE)
    if (abs(denominator_a) < 1e-6) {
      if (print.detail) warning(sprintf("Denominator for updating a[%d] is zero;
                                        keeping previous value", j))
      a.o[j] <- a.o[j]
    } else {
      a.o[j] <- 1 + soft(- sum(tmp1.gvec, na.rm = TRUE) - sum(tmp2.gvec, na.rm = TRUE), penal.ksi) / denominator_a
    }
    
    # update b[j]
    tmp3.gvec <- rep(0, G)
    tmp4.gvec <- rep(0, G)
    for (g in seq_len(G)) {
      A <- diag(a.o)
      B_gi <- matrix(rep(b.o, each = N1_gvec[g]), nrow = N1_gvec[g])
      X1.gi.cor <- X1.glist[[g]] %*% A + B_gi
      X.gi <- rbind(X0.glist[[g]], X1.gi.cor)
      
      X.gi.sca <- scale(X.gi, center = TRUE, scale = TRUE)
      X.gi.sca.attr <- attributes(X.gi.sca)
      Mu_g <- X.gi.sca.attr$`scaled:center`
      Sigma_g <- X.gi.sca.attr$`scaled:scale`
      
      if (is.na(Sigma_g[j]) || abs(Sigma_g[j]) < 1e-6) {
        Sigma_g[j] <- 1e-6
        if (print.detail) warning(sprintf("Sigma_g[%d] was zero or NA for group 
                                          %d in b-update; replaced with 1e-6", j, g))
      }
      Y_g <- X.gi.sca[(N0_gvec[g] + 1):N_gvec[g], ]
      
      # Update j-th column with the current a.o[j] and adjust by scaling
      if(is.null(nrow(Y_g))) {
        Y_g[j] <- (X1.glist[[g]][, j] * a.o[j] - Mu_g[j]) / Sigma_g[j]
      } else {
        Y_g[, j] <- (X1.glist[[g]][, j] * a.o[j] - Mu_g[j]) / Sigma_g[j]
      }
      Z_g <- Y_g
      
      tmp <- as.vector(Z_g %*% Theta.list[[g]][, j])
      tmp3.gvec[g] <- 2 / (N_gvec[g] * Sigma_g[j]) * sum(tmp, na.rm = TRUE)
      tmp4.gvec[g] <- 2 * N1_gvec[g] * Theta.list[[g]][j, j] / (N_gvec[g] * Sigma_g[j]^2)
    }
    denominator_b <- sum(tmp4.gvec, na.rm = TRUE)
    if (abs(denominator_b) < 1e-6) {
      if (print.detail) warning(sprintf("Denominator for updating b[%d] is 
                                        zero; keeping previous value", j))
      b.o[j] <- b.o[j]
    } else {
      b.o[j] <- soft(- sum(tmp3.gvec, na.rm = TRUE), penal.gamma) / denominator_b
    }
  }
  
  coef.out <- list(coef.a = a.o, coef.b = b.o)
  return(coef.out)
}



# -------------------------------------------------------------
# find best parameters for penalty selection via CV+BIC (grid search)
# -------------------------------------------------------------
findBestPara <- function(X0.glist, X1.glist, penal.rho, eps, print.detail) {
  G <- length(X0.glist)
  p <- ncol(X0.glist[[1]])
  
  N0_gvec <- rep(0, G)
  N1_gvec <- rep(0, G)
  N_gvec <- rep(0, G)
  for (g in 1:G) {
    N0_gvec[g] <- nrow(X0.glist[[g]])
    N1_gvec[g] <- nrow(X1.glist[[g]])
    N_gvec[g] <- N0_gvec[g] + N1_gvec[g]
  }
  
  Sel.ksi <- 0
  Sel.gamma <- 0
  
  ksi_candidates <- c(1, 0.5, 0.3, 0.1)
  gamma_candidates <- c(1, 0.5, 0.3, 0.1)
  
  MinAvedist <- Inf
  
  for (ksi in ksi_candidates) {
    for (gamma in gamma_candidates) {
      Allpara <- BEgLasso(X0.glist, X1.glist, penal.rho, ksi, gamma, eps, print.detail)
      X1.cor.glist <- Allpara$X1.cor
      Theta.list <- Allpara$Theta
      
      r <- 0.5
      
      ebic.gvec <- rep(0, G)
      for (i in 1:G) {
        X.gi <- rbind(X0.glist[[i]], X1.cor.glist[[i]])
        X.gi.sca <- scale(X.gi, center = TRUE, scale = TRUE)
        S_i <- cov(X.gi.sca)
        Theta_i <- Theta.list[[i]]
        E.num.gi <- sum(Theta_i[upper.tri(Theta_i, diag = FALSE)] != 0)
        
        # Check determinant of Theta_i before taking its log
        detTheta <- det(Theta_i)
        if (is.na(detTheta) || detTheta <= 0) {
          ebic.gvec[i] <- Inf
        } else {
          ebic.gvec[i] <- - N_gvec[i] * (log(detTheta) - tr(S_i %*% Theta_i)) +
            E.num.gi * log(N_gvec[i]) + 4 * E.num.gi * r * log(p)
        }
      }
      
      ebic <- sum(ebic.gvec)
      
      if (!is.na(ebic) && ebic < MinAvedist) {
        MinAvedist <- ebic
        Sel.ksi <- ksi
        Sel.gamma <- gamma
      }
    }
  }
  
  penterm <- list(penal.ksi = Sel.ksi,
                  penal.gamma = Sel.gamma,
                  MinAvedist = MinAvedist)
  return(penterm)
}



# -------------------------------------------------------------
# CV + BIC select rho
# -------------------------------------------------------------
#' @importFrom stats cov
#' @importFrom lava tr
selrho.useCVBIC <- function(X, print.detail = TRUE) {
  N <- nrow(X)
  fold <- selfoldforCV(N)
  CVset.size <- N / fold
  
  rhos <- seq(from = 0.1, to = 0.9, by = 0.1)
  CVerr1 <- matrix(0, fold, length(rhos))
  CVerr2 <- rep(0, length(rhos))
  
  for (r in seq_along(rhos)) {
    rho <- rhos[r]
    if (fold != 1) {
      for (i in 1:fold) {
        start.index <- (i-1) * CVset.size + 1
        end.index <- i * CVset.size
        X.cv <- X[start.index:end.index, ]
        X.tr <- X[-(start.index:end.index), ]
        
        c.mat <- graphicalLasso(X.tr, rho)
        Theta <- c.mat$Theta
        
        X.cv.sca <- scale(X.cv, center = TRUE, scale = TRUE)
        S.cv <- cov(X.cv.sca)
        
        k <- sum(Theta[upper.tri(Theta, diag = FALSE)] != 0)
        CVerr1[i, r] <- k * log(CVset.size) - CVset.size * (log(det(Theta)) - tr(S.cv %*% Theta))
      }
    } else {
      c.mat <- graphicalLasso(X, rho, print.detail)
      Theta <- c.mat$Theta
      # Added to handle matrix of 1 row
      if (nrow(X) > 1) { 
        X.sca <- scale(X, center = TRUE, scale = TRUE) 
        S <- cov(X.sca) 
      } else { 
        S <- matrix(0, ncol(X), ncol(X)) 
      }
      
      k <- sum(Theta[upper.tri(Theta, diag = FALSE)] != 0)
      CVerr2[r] <- k * log(N) - 2 * (log(det(Theta)) - tr(S %*% Theta))
    }
  }
  
  if (fold != 1) {
    CVerr1 <- colMeans(CVerr1)
    MinCVerr <- min(CVerr1)
    rho.cv <- rhos[CVerr1 == MinCVerr]
  } else {
    MinCVerr <- min(CVerr2)
    rho.cv <- rhos[CVerr2 == MinCVerr]
  }
  
  if (print.detail) {
    message('CVBIC: select rho =', rho.cv, '\n')
  }
  
  return(c(rho.cv, MinCVerr))
}


# -------------------------------------------------------------
# Delete outliers in data 
# -------------------------------------------------------------
#' @importFrom stats sd
#' @importFrom mixOmics pca
DelOutlier <- function(X) {
  pca.dat <- pca(X, ncomp = 3, center = TRUE, scale = TRUE)
  pca.dat.varX <- pca.dat$variates$X
  delsampIdx <- c()
  for (i in 1:3) {
    pc.i <- pca.dat.varX[, i]
    pc.i.m <- mean(pc.i)
    pc.i.sd <- sd(pc.i)
    pc.i.min <- pc.i.m - 3 * pc.i.sd
    pc.i.max <- pc.i.m + 3 * pc.i.sd
    delsampIdx <- c(delsampIdx, which(pc.i < pc.i.min | pc.i > pc.i.max))
  }
  delsampIdx <- unique(delsampIdx)
  
  if (length(delsampIdx) == 0) {
    X.out <- X
  } else {
    X.out <- X[-delsampIdx, ]
  }
  
  Del.result <- list(delsampIdx = delsampIdx, 
                     X.out = X.out)
  return(Del.result)
}


# -------------------------------------------------------------
# Impute outliers in data
# -------------------------------------------------------------
#' @importFrom stats sd
#' @importFrom DMwR2 knnImputation
ImputeOutlier <- function(X) {
  p <- ncol(X)
  X.out <- X
  for (i in 1:p) {
    dat.i <- X[, i]  # Fixed use column i instead of p
    dat.i.m <- mean(dat.i)
    dat.i.sd <- sd(dat.i)
    dat.i.max <- dat.i.m + 4 * dat.i.sd
    dat.i.min <- dat.i.m - 4 * dat.i.sd
    dat.i[dat.i < dat.i.min | dat.i > dat.i.max] <- NA
    X.out[, i] <- dat.i  # assign to column i
  }
  
  na.num <- sum(is.na(X.out))
  if (na.num != 0) {
    X.out <- as.data.frame(X.out)
    X.out <- knnImputation(X.out)
  }
  
  X.out <- as.matrix(X.out)
  
  return(X.out)
}
