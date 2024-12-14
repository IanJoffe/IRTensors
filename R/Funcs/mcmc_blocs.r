#### Starting values
set.seed(seed)
K <- length(dim(X)) - 1
C <- 3  # Number of blocs, adjust as needed

# Initialize bloc-specific variables
B_blocs <- vector("list", C)
S_blocs <- vector("list", C)
s2_blocs <- rep(0, C)
z <- sample(1:C, dim(Y)[1], replace = TRUE)  # Random bloc assignments

if (!rstart) {
  S <- list()
  
  for (c in 1:C) {
    bloc_indices <- which(z == c)
    X_c <- X[bloc_indices, bloc_indices, , , drop = FALSE]
    Y_c <- Y[bloc_indices, , , , drop = FALSE]
    B_blocs[[c]]<-BMLE<-mlm.ALS(Y_c,X_c,imax=5,verbose=TRUE) 
    YI<-Y_c ; YI[is.na(Y_c)]<-tprod(X_c,B_blocs[[c]])[is.na(Y_c)]
    R<-YI-tprod(X_c,B_blocs[[c]])
    
    S_blocs[[c]]<-list() 
    for(k in 1:K)
    {
      S_blocs[[c]][[k]]<-tcrossprod(mat(R,k))
      S_blocs[[c]][[k]]<-S_blocs[[c]][[k]]*nrow(S_blocs[[c]][[k]])/tr(S_blocs[[c]][[k]])
    }
    s2_blocs[[c]]<-mean(R^2)
  }
} 
if (rstart) {
  # The authors don't use rstarts so I haven't bothered to change this
  for (c in 1:C) {
    B_blocs[[c]] <- lapply(1:K, function(k) rsan(c(dim(Y)[k], dim(X)[k])))
    S_blocs[[c]] <- lapply(1:K, function(k) solve(rwish(diag(dim(Y)[k]))))
    s2_blocs[c] <- rgamma(1, 1, 1)
  }
  YI <- Y
  YI[is.na(Y)] <- tprod(X, B_blocs[[z]])
}

#### Output
MSE <- vector("list", C)
S2 <- NULL
BPS <- SPS <- vector("list", K)
for (k in 1:K) {
  BPS[[k]] <- array(dim = c(dim(Y)[k], dim(X)[k], 0))
  SPS[[k]] <- array(dim = c(dim(Y)[k], dim(Y)[k], 0))
}

#### MCMC
for (s in 1:(NS + NB)) {
  ## Reassign countries to blocs
  for (i in 1:nrow(Y)) {
    cat("Iteration ", s, " Reassigning Country ", i, "\n")
    current_bloc <- z[i]
    best_loss <- Inf
    best_blocs <- B_blocs
    
    for (c in 1:C) {
      tse_c <- 0
      z[i] <- c  # Temporarily assign to bloc c
      potential_b_blocs <- vector("list", C)
      
      for (c_inner in 1:C) {
        potential_b_blocs[[c]][[c_inner]] <- vector("list", 3)
        X_c <- X[which(z == c), which(z == c), , , drop = FALSE]
        Y_c <- Y[which(z == c), which(z == c), , , drop = FALSE]
        B_c<-BMLE<-mlm.ALS(Y_c,X_c,imax=5,verbose=FALSE)
        potential_b_blocs[[c]][[c_inner]] <- B_c
        XB_bloc <- tprod(X_c, B_c)
        tse_c <- tse_c + sum((Y_c - XB_bloc)^2)
      }
      if (tse_c < best_loss) {
        best_assignment <- c
        best_loss <- tse_c
        best_blocs <- potential_b_blocs
      }
    }
    z[i] <- best_assignment  # Assign to the best bloc
    B_blocs <- best_blocs
  }
  
  ## Update bloc-specific parameters
  for (c in 1:C) {
    bloc_indices <- which(z == c)
    if (length(bloc_indices) > 0) {
      X_c <- X[bloc_indices, , , , drop = FALSE]
      Y_c <- Y[bloc_indices, , , , drop = FALSE]
      B_c <- mlm.ALS(Y_c,X_c,imax=5,verbose=TRUE)
      BS <- rBSa_fc(Y_c, X_c, B_c, S_blocs[[c]])
      B_blocs[[c]] <- BS$B
      S_blocs[[c]] <- BS$S
      s2_blocs[c] <- BS$s2
    }
  }
  
  ## Impute missing data
  for (c in 1:C) {
    bloc_indices <- which(z == c)
    if (length(bloc_indices) > 0) {
      XB_bloc <- tprod(X[bloc_indices, , , , drop = FALSE], B_blocs[[c]])
      Y[block_indices, , , , drop = FALSE][is.na(Y[block_indices, , , , drop = FALSE])] <- 
        (XB_bloc + sqrt(s2_blocs[c]) * tprod(rsan(dim(Y)[bloc_indices, , , drop = FALSE]), 
                                             lapply(S_blocs[[c]], mhalf)))[is.na(Y[block_indices, , , , drop = FALSE])]
    }
  }
  
  ## Output
  for (c in 1:C) {
    if (length(which(z == c)) > 0) {
      XB_bloc <- tprod(X[which(z == c), , , drop = FALSE], B_blocs[[c]])
      mse_c <- mean((YI[which(z == c), , , , drop = FALSE] - XB_bloc)^2)
      MSE[[c]] <- c(MSE[[c]], mse_c)
      cat(MSE[[c]][length(MSE[[c]])])
    }
  }
  
  if (plot) {
    par(mfrow = c(2, 2))
    lapply(1:C, function(c) plot(MSE[[c]], type = "l", main = paste("Bloc", c)))
  }
  
  if (s %% sdens == 0) {
    save.image(file = fname)
  }
}
save.image(file = fname)