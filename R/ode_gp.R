p_Xn <- function(tn, Xn, phi_n, sigma_n) {
  
  N <- length(Xn)
  Id <- diag(N)
  K_XX <- QQ(tn, tn, phi_n) + sigma_n^2*Id
  K_XsX <- QQ(tn, tn, phi_n)
  K_XXs <- t(K_XsX)
  K_XsXs <- QQ(tn, tn, phi_n)
  
  mn <- K_XsX %*% solve(K_XX, Xn)
  Kn <- K_XsX - K_XsX %*% solve(K_XX,K_XXs)
  
  return(list(mn = mn, Kn = Kn))
}

# given a vector time series data, return mean and cov matrix
# of the time derivative. See page 3 of Accelerating Bayesian
# inference over nonlinear differential equations with Gaussian processes
p_dotXn <- function(tn, Xn, phi_n, sigma_n) {
  
  N <- length(Xn)
  Id <- diag(N)
  QQ <- QQ(tn, tn, phi_n)
  QR <- QR(tn, tn, phi_n)
  RQ <- t(QR)
  RR <- RR(tn, tn, phi_n)
  
  mn <- RQ %*% solve(QQ + sigma_n^2*Id, Xn)
  Kn <- RR - RQ %*% solve(QQ + sigma_n^2*Id,QR)
  
  return(list(mn = mn, Kn = Kn))
}

# simaltaneously do the same for a matrix of time series
p_dotX <- function(X, phi, sigma_sq) {
  
}

# creates function that returns normal mean and variance of the time derivative, dot{x}_n,
# evaluated at a new point x, given all our data X, and all the draws
# of \dot{x} we've already made at previous new (star) points
create_p_dotXnS <- function(Xn_list, mn, Kn, theta) {
  
  X <- cbind(Xn_list)
  N <- nrow(X)
  D <- ncol(X)
  
  i <- 0
  K_XX <- QQard(X,X,theta)
  K_Xs1X <- matrix(nrow = 0, ncol = N)
  K_XXs1 <- t(K_Xs1X)
  K_Xs1Xs1 <- matrix(nrow = 0, ncol = 0)
  x1 <- c()
  
  K_XX_qr <- qr.(K_XX, LAPACK = TRUE)
  K_XX_1_mn <- solve(K_XX_qr, mn)
  
  K_XX_Kn <- solve(K_XX_qr,Kn)
  
  p_dotXnS <- function(xs) {
    
    K_Xs2X <- QQard(xs,X,theta)
    K_XXs2 <- t(K_Xs2X)
    
    K_Xs1Xs2 <- QQard(x1,xs,theta)
    K_Xs2Xs1 <- t(K_Xs1Xs2)
    K_Xs2Xs2 <- QQard(xs,xs,theta)
    
    #if this is the first draw there's no K_Xs1X
    if(i == 0) {
      mu_s2 <- K_Xs2X %*% K_XX_1_mn
      
      K_XX_1_K_XXs2 <- solve(K_XX_qr, K_XXs2)
      Sigma <- QQard(xs,xs) - K_Xs2X %*% K_XX_1_K_XXs2
      sigma_s2 <- Sigma + K_Xs2X %*% K_XX_Kn %*% K_XX_1_K_XXs2
    }
    else {
      mu_1 <- K_Xs1X %*% K_XX_1_mn
      mu_2 <- K_Xs2X %*% K_XX_1_mn
      
      K_XX_1_K_XXs1 <- solve(K_XX_qr, K_XXs1)
      K_XX_1_K_XXs2 <- solve(K_XX_qr, K_XXs2)
      
      Sigma_11 <- K_Xs1Xs1 - K_Xs1X %*% solve(K_XX_qr, K_XXs1) + K_Xs1X %*% K_XX_Kn %*% K_XX_1_K_XXs1
      Sigma_12 <- K_Xs1Xs2 - K_Xs1X %*% solve(K_XX_qr, K_XXs2) + K_Xs1X %*% K_XX_Kn %*% K_XX_1_K_XXs2
      Sigma_21 <- K_Xs2Xs1 - K_Xs2X %*% solve(K_XX_qr, K_XXs1) + K_Xs2X %*% K_XX_Kn %*% K_XX_1_K_XXs1
      Sigma_22 <- K_Xs2Xs2 - K_Xs2X %*% solve(K_XX_qr, K_XXs2) + K_Xs2X %*% K_XX_Kn %*% K_XX_1_K_XXs2
      
      mu_s2 <- mu_2 + Sigma_21 %*% solve(Sigma_11, x1 - mu_1)
      sigma_s2 <- Sigma_22 - Sigma_21 %*% solve(Sigma_11, Sigma_12)
    }
    
    dotxs <- rnorm(mu_s2, sigma_s2)
    
    i <<- i + 1
    K_Xs1X <<- rbind(K_Xs1X, K_Xs2X)
    K_XXs1 <<- t(K_Xs1X)
    K_Xs1Xs1 <<- cbind(rbind(K_Xs1Xs1,K_Xs1Xs2), rbind(K_Xs2Xs1,K_Xs2Xs2))
    x1 <<- c(x1, dotxs)
    
    return(list(mu = mu_s2, sigma = sigma_s2, dotxs = dotxs))
  }
}

