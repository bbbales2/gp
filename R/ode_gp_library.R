library(condMVNorm)

p_Xn <- function(tn, Xn, phi_n, sigma_n) {
  
  N <- length(Xn)
  Id <- diag(N)
  
  K_XX <- UU(tn, tn, phi_n)
  K_XsX <- UU(tn, tn, phi_n)
  K_XXs <- t(K_XsX)
  K_XsXs <- UU(tn, tn, phi_n)
  
  m <- rep(0,2*N)
  K <- rbind(cbind(K_XX + sigma_n^2*Id, K_XXs),
             cbind(K_XXs, K_XsXs)) + 1e-6*diag(2*N)
  
  condMVN(m, K, (N+1):(2*N), 1:N, Xn)
}

# given a vector time series data, return mean and cov matrix
# of the time derivative. See page 3 of Accelerating Bayesian
# inference over nonlinear differential equations with Gaussian processes
p_dotXn <- function(tn, Xn, phi_n, sigma_n) {
  
  N <- length(Xn)
  Id <- diag(N)
  
  m <- rep(0,2*N)
  K <- rbind(cbind(UU(tn, tn, phi_n) + sigma_n^2*Id, UD(tn, tn, phi_n)),
             cbind(t(UD(tn, tn, phi_n)), DD(tn, tn, phi_n))) + 1e-6*diag(2*N)
  
  condMVN(m, K, (N+1):(2*N), 1:N, Xn)
}

# simaltaneously do the same for a matrix of time series
p_dotX <- function(X, phi, sigma_sq) {
  
}

# creates function that returns normal mean and variance of the time derivative, dot{x}_n,
# evaluated at a new point x, given all our data X, and all the draws
# of \dot{x} we've already made at previous new (star) points
create_p_dotXnS <- function(Xn_list, mn, Kn, theta) {
  
  X <- do.call(cbind, Xn_list)
  N <- nrow(X)
  D <- ncol(X)
  
  i <- 1
  K_XX <- QQard(X,X,theta)
  K_XsX <- matrix(nrow = 0, ncol = N)
  K_XsXs <- matrix(nrow = 0, ncol = 0)
  
  # pre-factorize inverted matrices
  K_XX_qr <- qr(K_XX + 1e-6*diag(N))
  K_XX_1_mn <- solve(K_XX_qr, mn)
  K_XX_1_Kn <- solve(K_XX_qr,Kn)
  
  # input points we already looked at (x) and corresponding
  # derivative points we already drew
  Xs <- matrix(nrow = 0, ncol = D)
  dot_Xs <- matrix(nrow = 0, ncol = D)
  
  #xs must be 1 X D matrix
  p_dotXnS <- function(xs_vec) {
    
    # new input state to draw derivative from
    xs <- matrix(xs_vec, nrow = 1)
    
    # update gram matrix to include new point
    K_XsX <<- rbind(K_XsX, QQard(xs,X,theta))
    K_XsXs <<- rbind(cbind(K_XsXs, QQard(Xs,xs,theta)),
                     cbind(t(QQard(Xs,xs,theta)), QQard(xs,xs,theta)))
    
    # reset new joint mean and covariance of all star points
    m <- K_XsX %*% K_XX_1_mn
    K <- K_XsXs - K_XsX %*% solve(K_XX_qr, t(K_XsX)) + K_XsX %*% K_XX_1_Kn %*% solve(K_XX_qr, t(K_XsX))
    K <- (K + t(K))/2 + 1e-6*diag(dim(K)[1])
    
    # get conditional dist
    if(i == 1) cond_dist <- condMVN(m, K, i)
    else cond_dist <- condMVN(m, K, i, 1:(i-1), c(dot_Xs))
    
    dot_xs <- rnorm(1,as.numeric(cond_dist$condMean), as.numeric(cond_dist$condVar))
    
    # update count and list of drawn locations
    i <<- i + 1
    Xs <<- rbind(Xs, xs)
    dot_Xs <<- rbind(dot_Xs, matrix(dot_xs, nrow = 1))
    
    return(list(mu = as.numeric(cond_dist$condMean), sigma = as.numeric(cond_dist$condVar), dot_xs = dot_xs))
  }
}