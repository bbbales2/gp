#convert and NxD matrix in to a list with N vectors of length D
mat_to_obs_list <- function(X) lapply(1:nrow(X), function(i) X[i,])

#apply the function K(x,y), which takes in two vectors, to every combination of two lists X,Y
obs_list_outer <- function(X,Y,K) outer(X, Y, function(X,Y) vapply(seq_along(X), function(i) K(X[[i]], Y[[i]]), numeric(1))) 

#given a function K(x,y,phi) that takes in two vectors and a list
#of parameters, create version that can take in two matrices X,Y
#and return a matrix with the function applied to all combinations
#of a row vectors from X and a row vector from Y
create_kernel_function <- function(K) {
  function(X,Y,phi) {
    obs_list_outer(mat_to_obs_list(X),
                   mat_to_obs_list(Y),
                   function(x,y) K(x,y,phi))
  }
}

QQard <- create_kernel_function(function(x,y,phi) phi[[1]]^2*exp(-(1/2)*sum(((x-y)/phi[[2]])^2)))

#for 1-D input we use these kernels
QQ <- function(x,y,phi) {
  outer(x,y, function(x,y) phi[[1]]^2*exp(-((x - y)^2/(2 * phi[[2]]^2))))
}

QR <- function(x,y,phi) {
  outer(x,y, function(x,y) phi[[1]]^2*(exp(-((x - y)^2/(2 * phi[[2]]^2))) * (x - y))/phi[[2]]^2)
}

RR <- function(x,y,phi) {
  outer(x,y, function(x,y) phi[[1]]^2*(exp(-((x - y)^2/(2 * phi[[2]]^2))))/phi[[2]]^2 - (exp(-((x - y)^2/(2 * phi[[2]]^2))) * (x - y)^2)/phi[[2]]^4)
}
