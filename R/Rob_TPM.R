#' @title Robust Power Tensor Method for Decomposition of A Symmetric Moment Tenosr
#' @description Decomposes a higher order symmetric tensors (e.g. higher order centered and
#' standardized tensorian moments) using Robut Power Tensor Method (RPTM) based on the vectorization
#' provided by Di et al.
#'
#' @param Y the symmetric tensorian stored as an array
#' @param center column-center the data first, default is \code{TRUE}.
#' @param standardize standardize the multivariate data, i.e. convert it to the left singular matrix
#' @param N number of inner loop, i.e. steps for each power tensor iteration to converge
#' @param L number of outloop, i.e. number of power tensor iterations
#' @param rank rank of the decomposition
#' @param order order of the decomposition, 3 or 4
#' @param tol tolerance for the inner loop, i.e power tensor iteration, default to be 1e-06
#'
#' @importFrom stats rnorm
#'
#' @importFrom rTensor as.tensor k_unfold ttl
#'
#' @return A list with eliments
#' \item{eigenv}{Estimated eigen vectors, sorted by based on eigenl}
#' \item{eigenl}{Estimated eigen values, sorted from largest to smallest}
#' \item{w}{The whitening matrix w if \code{standardize = TRUE}}
#' @export
#'
#' @details See Di et al. 2018. This function does not contruct tensor first. It can only be used directly to the
#' data matrix, and it produces decomposition of the (standardized) moment tensors
#'
#' @references Di et al.
#' @references A Anandkumar et al. Tensor decompositions for learning latent variable models, 2012.
#'
#' @examples
#' data(dat)
#' rob_3 = Rob_TPM(Y = dat, order = 3, rank = 10, tol = 1e-6)
#'


Rob_TPM = function(Y, center = TRUE, standardize = TRUE, L = 10, N = 10, order = c(3,4), rank = ncol(Y), tol = 1e-6){
  Y = as.matrix(Y)

  if(!order %in% c(3,4)){
    stop("order has to be either 3 or 4")
  }
  if(center){
    Y = scale(Y, center = TRUE, scale = FALSE)
  }

  if(standardize){
    # u = svd(Y)$u
    v = svd(Y)$v
    s = diag(svd(Y)$d)
    w = v %*% solve(s)
    # Y = Y %*% w
    Y = svd(Y)$u
    rm(list = c("v","s"))
  }


  result.1 = Est_PI1(Y = Y, L = L, N = N, order = order, tol = tol)
  eigenv = matrix(result.1$theta_hat,nrow = length(result.1$theta_hat))
  eigenl = result.1$lambda_hat
  print("Eigen pairs number: 1")

  for( m in 2:rank){
    result.m = Est_PIk(Y = Y,eigenv = eigenv, eigenl = eigenl, L = L, N = N, order = order, tol = tol)
    eigenv = cbind(eigenv, result.m$theta_hat)
    eigenl = c(eigenl, result.m$lambda_hat)
    print(paste("Eigen pairs number:",m, sep = " "))
  }

  eigenv = eigenv[,order(abs(eigenl),decreasing = T)]
  eigenl = eigenl[order(abs(eigenl),decreasing = T)]

  if(!standardize){
    w = NA
  }

  result = list("eigenv" = eigenv, "eigenl" = eigenl,"w" = w)
}

# power iteration for the first eigen vecotr
power_itr1 = function(theta, Y, N, order, tol){
  n = nrow(Y)
  d = 1/n

  for(k in 1:N){
    w = (Y %*% theta)^(order-1)
    next_itr = as.vector(t(w) %*% Y)
    next_itr = d * next_itr

    thetai = next_itr/sqrt(sum(next_itr^2))
    err = sum((thetai - theta)^2)
    theta = thetai
    if( err <= tol){break}
  }
  return(theta)
}

# power iteration for all the following eigen vectors starting from 2nd
power_itrk = function(theta, eigenv, eigenl, Y, N, order, tol){
  n = nrow(Y)
  d = 1/n

  for(k in 1:N){
    w1 = (Y %*% theta)^(order-1)
    next_itr1 = as.vector(t(w1) %*% Y)

    w2 = (t(eigenv) %*% theta)^(order - 1)
    if(length(eigenl) == 1){next_itr2 = as.vector(eigenv %*% eigenl %*% w2)}
    if(length(eigenl) >1){next_itr2 = eigenv %*% diag(eigenl) %*% w2}

    next_itr = d * next_itr1 - next_itr2
    thetai = next_itr/sqrt(sum(next_itr^2))
    err = sum((thetai - theta)^2)
    theta = thetai
    if( err <= tol ){break}
  }
  return(theta)
}

# get first lambda
Est_PI1 = function(Y, L, N, order, tol){
  p = ncol(Y)
  theta_list = list()
  for ( t in 1:L){
    vc = rnorm(p)
    theta_0 = vc/sqrt(sum(vc^2))
    theta_list[[t]] = power_itr1(theta = theta_0, Y = Y, N=N, order = order, tol = tol)
  }

  n = nrow(Y)
  d = 1/n
  lambda_list = NULL

  for(t in 1:L){
    a = (Y %*% theta_list[[t]])^order
    lambda_t = d * sum(a)
    lambda_list = c(lambda_list,lambda_t)
  }

  ind = which.max(abs(lambda_list))
  theta_tau = theta_list[[ind]]

  theta_hat = power_itr1(theta = theta_tau,Y = Y, N = N, order = order, tol = tol)

  a = (Y %*% theta_hat)^order
  lambda_hat = d * sum(a)

  result = list("theta_hat" = theta_hat, "lambda_hat" = lambda_hat)
  return(result)
}

# get all following lambda
Est_PIk = function(Y, eigenv, eigenl, L=10, N=10, order, tol){
  p = ncol(Y)
  theta_list = list()
  for ( t in 1:L){
    vc = rnorm(p)
    theta_0 = vc/sqrt(sum(vc^2))
    theta_list[[t]] = power_itrk(theta = theta_0, eigenv, eigenl, Y = Y, N=N, order = order, tol = tol)
  }

  n = nrow(Y)
  d = 1/n
  lambda_list = NULL

  for(t in 1:L){
    a =  sum((Y %*% theta_list[[t]])^order)
    if(length(eigenl) == 1){ b =  sum(eigenl * ((t(eigenv) %*% theta_list[[t]])^order))}
    if(length(eigenl) >1){b =  sum(diag(eigenl) %*% ((t(eigenv) %*% theta_list[[t]])^order))}
    lambda_t = d * a - b
    lambda_list = c(lambda_list,lambda_t)
  }

  ind = which.max(abs(lambda_list))
  theta_tau = theta_list[[ind]]
  theta_hat = power_itrk(theta = theta_tau,eigenv, eigenl,Y = Y, N = N, order = order, tol = tol)

  a =  sum((Y %*% theta_hat)^order)
  if(length(eigenl) == 1){ b =  sum(eigenl * ((t(eigenv) %*% theta_hat)^order))}
  if(length(eigenl) >1){b =  sum(diag(eigenl) %*% ((t(eigenv) %*% theta_hat)^order))}

  lambda_hat = d * a - b

  result = list("theta_hat" = theta_hat, "lambda_hat" = lambda_hat)
  return(result)
}
