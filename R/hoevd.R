#' @title Higher Order Eigen Value Decomposition of A Symmetric Tensor
#' @description Decomposes a higher order symmetric tensors (e.g. higher order centered and
#' standardized tensorian moments) using HOEVD. It is a special case of HOSVD.
#'
#' @param x the symmetric tensorian stored as an array
#' @param rank rank of the decomposition
#'
#' @importFrom rTensor as.tensor k_unfold ttl
#'
#' @return A list with eliments
#' \item{u}{Estimated eigen vectors}
#' \item{z}{Estimated core tensor}
#' @export
#'
#' @details x = Z \code{*} [u,..u]
#'
#' @examples
#' data(dat)
#' moment_center_sd_3 = tensor_moment(Y = dat, center = TRUE, standardize = TRUE, p = 3)
#' decomp_3 = hoevd(x = moment_center_sd_3, rank = 32)
#'


hoevd = function(x,rank = p){
  p = unique(dim(x))
  if(length(p) > 1){
    stop("Input tensor is not symmetric")
  }
  if(rank > p){
    stop("Rank cannot exceed dimension in each mode")
  }

  tnsr = as.tensor(x)
  unfold = k_unfold(tnsr,m = 1)@data
  u = svd(unfold, nu = rank)$u
  ulist = list(u)
  ulist = rep(ulist,tnsr@num_modes)
  z = ttl(tnsr, lapply(ulist, t), ms = 1:tnsr@num_modes)@data
  return(list( u = u, z = z))
}
