#' @title Construct Higher Order Tensorian Moments
#' @description 3rd and 4th order (centered and standardized) higher order tensorian moments
#' @param Y data.frame or matrix of multivariate data of dimension n*p
#' @param center column-center the data first, default is \code{TRUE}.
#' @param standardize standardize the multivariate data, i.e. convert it to the left singular matrix
#' @param p what is the order to construct
#'
#'
#' @return 3rd or 4th order moments
#'
#' @export
#'
#' @details Suppose Y = USV' is the SVD, standardization is done by taking U only
#'
#' @examples
#' data(dat)
#' moment_center_sd_3 = tensor_moment(Y = dat, center = TRUE, standardize = TRUE, p = 3)


tensor_moment = function (Y, center = TRUE, standardize = TRUE, p = c(3,4)) {
  Y = as.matrix(Y)

  if(!p %in% c(3,4)){
    stop("p has to be either 3 or 4")
  }
  if(center){
    Y = scale(Y, center = TRUE, scale = FALSE)
  }

  if(standardize){
    # u = svd(Y)$u
    # v = svd(Y)$v
    # s = diag(svd(Y)$d)
    # w = v %*% solve(s)
    # Y = Y %*% w
    Y = svd(Y)$u
  }

  n = nrow (Y)
  m = ncol (Y)
  mm = 1 : m
  nn = 1 : n
  if(p == 3){
    mt = array (0, c (m, m, m))
    for (i in nn) {
      mt = mt + outer (outer (Y[i, ], Y[i, ]), Y[i,])
    }
  }
  if(p == 4){
    mt = array (0, c (m, m, m, m))
    for (i in nn) {
      mt = mt + outer (outer (Y[i, ], Y[i, ]), outer (Y[i, ], Y[i, ]))
    }
  }

  return(mt)
}
