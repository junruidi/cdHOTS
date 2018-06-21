#' @title Construct Higher Order Tensorian Cumulants
#' @description 3rd and 4th order tensorian cumulants
#' @param Y data.frame or matrix of multivariate data of dimension n*p
#' @param p what is the order to construct
#'
#'
#' @return 3rd or 4th order cumulant tensors
#'
#' @export
#'
#'
#' @examples
#' data(dat)
#' tensor_3 = tensor_cumulant(Y = dat, p = 3)


tensor_cumulant = function (Y, p) {
  Y = as.matrix(Y)

  if(!p %in% c(3,4)){
    stop("p has to be either 3 or 4")
  }


  n = nrow (Y)
  m = ncol (Y)
  mm = 1 : m
  nn = 1 : n
  r1 = colSums (Y) / n
  r2 = crossprod (Y) / n

  if(p == 3){
    r3 = array (0, c (m, m, m))
    for (i in nn) {
      r3 = r3 + outer (outer (Y [i, ], Y [i, ]), Y [i,])
    }

    r3 = r3 / n
    c3 = r3
    for (i in mm) for (j in mm) for (k in mm) {
      s3 = r3 [i, j, k]
      s21 = r2 [i, j] * r1 [k] + r2 [i, k] * r1 [j] + r2[j, k] * r1 [i]
      s111 = r1 [i] * r1 [j] * r1 [k]
      c3 [i, j, k] = s3 - s21 + 2 * s111
    }
    cmu = c3
  }

  if(p == 4){
    r4 = array (0, c (m, m, m, m))
    for (i in nn) {
      r3 = r3 + outer (outer (Y [i, ], Y [i, ]), Y [i,])
      r4 = r4 + outer (outer (Y [i, ], Y [i, ]), outer (Y [i, ], Y [i, ]))
    }

    r4 = r4 / n
    c4 = r4
    for (i in mm) for (j in mm) for (k in mm) for (l in mm){
      s4 = r4 [i, j, k, l]
      s31 = r3 [i, j, k] * r1 [l] + r3 [i, j, l] * r1 [k] + r3 [i, k, l] * r1 [j] + r3 [j, k, l] * r1 [i]
      s22 = r2 [i, j] * r2 [k, l] + r2 [i, k] * r2 [j, l] + r2 [j, k] * r2 [i, l]
      s211 = r2 [i, j] * r1 [k] * r1 [l] + r2 [i, k] *
        r1 [j] * r1 [l] + r2 [i, l] * r1 [k] * r1 [j] +
        r2 [j, k] * r1 [i] * r1 [l] + r2 [j, l] * r1 [i] *
        r1 [k] + r2 [k, l] * r1 [i] * r1 [j]
      s1111 = r1 [i] * r1 [j] * r1 [k] * r1 [l]

      c4 [i, j, k, l] = s4 - s31 - s22 + 2 * s211 - 6 * s1111
    }

    cmu = c4
  }

  return (cmu)
}
