#' third order mgf
#' @param x matrix
#' @return third order
#' @example
#' mg_3 = mg3(x)
#' @export
#'
mg3<- function (x) {
  n <- nrow (x)
  m <- ncol (x)
  mm <- 1 : m
  nn <- 1 : n
  c3 <- array (0, c (m, m, m))
  for (i in nn) {
    c3 <- c3 + outer (outer (x [i, ], x [i, ]), x [i,])
  }
  return (c3)
}
