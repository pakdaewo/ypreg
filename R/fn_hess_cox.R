fn_hess_cox <- function(best, oy, od, owz, ...) {

  kc <- ncol(owz)
  n <- nrow(owz)

  m <- kc + 1
  nb <- matrix(0, 1, kc)
  ob <- nb + 0.1
  jh <- 1E-8

  rf0 <- exp(owz %*% best)
  s0 <- as.numeric(matrix(colSums(rf0), nrow = 1)[rep(1, n), ]  - apply(rf0, 2, cumsum) + rf0)
  s1 <- matrix(0, n, kc)
  for (i in 1:kc) {
    rf1 <- owz[,i] * rf0
    s1[,i] <- matrix(colSums(rf1), nrow = 1)[rep(1, n), ]  - apply(rf1, 2, cumsum) + rf1
  }

  s2mat <- array(0, dim = c(n, kc, kc))
  for (i in 1:kc) {
    for (j in 1:kc) {
      tem <- as.numeric(owz[,i] * owz[,j] * rf0)
      s2mat[, i, j] <- sum(tem) - cumsum(tem) + tem
    }
  }

  sig <- diag(kc)

  for (i in 1:kc) {
    for (j in 1:kc) {
      sig[i, j] <- sum(od * (s2mat[,i, j]/s0 - (s1[,i]/s0) * (s1[,j]/s0)))

    }
  }

  return(sig)

}
