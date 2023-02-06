fn_est_cox <- function(oy, od, owz, ...) {

  kc <- ncol(owz)
  n <- nrow(owz)

  m <- kc + 1
  nb <- matrix(0, 1, kc)
  ob <- nb + 0.1
  jh <- 1E-8

  mod <- mean(od)
  nm <- 50*(mod < .05) + 25*(mod >= .05)*(mod < .1) +
    15*(mod >= .1)*(mod < .2) + 10*(mod >= .2)*(mod < .5) + 7*(mod >= .5);

  ittm <- 100
  cn <- 0
  stoprule <- FALSE
  while (stoprule == FALSE & cn < ittm) {

    b <- as.numeric(nb) + cbind(0, diag(jh, 4))
    rf <- exp(owz %*% b)
    kk0 <- matrix(colSums(rf), nrow = 1)[rep(1, n), ]  - apply(rf, 2, cumsum) + rf
    u <- matrix(NA, kc, m)
    for (i in 1:kc) {
      rf1 <- owz[,i] * rf
      kk1 <- matrix(colSums(rf1), nrow = 1)[rep(1, n), ]  - apply(rf1, 2, cumsum) + rf1
      u[i,] <- as.numeric(t(-owz[,i]) %*% od)  + t(od) %*% (kk1/kk0)
    }

    h <- diag(kc)
    for (i in 1:kc) {
      h[,i] <- (u[,i+1] - u[,1])/jh
    }

    ep <- 0.0001
    cr <- 0

    while ((min(abs(svd(h)$d)) < 1.e-4) & (cr<100)) {
      ep <- 2 * ep
      h <- h + ep * diag(kc)
    }

    po <- -as.numeric(solve(h) %*% u[, 1])

    cn <- cn + 1
    t <- 1
    bc <- nb + t*po
    idb <- sum((abs(bc) > log(nm)))
    tb <- 0

    while ((idb>0)&(tb<20)) {
      tb <- tb + 1
      t <- t * 0.618
      bc <-  nb + t*po
      idb <- sum((abs(bc) > log(nm)))
    }

    nb <- bc
    stoprule <- sum(abs(round(1000*ob) - round(1000*nb))) == 0

    ob <- nb

    # os <- min(sum(u^2))
    # of <- os
  }


  return(as.numeric(nb))

}
