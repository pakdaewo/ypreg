#T1 = 1; T2 = 7; L = 0; U = 8; h = 10E-8
fn_summ <- function(best, oy, od, owz, alpha, L, U, h = 10E-8, repnum, tau, ...) {

  zalpha <- qnorm(1 - alpha/2, lower.tail = T)
  oz <- owz[, 1]
  ox <- owz[, -1]
  n <- length(oy)

  k <- ncol(owz) - 1
  m <- k + 3
  imat <- cbind(rep(0, k + 2), h * diag(k + 2))
  b <- best + imat

  u <- fn_ux(b, od, owz)
  gama1 <- exp(-oz %x% t(b[1,]) - ox %*% b[3:(k+2), ])
  gama2 <- exp(-oz %x% t(b[2,]) - ox %*% b[3:(k+2), ])

  sm <- n:1
  dlam <- (od/sm) %x% t(rep(1, m))

  lam2 <- apply(dlam * gama2, 2, cumsum)
  p <- exp(-lam2)
  pl <- rbind(matrix(1, ncol = m), p[1:(n-1), ])
  ru <- apply(pl * dlam * gama1, 2, cumsum)/p

  dd <- which(od != 0)
  tk <- length(dd)

  pq <- matrix(0, k + 2, m - 1)
  pr <- matrix(0, tk, m - 1)

  for (i in 1:(m - 1)) {
    pq[, i] <- (u[, i+1] - u[, 1])/h
    pr[, i] <- (ru[dd, i+1] - ru[dd, 1])/h
  }

  oyl <- c(0, oy[1:(n - 1)])
  po  <-p[dd, 1]
  plo <- pl[dd, 1]
  dp <- po - plo
  r <- ru[dd, 1]
  rl <- c(0, r[1:(tk - 1)])
  dr <- r - rl
  g1 <- gama1[, 1]
  g2 <- gama2[, 1]
  g1d <- g1[dd]
  g2d <- g2[dd]

  bt <- exp(-best)
  den <- bt[1] + bt[2] * r
  hr <- (1 + r)/den
  tem <- (bt[1] - bt[2])/den^2
  a <- matrix(1, tk, k + 2)
  a[, 1] <- tem * pr[, 1] + bt[1] * (1 + r)/den^2
  a[, 2] <- tem * pr[, 2] + bt[2] * r * (1 + r)/den^2
  for (i in 3:(k+2)) {
    a[, i] <- tem * pr[, i]
  }
  br <- tem/po

  if(is.null(tau)) tau <- oy[n]
  ttk <- sum(oy[dd] <= tau)
  if(ttk == 0) {
    print("There are no events before tau")
  }

  dyt <- c(oy[dd[1:ttk]], tau) - c(0, oy[dd[1:ttk]])
  cua <- colSums(rbind(a[1, ], a[1:ttk, ]) * dyt)
  avhr <- sum(c(hr[1], hr[1:ttk]) * dyt/tau)
  cbr0 <- cumsum(c(br[1], br[1:ttk]) * dyt)
  cbr <- cbr0[ttk + 1] - cbr0[1:ttk]

  kt <- sm[dd]
  inrs <- matrix(0, tk, k+2)
  inrsw <- matrix(0, tk, k+2)
  inrw <- matrix(0, tk, 1)

  for (ti in 1:tk) {
    yk <- oy >= oy[dd[ti]]
    dk <- g1 + r[ti] * g2
    tek <- yk/dk
    inrs[ti, 1] <- sum(oz * (g1 * tek/dk))
    inrs[ti, 2] <- r[ti] * sum(oz * g2 * tek/dk)
    inrsw[ti,1] <- sum(oz*(g1 * g2 * tek/dk^2))
    inrsw[ti,2] <- r[ti] * sum(oz * (g2^2 *tek/dk^2))

    for (i in 1:k) {
      inrs[ti, i+2] <- sum(ox[, i] * tek)
      inrsw[ti, i+2] <- sum(ox[, i] * g2 * tek/dk)
    }
  }


  inr <- matrix(0, tk, k + 2)
  rmul <- plo/kt

  for (i in 1:(k + 2)) {
    inr[,i] <- inrsw[,i]*dr/po + inrs[,i] * dp/po/plo
    inr[,i] <- inr[,i] + sum(inr[,i]) - cumsum(inr[,i])
    inr[,i] <- inr[,i] * rmul
  }

  di <- g1d + g2d * r
  ozd <- oz[dd]
  xid <- matrix(0, tk, k + 2)
  xid[, 1] <- ozd * g1d / di - inrs[, 1] * di/kt + inr[, 1] * di
  xid[, 2] <- ozd * g2d * r/ di - inrs[, 2] * di/kt + inr[, 2] * di

  for (i in 1:k) {
    xid[, i+2] <- ox[dd, i] - inrs[, i+2] * di/kt + inr[, i+2] * di
  }

  cid <- plo/kt * di
  inq <- solve(-pq/n)

  vmb <- matrix(1, k+2, repnum)
  wtild0 <- matrix(0, tk, repnum)

  for (irep in 1:repnum) {
    g <- rnorm(tk)
    vmb[, irep] <- as.numeric(inq %*% (t(xid) %*% g)/sqrt(n))
    wtild0[,irep] <- as.numeric(a %*% vmb[,irep]) + sqrt(n) * br * cumsum(cid * g)
  }

  wt0 <- c()
  for (tki in 1:tk) {
    wt0[tki] <- sd(wtild0[tki, ])
  }

  ovmb <- apply(vmb, 1, sort)
  stv <- apply(vmb, 1, function(x, ...) sd(x)/sqrt(n))
  exbeciva <- cbind(exp(best), exp(best - zalpha * stv), exp(best + zalpha * stv))

  zv <- best/stv
  pval <- 2 * pnorm(abs(zv), lower.tail = F)
  exbcip <- cbind(exbeciva, stv, zv, pval)

  ivn <- 50
  yt <- oy[dd]
  tp <- matrix(0, ivn, 1)
  for (i in 1:ivn) {
    tp[i] <- sum(yt <= tau*i/ivn)
  }
  lotp <- which(tp > 0)
  ntp <- tp[lotp]
  tpk <- length(ntp)
  nwt0 <- wt0[ntp]

  nwt <- c()
  for (i in 1:tk) {
    ctk <- sum(yt[ntp] <= oy[dd[i]])
    if(ctk <= 1) {
      nwt[i] <- nwt0[1]
    } else {
      if (ctk < tpk) {
        nwt[i] <- (yt[i] - yt[ntp[ctk]]) *
          (nwt0[ctk+1] - nwt0[ctk])/(yt[ntp[ctk+1]] - yt[ntp[ctk]]) + nwt0[ctk];
      } else {
        nwt[i] <- nwt0[ctk]
      }
    }
  }

  if (is.null(L)) L <- min(yt)
  if (is.null(U)) U <- max(yt)

  ld2 <- sum(yt <= L)
  ud2 <- sum(yt <= U)

  if (ld2 == 0) ld2 = 1

  mb22 <- matrix(1, repnum, 1)
  avhs <- matrix(1, repnum, 1)

  for (irep in 1:repnum) {
    g <- rnorm(tk)
    vmb[, irep] <- as.numeric(inq %*% (t(xid) %*% g)/sqrt(n))
    wtild <- a %*% vmb[, irep] + sqrt(n) * br * cumsum(cid * g)
    avhs[irep] <- sum(cua * vmb[, irep]) + sqrt(n) * sum(cbr[1:ttk] * (cid[1:ttk] * g[1:ttk]))
    bn2 <- wtild/nwt
    mb22[irep] <- max(abs(bn2[ld2:ud2]))
  }

  stanav <- sd(avhs)/tau/sqrt(n)
  zv1 <- (1 - avhr)/stanav
  pav <- 2 * pnorm(abs(zv1), lower.tail = F)
  avhrcip <- c(avhr, stanav, avhr - zalpha * stanav, avhr + zalpha * stanav, zv1, pav)

  m22 <- sort(mb22)
  ca22 <- m22[repnum * (1 - alpha)] #0.95

  colnames(exbcip) <- c("expb", "lower.expb", "upper.expb", "SE", "z", "pvalue")
  names(avhrcip) <- c("avhr", "SE", "lower.avhr", "upper.avhr", "z", "pvalue")

  results <- list()
  results$L <- L
  results$U <- U
  results$n <- n
  results$yt <- yt
  results$hr <- hr
  results$zalpha <- zalpha
  results$exbcip <- exbcip
  results$avhrcip <- avhrcip
  results$ca22 <- ca22
  results$nwt <- nwt
  results$ld2 <- ld2
  results$ud2 <- ud2

  return(results)

}
