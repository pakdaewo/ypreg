#bound = 50; alpha = 0.05; tau = NULL; repnum = 5000;time.hr = NULL; L = NULL; U = NULL
# res <- ypreg(time, event, group, X, tie = FALSE)
# plot(res)
ypreg.default <- function(time, event, group, X = NULL, tie = TRUE, alpha = 0.05, time.hr = NULL, L = NULL, U = NULL, bound = 50, repnum = 5000, tau = NULL, ...) {

  match.call()

  y <- as.numeric(time)
  d <- as.numeric(event)
  z <- as.numeric(group)
  x <- as.matrix(X)
  xnames <- colnames(X)

  n <- length(y)
  p <- ncol(X)

  random.add <- FALSE
  if (tie == TRUE) {
    if (length(y) != length(unique(y))) {
      JitterValue <- runif(n, 0, 1e-4)
      y <- y + JitterValue
      random.add <- TRUE
    }
  }

  I <- order(y)
  oy <- y[I]
  od <- d[I]
  oz <- z[I]
  ox <- x[I,]
  owz <- cbind(oz, ox)

  k <- ncol(owz)

  mean.x <- colMeans(x)
  sd.x <- apply(x, 2, sd)

  if (any(c(mean.x > 1, sd.x > 1.5))) stop("Please standadize covariates X first")

  # cox regression
  best_cox <- fn_est_cox(oy, od, owz)
  hess_cox <- fn_hess_cox(best_cox, oy, od, owz)
  summ_cox <- fn_summ_cox(best_cox, hess_cox, oy, od, owz, alpha)

  res0_yp0 <- fn_est0_ypmodel0(oy, od, oz, b0 = c(0, 0))
  ini_b0 <- res0_yp0$minb

  res_yp0 <- fn_est_ypmodel0(oy, od, oz, ini_b0)
  best_b0 <- res_yp0$minb

  ini_b <- c(best_b0, best_cox[-1])
  res_ypx <- fn_est_ypmodelx(oy, od, owz, ini_b)
  best_ypx <- res_ypx$minb
  res_summ <- fn_summ(best_ypx, oy, od, owz, alpha, L, U, h = 10E-8, repnum, tau)

  if (!is.null(time.hr)) {
    res_hrci <- fn_hrci(time.hr, res_summ)
  } else {
    res_hrci <- NULL
  }


  results <- c()
  results$summ_cox <- summ_cox
  results$ini_b0 <- ini_b0
  results$best_b0 <- best_b0
  results$best_ypx <- best_ypx
  results$res_summ <- res_summ
  results$res_hrci <- res_hrci
  results$tau <- tau
  results$xnames <- xnames
  results$p <- p
  results$n <- n
  results$time.hr <- time.hr
  results$alpha <- alpha
  results$random.add <- random.add

  class(results) <- "ypreg"

  results

}

#
# library(survival)
# coxph(Surv(oy, od) ~ owz)
# time <- test1$time
# status <- test1$status
# x <- as.matrix(test1$x)
#
# rootSolve::multiroot(f = estfn, start = rep(0, 1), y = time, d = status, x = x)$root
