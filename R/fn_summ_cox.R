fn_summ_cox <- function(best, hess_cox, oy, od, owz, alpha, ...) {

  bmat <- solve(hess_cox)
  betav <- diag(bmat)
  za <- qnorm(1 - alpha/2)
  bse <- sqrt(betav)
  up <- exp(best + za * bse)
  dn <- exp(best - za * bse)
  zv <- best/bse
  pval <- 2 * pnorm(abs(zv), lower.tail = FALSE)

  tm0 <- as.numeric(owz %*% best)
  tem <- exp(tm0)
  tem2 <- log(tem + sum(tem) - cumsum(tem))
  logl <- sum(od * (tm0 - tem2))

  mres <- cbind(est = best, se = bse, exp_best = exp(best), dn, up, za, pval)
  results <- list()
  results$mres <- mres
  results$logl <- logl

  return(results)
}
