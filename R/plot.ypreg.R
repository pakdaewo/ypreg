plot.ypreg <- function(x, legend.loc = "bottomleft", ...) {

  match.call()
  dotlist <- list(...)

  hr <- x$res_summ$hr
  yt <- x$res_summ$yt
  nwt <- x$res_summ$nwt
  ca22 <- x$res_summ$ca22
  ld2 <- x$res_summ$ld2
  ud2 <- x$res_summ$ud2
  n <- x$n
  alpha <- x$alpha

  #original scale
  upp22 <- hr + ca22/sqrt(n) * nwt
  low22 <- hr - ca22/sqrt(n) * nwt

  # log-transformation
  upp3 <- exp(ca22/sqrt(n) * nwt/hr) * hr
  low3 <- exp(-ca22/sqrt(n) * nwt/hr) * hr

  xmax <- max(yt)
  ymax <- max(c(low22[ld2:ud2], upp3[ld2:ud2]))

  main.text <- paste("Hazard Ratio with", paste((1 - alpha) * 100, "% Confidence Intervals", sep = ""), collapse = "")

  if(is.null(dotlist$xlab)) dotlist$xlab = "Time"
  if(is.null(dotlist$ylab)) dotlist$ylab = "Hazard Ratio"
  if(is.null(dotlist$main)) dotlist$main = main.text
  if(is.null(dotlist$xlim)) dotlist$xlim = c(0, xmax)
  if(is.null(dotlist$ylim)) dotlist$ylim = c(0, ymax)
  if(is.null(dotlist$col)) dotlist$col = "blue"
  if(is.null(dotlist$lwd)) dotlist$lwd = 1.5
  if(is.null(dotlist$lty)) dotlist$lty = 2
  if(is.null(dotlist$type)) dotlist$type = 'l'
  if(is.null(dotlist$axes)) dotlist$axes = FALSE
  dotlist$x = yt[ld2:ud2]
  dotlist$y = hr[ld2:ud2]

  do.call(what = plot, args = dotlist)
  abline(h = 1, lty = 2, col = "black")
  lines(yt[ld2:ud2], upp22[ld2:ud2], type = 'l', lty = 1, lwd = 1.5, col = "red")
  lines(yt[ld2:ud2], low22[ld2:ud2], type = 'l', lty = 1, lwd = 1.5, col = "red")
  lines(yt[ld2:ud2], upp3[ld2:ud2], type = 'l', lty = 3, lwd = 1.5, col = "pink")
  lines(yt[ld2:ud2], low3[ld2:ud2], type = 'l', lty = 3, lwd = 1.5, col = "pink")


  x.max <- dotlist$xlim[2]
  y.max <- dotlist$ylim[2]

  if(dotlist$axes == FALSE) {
    x.axis0 <- signif(seq(from = 0, to = x.max, length.out = 6),2)
    x.axis <- c(x.axis0[x.axis0 <= signif(x.max, 2)], signif(x.max, 2))
    y.axis0 <- signif(seq(from = 0, to = y.max, length.out = 6),2)
    y.axis <- c(y.axis0[y.axis0 <= signif(y.max, 2)], signif(y.max, 2))
    axis(1, x.axis)
    axis(2, c(1, y.axis))
  }

  legend(legend.loc, legend = c("hazard ratio", "CI(original scale)", "CI(log-transformed)"),
         col = c(dotlist$col, "red", "pink"), lty = c(dotlist$lty, 1, 3),
         lwd = c(dotlist$lwd, 1, 1), box.lty = 0)


}
