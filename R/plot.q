###  Copyright (C) 1998
###  Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
###  Bill Venables <wvenable@stats.adelaide.edu.au>
### --> ../COPYRIGHT for more details
plot.l1celist <-
    function(x, plot=TRUE, all=TRUE, constrained=FALSE,
             type = "b", xlab = "bounds", ylab = "coeff | bounds", ...)
{
  bounds <- aux(x)
  mat.of.coef <- coef(x, all=all, constrained = constrained)

  if(plot) {
    matplot(bounds[, 1], mat.of.coef,
            type = type, xlab = xlab, ylab = ylab, ...)
    rug(bounds[, 1])
  }
  invisible(list(bounds=bounds, mat.of.coef=mat.of.coef))
}
