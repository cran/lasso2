###  Copyright (C) 1998
###  Justin Lokhorst <jlokhors@stats.adelaide.edu.au>
###  Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
###  Bill Venables <wvenable@stats.adelaide.edu.au>
###
###  This library is free software; you can redistribute it and/or
###  modify it under the terms of the GNU Library General Public
###  License as published by the Free Software Foundation; either
###  version 2 of the License, or (at your option) any later version.
###
###  This library is distributed in the hope that it will be useful,
###  but WITHOUT ANY WARRANTY; without even the implied warranty of
###  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
###  Library General Public License for more details.
###
###  You should have received a copy of the GNU Library General Public
###  License along with this library; if not, write to the Free Software
###  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
###  MA 02111-1307 USA
print.l1ce <- function(x, ...){

  if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }

  cat("\nCoefficients:\n")
  print(coef(x), ...)

  if(!is.null(rel <- x$relative.bound)){
    cat("\nThe relative L1 bound was      :", format(rel), "\n")
  }else{
    cat("\n")
  }
  cat("The absolute L1 bound was      :", format(x$bound), "\n")
  cat("The Lagrangian for the bound is: ",
      format(x$Lagrangian), "\n")

  if(!is.null(attr(x, "na.message")))
    cat("\n", attr(x, "na.message"), "\n")

  invisible(x)
}

print.l1celist <- function(x, ...){

  if(!is.null(cl <- attr(x, "call"))) {
    cat("Call:\n")
    dput(cl)
  }

  mat.of.coef <- coef(x)
  cat("\nCoefficients:\n")
  print(mat.of.coef, ...)

  x.aux <- aux(x)
  if( ncol(x.aux) == 2 ){
    cat("\nAbsolute L1 bounds and the Lagrangians:\n")
    print(x.aux, ...)
  }else{
    cat("\nRelative and absolute L1 bounds and the Lagrangians:\n")
    print(x.aux, ...)
  }

  if(!is.null(attr(x, "na.message")))
    cat("\n", attr(x, "na.message"), "\n")

  invisible(x)
}

print.summary.l1ce <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:\n ")
  dput(x$call)
  resid <- x$residuals
  df <- x$df
  rdf <- df[2]
  if(rdf > 5) {
    cat("Residuals")
    if(any(na <- is.na(resid)))
        cat(" (",sum(na), "NA's)")
    cat(":\n")
    qnam <- c("Min", "1Q", "Median", "3Q", "Max")
    if(length(dim(resid)) == 2) {
      rq <- apply(t(resid), 1, quantile, na.rm = TRUE)
      dimnames(rq) <- list(qnam, dimnames(resid)[[2]])
    }
    else {
        rq <- quantile(resid[!na])
        names(rq) <- qnam
    }
    print(rq, digits = digits, ...)
  }
    else if(rdf > 0) {
      cat("Residuals:\n")
      print(resid, digits = digits, ...)
    }
  cat("\nCoefficients:\n")
  print(format(round(x$coef, digits = digits)), quote = FALSE, ...)
  if(x$sigma.provided){
    cat("\nResidual standard error:",
        format(signif(x$sigma, digits)),
        "was provided to \"summary.l1ce\"\n")
  }else{
    cat("\nResidual standard error:",
        format(signif(x$sigma, digits)),
        "on", format(signif(rdf,digits)), "degrees of freedom\n")
  }
  if(!is.null(rel <- x$relative.bound))
    cat("The relative L1 bound was      :", format(rel), "\n")
  cat("The absolute L1 bound was      :", format(x$bound), "\n")
  cat("The Lagrangian for the bound is: ",
      format(x$Lagrangian), "\n")

  correl <- x$correlation
  if(!is.null(correl)) {
    p <- dim(correl)[2]
    if(p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      ll <- lower.tri(correl)
      correl[ll] <- format(round(correl[ll], digits))
      correl[!ll] <- ""
      print(correl[-1, -p, drop = FALSE], quote = FALSE, digits = digits, ...)
    }
  }
  invisible(x)

}

print.gl1ce <- function(x, ...) {
  if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }

  cat("\nCoefficients:\n")
  print(coef(x), ...)

  cat("\nFamily:\n")
  print(x$family, ...)

  cat("\nThe absolute L1 bound was       : ", format(x$bound), "\n")
  cat("The Lagrangian for the bound is : ", format(x$Lagrangian), "\n")

  if(!is.null(attr(x, "na.message")))
    cat("\n", attr(x, "na.message"), "\n")

  invisible(x)
}

print.summary.gl1ce <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:\n")
  dput(x$call)
  resid <- residuals(x)
  cat("\nDeviance Residuals")
  qnam <- c("Min", "1Q", "Median", "3Q", "Max")
  if(any(na <- is.na(resid)))
      cat(" (",sum(na), "NA's)")
  cat(":\n")
  qnam <- c("Min", "1Q", "Median", "3Q", "Max")
  if(length(dim(resid)) == 2) {
      rq <- apply(t(resid), 1, quantile, na.rm = TRUE)
      dimnames(rq) <- list(qnam, dimnames(resid)[[2]])
  }
  else {
      rq <- quantile(resid[!na])
      names(rq) <- qnam
  }
  print(rq, digits = digits, ...)

  ccf <- format(round(x$coeff, digits = digits))
  ccf[x$coeff == 0.] <- "0" # visualize exact zeroes
  cat("\nCoefficients\n")
  print(ccf, quote = FALSE, ...)

  cat("\nFamily:\n")
  print(x$family, ...)

  cat("\nThe absolute L1 bound was     : ", format(x$bound), "\n")
  cat("The Lagrangian for the bound is : ", format(x$Lagrangian), "\n")

  cat("\nNumber of Iterations:", format(trunc(x$iter)), "\n")

  invisible(x)
}

