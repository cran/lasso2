###  Copyright (C) 1998, 1999
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
summary.l1ce <- function(object,
                         correlation = TRUE,
                         type = c("OPT", "Tibshirani"),
                         gen.inverse.diag = 0,
                         sigma = NULL, ...)
{
  type <- match.arg(type)

  coef <- coef(object)
  cnames <- labels(coef)
  resid <- resid(object)
  fv <- fitted(object)

  covdf <- vcov(object,type,gen.inverse.diag)
  sigma.provided <- !missing(sigma)
  if(!sigma.provided)
    sigma <- sqrt(deviance(object)/covdf$df[2])

  se <- as.vector(sqrt(diag(covdf$cov.unscaled)))

  correl <-
      if(correlation) {
          p  <- length(se)
          correl <- covdf$cov.unscaled/se
          correl/se[rep(1:p,rep(p,p))]
      } ## else NULL

  coef <- array(coef, c(p, 4))
  dimnames(coef) <- list(cnames,
                         c("Value", "Std. Error", "Z score", "Pr(>|Z|)"))
  coef[, 2] <- se %o% sigma
  coef[, 3] <- coef[, 1]/coef[, 2]
  coef[, 4] <- 2*(1-pnorm(abs(coef[,3])))

  keep <- c("call", "terms", "bound", "relative.bound", "Lagrangian")
  object <- object[keep[!is.na(match(keep,names(object)))]]
  object$residuals <- resid
  object$coefficients <- coef
  object$sigma <- sigma
  object$sigma.provided <- sigma.provided
  object$df <- covdf$df
  object$cov.unscaled <- covdf$cov.unscaled
  object$correlation <- correl
  class(object) <- "summary.l1ce"
  object
}

summary.gl1ce <- function(object, dispersion = NULL, correlation=FALSE, ...)
{
  if(correlation)
    stop("The `correlation' argument is not yet implemented for gl1ce objects")
  ## else
  coef <- coef(object)
  if(is.null(cnames <- names(coef)))
      cnames <- c("(Intercept)", labels(object))
  resid <- residuals(object, type="deviance")
  fv <- fitted(object)
  family <- object$family
  iter <- object$iter
  coef <- array(coef, length(coef))
  dimnames(coef) <- list(cnames)

  keep <- c("call", "terms", "bound", "Lagrangian")
  object <- object[keep[!is.na(match(keep, names(object)))]]
  object$residuals <- resid
  object$coefficients <- coef
  object$family <- family
  object$iter <- iter
  class(object) <- "summary.gl1ce"

  object
}
