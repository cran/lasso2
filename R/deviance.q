###  Copyright (C) 1998, 1999
###  Justin Lokhorst <jlokhors@stats.adelaide.edu.au>
###  Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
###  Bill Venables <wvenable@stats.adelaide.edu.au>
### --> ../COPYRIGHT for more details

deviance.l1ce <- function(object, ...){
  if(is.null(w <- object$weights)){
    sum(object$residuals^2)
  }else{
    sum(w * object$residuals^2)
  }
}

deviance.l1celist <- function(object, ...)
{
  res2 <- do.call("cbind", lapply(object, "[[", "residuals"))^2
  weights <- attr(object, "weights")
  apply(if(is.null(weights)) res2 else weights*res2, 2, sum)
}

deviance.gl1ce <- function(object, ...) object$deviance
