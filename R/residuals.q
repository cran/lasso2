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
residuals.l1ce <-
    function(object, type = c("working", "pearson", "deviance"), ...)
{
  type <- match.arg(type)
  switch(type,
         working = object$residuals,
         pearson =, deviance =
         if(is.null(object$weights)) object$residuals else
         residuals.glm(object, "pearson")
         )
}

residuals.l1celist <-
    function(object, type = c("working", "pearson", "deviance"), ...)
{
  type <- match.arg(type)
  resid <- do.call("cbind", lapply(object, "[[", "residuals"))
  weights <- attr(object, "weights")
  switch(type,
         working = resid,
         pearson =, deviance =
         if(is.null(weights)) resid else sqrt(weights)*resid
         )
}

residuals.gl1ce <-
    function(object, type = c("deviance", "pearson", "working", "response"),
             ...)
{
  type <- match.arg(type)
  switch(type,
         working = object$residuals,
         pearson = sqrt(object$weights) * object$residuals,
         ## MM: the above does *not* seem the same as for summary.glm()
         deviance = {
           y <- object$y
           mu <- object$fitted
           family <- object$family
           w <- object$prior.weights
           if(is.null(w))
             w <- rep(1, length(mu))
           family$deviance(mu, y, w, residuals=TRUE)
         },
         response = object$y - fitted(object) )
}
