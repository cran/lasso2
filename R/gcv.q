###  Copyright (C) 1998
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
gcv <- function(object, ...) UseMethod("gcv")

gcv.l1ce <- function(object,
                     type = c("OPT", "Tibshirani"),
                     gen.inverse.diag = 0, ...)
{
  type <- match.arg(type)

  covdf <- vcov(object,type,gen.inverse.diag)
  rss  <- deviance(object)

  n <- length(object$fitted)
  res <- rss/((1-covdf$df[1]/n)^2*n)

  cbind(bound=object$bound, gcv=res)
}

gcv.l1celist <- function(object, type = c("OPT", "Tibshirani"),
                         gen.inverse.diag = 0, ...)
{
  type <- match.arg(type)

  covdf <- vcov(object,type,gen.inverse.diag)
  rss  <- deviance(object)

  rssdf <- do.call("rbind",lapply(covdf, "[[", "df"))[,1]
  n <- length(object[[1]]$fitted)
  res <- rss/((1-rssdf/n)^2*n)

  cbind(bound=aux(object), gcv = res)
}

