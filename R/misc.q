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

"[[.l1celist" <- function(x, ..., drop = TRUE)
{
  class(x) <- NULL
  val <- NextMethod("[[")
  if( !length(val) ) return (NULL)

  xattr <- attributes(x)
  nattr <- names(xattr)
  val[nattr] <- xattr[nattr]

  class(val) <- "l1ce"
  val
}

"[.l1celist" <- function(x, ..., drop = TRUE)
{
  class(x) <- NULL
  val <- NextMethod("[")
  junk <- do.call("c", lapply(val,length))
  val <- val[junk!=0]
  if(!length(val)) return (NULL)
  if(drop && length(val) == 1) {
    val <- val[[1]]
    xattr <- attributes(x)
    nattr <- names(xattr)
    val[nattr] <- xattr[nattr]
    class(val) <- "l1ce"
  } else {
    attributes(val) <- attributes(x)
    class(val) <- "l1celist"
  }
  val
}

aux <- function(object, ...) UseMethod("aux")

aux.l1celist <- function(object, ...)
{
  rbnd <- do.call("c", lapply(object, "[[", "relative.bound"))
  structure(cbind(rbnd,
                  do.call("c", lapply(object, "[[", "bound")),
                  do.call("c", lapply(object, "[[", "Lagrangian")) ),
            dimnames=list(NULL,
            c(if(!is.null(rbnd))"rel.bound", "abs.bound", "Lagrangian")))
}

## Trace
tr <- function(mat)
{
  dims <- dim(mat)
  if((length(dims) != 2) || dims[1] != dims[2])
    stop("This function is only defined for square matrices")

  sum(diag(mat))
}
