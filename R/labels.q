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
labels.l1ce <- function(object, ...)
{
    TL <- object$terms
    if(!is.null(TL)) {
        TL <- attr(TL, "term.labels")
        TA <- object$assign
        if(!is.null(TA)) {
            TA <- names(TA)#c
            TL <- TL[match(TA, TL, 0)]
        }
    }
    TL
}

labels.l1celist <- function(object, ...)
{
  TL <- attr(object, "terms")
  if(!is.null(TL)) {
    TL <- attr(TL, "term.labels")
    TA <- object$assign
    if(!is.null(TA)) {
      TA <- names(TA)
      TL <- TL[match(TA, TL, 0)]
    }
  }
  TL
}


