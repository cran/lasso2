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
predict.l1ce <-
    function(object, newdata, type=c("response"), se.fit = FALSE, ...)
{
    type <- match.arg(type)
    if(missing(newdata))
        return(object$fitted)
    if(se.fit)
        stop("The `se.fit' argument is not currently implemented for l1ce objects")

    tt <- object$terms
    if(!inherits(tt, "terms"))
        stop("invalid terms component of  fit")
    offset <- attr(tt, "offset")
    intercept <- attr(tt, "intercept")

    if(missing(newdata)) {
        x <- model.matrix(object)
    } else if(!((is.atomic(newdata) && length(newdata) == 1
                && length(object$coef) != 1 && newdata > 0
                && (newdata - trunc(newdata) < .Machine$single.eps))
               | is.list(newdata))) {
        ## try and coerce newdata to look like the x matrix
        if (!is.null(offset)) {
            warning("Offset not included")
            offset <- NULL
        }
        TT <- length(object$coef)
        if(is.matrix(newdata) && ncol(newdata) == TT)
            x <- newdata
        else if(length(newdata) == TT)
            x <- matrix(newdata, 1, TT)
        else stop("Argument `newdata' is not a data frame, and cannot be coerced to an appropriate model matrix")
    } else {
        ## newdata is a list, data frame or frame number
        vv <- attr(tt, "term.labels")
        attr(tt, "term.labels") <- vv[ - attr(tt, "response")]
        x <- model.matrix(tt, newdata, object$contrasts)
        if(!is.null(offset))
            offset <- eval(attr(tt, "term.labels")[offset], newdata)
    }

    coefs <- coef(object)
    pred <- drop( x%*% coefs )
    if(!is.null(offset)) {
        if(missing(newdata)) {
            warning("Offset not included")
        }
        else {
            pred <- pred + offset
        }
    }
    pred
}

predict.gl1ce <-
    function(object, newdata, type=c("link", "response"), se.fit = FALSE, ...)
{
    type <- match.arg(type)
    if(!se.fit){
        if (missing(newdata)) {
            switch(type,
                   link = object$linear.predictors,
                   response = object$fitted)
        } else {
            switch(type,
                   response = {
                       linkinv <- family(object)[[if(is.R())"linkinv" else "inverse"]]
                       linkinv(NextMethod("predict"))
                   },
                   link = {
                       type <- "response"
                       NextMethod("predict")
                   })
        }
    } else {
        stop("The `se.fit' argument is not currently implemented for gl1ce objects")
    }
}
