###  Copyright (C) 1998
###  Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
###  Bill Venables <wvenable@stats.adelaide.edu.au>
### --> ../COPYRIGHT for more details

coef.l1ce <- function(object, all=TRUE, constrained=FALSE, ...)
{
    label <- if(constrained) "constrained.coefficients" else "coefficients"
    if(all)
        object[[label]]
    else
        object[[label]][abs(object[[label]]) > sqrt(.Machine$double.eps)]
}

coef.l1celist <- function(object, all=TRUE, constrained=FALSE, ...)
{
    label <- if(constrained) "constrained.coefficients" else "coefficients"
    if(all) {
        do.call("rbind", lapply(object, "[[", label))
    }
    else {
        junk <- do.call("rbind", lapply(object, "[[", label))
        ind <- apply(abs(junk), 2, sum)
        junk[, ind != 0]
    }
}

