###  Copyright (C) 1998
###  Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
###  Bill Venables <wvenable@stats.adelaide.edu.au>
### --> ../COPYRIGHT for more details

is.formula <- function(x) inherits(x, "formula")

merge.formula <- function(x,y, ...)
{
    if(!(is.formula(x) && length(x) == 3))
        stop("First argument is invalid")
    if(!is.formula(y)) stop("Second argument is invalid")
    if(length(list(...))) warning("extraneous arguments discarded")

    str <- paste(c(deparse(x[[2]]), "~",
                   deparse(x[[3]]), "+",
                   deparse(y[[length(y)]])), collapse = "")
    as.formula(str)
}
