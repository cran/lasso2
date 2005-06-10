###  Copyright (C) 1998
###  Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
###  Bill Venables <wvenable@stats.adelaide.edu.au>
### --> ../COPYRIGHT for more details

if(is.R()) {
### -- this is also used in gl1ce() hence ''package global'' :

    ## Orig 2.1 version (for S+) calls  qr.rtr.inv(.) which is not in R;
    ## Is only used to return the (R'R)^{-1} where MM thinks should
    ## rather return the QR object (or just its  $qr and $rank) !!
    ## It's only summary() or vcov() which needs this, and
    ## they really can compute it *then* instead of *now*
    qr.rtr.inv <- function(qr)
    {
        if(is.null(R <- qr$qr))
            stop("argument is not a valid \"qr\" object")
        p <- qr$rank
        rinv <- backsolve(R[1:p, 1:p, drop = FALSE], diag(p))
        r <- rinv %*% t(rinv)
        nm <- (dimnames(R)[[2]])[1:p]
        dimnames(r) <- list(nm, nm)
        r
    }
}

l1ce <- function(formula, data = sys.parent(), weights, subset, na.action,
                 sweep.out = ~ 1,
		 x = FALSE, y = FALSE, contrasts = NULL, standardize = TRUE,
		 trace = FALSE,
                 guess.constrained.coefficients = double(p),
		 bound = 0.5, absolute.t = FALSE)
{
    call <- match.call()
    m <- match.call(expand = FALSE)
    m$sweep.out <- m$x <- m$y <- m$contrasts <-
        m$standardize <- m$guess.constrained.coefficients <- m$trace <-
            m$bound <- m$absolute.t <- NULL

    something.to.sweep.out <- !is.null(sweep.out)
    if(something.to.sweep.out)
        m$formula <- merge.formula(formula,sweep.out)

    m[[1]] <- as.name("model.frame")

    if(!missing(data) && !is.data.frame(data)) {
        m$data <- data <- as.data.frame(data)
        warning(paste(deparse(substitute(data)), "is not a dataframe"))
    }

    m <- eval(m, parent.frame())# not just 'data'
    weights <- model.extract(m, weights)
    Y <- model.extract(m, response)
    Terms <- terms(formula, data = data)
    X <- model.matrix(Terms, m, contrasts)
    X.names <- dimnames(X)[[2]]

    trace <- as.logical(trace)
    X.to.C <- X
    Y.to.C <- Y

    if(weighted <- length(weights)) {
        Y.to.C <- Y.to.C * (w <- sqrt(weights))
        X.to.C <- X.to.C * w
    }

    if(something.to.sweep.out) {
        sweep.out[[3]] <- sweep.out[[2]]
        sweep.out[[2]] <- formula[[2]]
        X.sweep.out <- model.matrix(terms(sweep.out, data = data),
                                    m, contrasts)
        X.sweep.out.names <- dimnames(X.sweep.out)[[2]]
        if(weighted)
            X.sweep.out <- X.sweep.out * w

        name.matches <- match(X.sweep.out.names,X.names)
        all.matched <- !any(is.na(name.matches))
        if(!all.matched)
            warning("Variables in 'sweep.out' are not a subset of variables in 'formula'")

        name.matches <- name.matches[!is.na(name.matches)]
        if(some.matched <- length(name.matches)) {
            X.to.C <- X.to.C[,-name.matches,drop = FALSE]
            X.names <- X.names[-name.matches]
            if(!length(X.to.C))
                stop("you cannot sweep out all the variables")
        }

        X.so.qr <- qr(X.sweep.out)
        if( X.so.qr$rank != ncol(X.sweep.out) )
            warning("Matrix built from variables in 'sweep.out' is rank deficient")

        X.so.coefficients  <- qr.coef  (X.so.qr, Y.to.C)
        X.so.X     <- qr.coef  (X.so.qr, X.to.C)
        X.so.Y.fit <- qr.fitted(X.so.qr, Y.to.C)
        X.to.C     <- qr.resid (X.so.qr, X.to.C)
        Y.to.C     <- qr.resid (X.so.qr, Y.to.C)
    }

    if(standardize) {
        X.to.C.stds <- sqrt(apply(X.to.C,2,var))
        if(any(i <- X.to.C.stds < sqrt(100 * .Machine$double.eps) *
               max(X.to.C.stds)))
            stop("Transformed variable matrix has constant column ", which(i),
                 "; set standardize = FALSE")
        X.to.C <- sweep(X.to.C, 2, X.to.C.stds, "/")
    }

    n <- nrow(X.to.C)
    p <- ncol(X.to.C)

    if(!absolute.t) {
	rnk <- (X.to.C.qr <- qr(X.to.C))$rank
	if(rnk != p && p < n)
	    warning("X Matrix (transformed variables) has rank ",rnk,
		    " < p = ",p,", i.e., is deficient")
	else if (rnk == 0)
	    stop("Matrix built from transformed variables is null matrix")
	t0 <- sum(abs(qr.coef(X.to.C.qr, Y.to.C))[1:rnk])
	if (any(bound > 1))
	    stop("'bound'(s) must be between 0 and 1 if 'absolute.t' is false")

	bound <- (relative.bound <- bound) * t0
    }

    if(any(bound < 0))
        stop("'bound'(s) must be non negative")

    if( length(guess.constrained.coefficients) != p )
        stop("invalid argument for 'guess.constrained.coefficients'")

    keep <- c("coefficients", "fitted.values", "residuals", "success",
              "Lagrangian", "bound")

    if (1 == (num.bound <- length(bound)) ) { ## 1 bound ----------------------

        fit <- .C("lasso",
                  X = as.double(X.to.C),
                  n = n, p = p,
                  bound = as.double(bound),
                  coefficients = as.double(guess.constrained.coefficients),
                  Y = as.double(Y.to.C),
                  fitted.values = double(n),
                  residuals     = double(n),
                  Lagrangian    = double(1),
                  success = integer(1),
                  trace   = trace,
                  assub   = FALSE,
                  PACKAGE = "lasso2")[keep]

        if (fit$success < 0)
            stop("Oops, something went wrong in .C(\"lasso\",..): ",fit$success)
        ## else drop it:
        fit$success <- NULL

        fit$xtx <- crossprod(X.to.C)
        fit$xtr <- crossprod(X.to.C,fit$residuals)
        fit$constrained.coefficients <- fit$coefficients
        names(fit$coefficients) <- names(fit$constrained.coefficients) <-
            dimnames(X.to.C)[[2]]

        if(standardize) {
            fit$X.stds <- X.to.C.stds
            fit$coefficients   <- fit$coefficients / X.to.C.stds
        }

        if(something.to.sweep.out) {
            fit$fitted.values <- fit$fitted.values + X.so.Y.fit
            X.so.coefficients <- X.so.coefficients - X.so.X %*% fit$coefficients

            tmp <- fit$coefficients
            fit$coefficients <- rep(0,ncol(X))
            names(fit$coefficients) <- dimnames(X)[[2]]
            fit$coefficients[X.names] <- tmp
            names(X.so.coefficients) <- X.sweep.out.names
            ind <- names(fit$coefficients[name.matches])
            ##orig 1999 code : fit$coefficients[name.matches] <- X.so.coeff[ind]
            fit$coefficients[name.matches] <- X.so.coefficients[ind]

            fit$sweep.out <- list(sweep.out = sweep.out,
                                  X.sweep.out.names = X.sweep.out.names,
                                  name.matches = name.matches,
                                  all.matched = all.matched,
                                  some.matched = some.matched,
                                  X.so.rtr.inv = qr.rtr.inv(X.so.qr),
                                  X.so.X = X.so.X)
        }

        if(weighted) {
            fit$weights <- weights
            if( any(weights == 0) ) {
                fit$fitted.values <- X %*% fit$coefficients
                fit$residuals <- Y - fit$fitted.values
            } else {
                fit$fitted.values <- fit$fitted.values / w
                fit$residuals <- fit$residuals / w
            }
        }

        fit$terms <- Terms
        fit$call <- call
        fit$contrasts <- attr(X, "contrasts")
        fit$assign <- attr(X, "assign")
        if(x) fit$x <- X
        if(y) fit$y <- Y
        if(!absolute.t)
            fit$relative.bound <- relative.bound
        structure(fit, class = "l1ce",
                  na.message = attr(m, "na.message"))

    }
    else { ##--------- more than 1 bound ---------------------------------

        ordered.bound <- order(bound)
        guess.constraint.coefficients <-
            rep(guess.constrained.coefficients,num.bound)

        res <- .C("mult_lasso",
                  X = as.double(X.to.C), n = n, p = p,
                  bound = as.double(bound[ordered.bound]),
                  l = num.bound,
                  coefficients = as.double(guess.constraint.coefficients),
                  Y = as.double(Y.to.C),
                  fitted.values = double(n*num.bound),
                  residuals     = double(n*num.bound),
                  Lagrangian    = double(num.bound),
                  success = integer(1),
                  trace   = trace,
                  PACKAGE = "lasso2")[keep]

        if (res$success < 0)
            stop("Oops, something went wrong in .C(\"mult_lasso\",..): ",
                 res$success)
        ## else drop it:
        res$success <- NULL

        total.fit <- vector("list", num.bound)
        ind1 <- 1:n
        ind2 <- 1:p
        for(i in 1:num.bound) {

            resi <- res$residuals[ind1]
            fit <- list(coefficients  = res$coefficients[ind2],
                        fitted.values = res$fitted.values[ind1],
                        residuals     = resi,
                        bound         = res$bound[i],
                        Lagrangian    = res$Lagrangian[i],
                        xtr           = crossprod(X.to.C, resi))#FIXED!
            names(fit$coefficients) <- dimnames(X.to.C)[[2]]
            fit$constrained.coefficients <- fit$coefficients

            if(standardize)
                fit$coefficients   <- fit$coefficients / X.to.C.stds

            if(something.to.sweep.out) {
                fit$fitted.values <- fit$fitted.values + X.so.Y.fit
                X.so.coef <- X.so.coefficients - X.so.X %*% fit$coefficients

                tmp <- fit$coefficients
                fit$coefficients <- rep(0,ncol(X))
                names(fit$coefficients) <- dimnames(X)[[2]]
                fit$coefficients[X.names] <- tmp
                names(X.so.coef) <- X.sweep.out.names
                ind <- names(fit$coefficients[name.matches])
                fit$coefficients[name.matches] <- X.so.coef[ind]
            }

            if(weighted) {
                if( any(weights == 0) ) {
                    fit$fitted.values <- X %*% fit$coefficients
                    fit$residuals <- Y - fit$fitted.values
                } else {
                    fit$fitted.values <- fit$fitted.values / w
                    fit$residuals <- fit$residuals / w
                }
            }

            if(!absolute.t)
                fit$relative.bound <- relative.bound[ordered.bound[i]]

            total.fit[[ordered.bound[i]]] <- fit
            ind1 <- ind1 + n
            ind2 <- ind2 + p
        }

        if(standardize) attr(total.fit, "X.stds") <- X.to.C.stds

        if(something.to.sweep.out) {
            attr(total.fit, "sweep.out") <-
                list(sweep.out = sweep.out,
                     X.sweep.out.names = X.sweep.out.names,
                     name.matches = name.matches,
                     all.matched = all.matched,
                     some.matched = some.matched,
                     X.so.rtr.inv = qr.rtr.inv(X.so.qr),
                     X.so.X = X.so.X)
        }
        attr(total.fit, "xtx") <- crossprod(X.to.C)
        if(weighted) attr(total.fit, "weights") <- weights
        attr(total.fit, "terms") <- Terms
        attr(total.fit, "call") <- call
        attr(total.fit, "contrasts") <- attr(X, "contrasts")
        attr(total.fit, "assign") <- attr(X, "assign")
        if(x) attr(total.fit, "x") <- X
        if(y) attr(total.fit, "y") <- Y
        structure(total.fit, class = "l1celist",
                  na.message = attr(m, "na.message"))
    }## multi-bound case
}
