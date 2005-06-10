###  Copyright (C) 1999
###  Justin Lokhorst <jlokhors@stats.adelaide.edu.au>
###  Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
###  Bill Venables <wvenable@stats.adelaide.edu.au>
### --> ../COPYRIGHT for more details


###--- FIXME :    data = sys.parent() for S+
###--- =====    glm() in R uses a more general eval() !
gl1ce <- function(formula, data = sys.parent(), weights, subset, na.action,
                  family = gaussian, control = glm.control(...),
                  sweep.out = ~ 1,
                  x = FALSE, y = TRUE,
                  contrasts = NULL,
                  standardize = TRUE,
                  guess.constrained.coefficients = double(p),
                  bound = 0.5, ...)
{
  ret.x <- x
  ret.y <- y
  ## R's family functions partly *must* work with 'y' and
  ## binomial()$initialize has its own 'n' used in $aic() -- YUCK!! << Martin
  call <- match.call()
  mf <- match.call(expand = FALSE)

  mf$sweep.out <- mf$x <- mf$y <- mf$contrasts <- mf$standardize <-
    mf$guess.constrained.coefficients <- mf$trace <- mf$bound <-
      mf$family <- mf$control <- mf$... <- NULL

  something.to.sweep.out <- !is.null(sweep.out)

  if(something.to.sweep.out)
    mf$formula <- merge.formula(formula, sweep.out)

  mf[[1]] <- as.name("model.frame")
  if(!missing(data) && !is.data.frame(data)) {
    mf$data <- data <- as.data.frame(data)
    warning(paste(deparse(substitute(data)), "is not a dataframe"))
  }

  mf <- eval(mf, parent.frame())
  weights <- model.extract(mf, weights)
  y <- model.extract(mf, response)
  Terms <- terms(formula, data = data)
  X <- model.matrix(Terms, mf, contrasts)
  nobs <- nrow(X)# needed for R's family functions (do not use 'n', see above!)
  X.names <- dimnames(X)[[2]]
  offset <- model.extract(mf,offset)
  if(!is.numeric(offset)) offset <- 0 # for R (at least)

  if(!length(weights))
      weights <- rep(1, nobs)# not length(y) which might be 2 * nobs {binomial!}
  else if(any(weights < 0))
      stop("Negative weights not allowed")
  else if(any(weights == 0))
      stop("Some weigthts are 0.\n",
           "Rather, use 'subset= *' for observations that need to be excluded")

  do.trace <- as.logical(control$trace)

  if(!inherits(family, "family"))
      family <- family(family)

  ## cat("DEBUGGING gl1ce(): family =\n"); str(family)

  variance.function <- family$variance

  if(is.R()) {
      link     <- family$linkfun
      inv.link <- family$linkinv
      mu.eta   <- family$mu.eta    # == mu.eta(eta) = 1/deriv(mu(eta))
      ## this is straight from glm.fit() :
      valideta <- family$valideta
      if (is.null(valideta))
          valideta <- function(eta) TRUE
      validmu <- family$validmu
      if (is.null(validmu))
          validmu <- function(mu) TRUE
      ## eval(family$initialize)
      ## calculates mustart and may change y and weights and set n (!)
      eval(family$initialize)

      family$deviance <- function(mu, y, weights, residuals = FALSE)
      {
          dr <- dev.res(y, mu, weights)
          if(residuals) {
              d.res <- sqrt(pmax(dr, 0))
              ifelse(y > mu, d.res, -d.res)
          } else sum(dr)
      }
      ee <- new.env()
      assign("dev.res", family$dev.res, env = ee)
      environment(family$deviance) <- ee

      if (NCOL(y) > 1)
          stop("y must be univariate unless binomial")

      eta <-
          ## 'etastart' and 'start' (from glm.fit()) *are* NULL
          family$linkfun(mustart)

      mu <- inv.link(eta)
      eta <- eta - offset

  } else { ##-- S-plus --

      link <- family$link
      inv.link <- family$inverse
      deriv.link <- family$deriv

      mu <- # binomial hack -- does NOT work when y is cbind(k, n-k) !!
          if(family$family[[1]] == "Binomial")
              (y + (1/2))/(weights + 1)
          else (y + mean(y))/2

      eta <- link(mu) - offset
  }

  if(something.to.sweep.out) {
      X.sweep.out <- model.matrix(terms(sweep.out, data = data), mf, contrasts)
      X.sweep.out.names <- dimnames(X.sweep.out)[[2]]
      name.matches <- match(X.sweep.out.names, X.names)
      all.matched <- !any(is.na(name.matches))

      if(!all.matched)
      warning("Variables in 'sweep.out' are not a subset of those in 'formula'")

      name.matches <- name.matches[!is.na(name.matches)]
      if(some.matched <- length(name.matches)) {
          X.to.C <- X[, - name.matches, drop = FALSE]
          X.names <- X.names[ - name.matches]
          if(length(X.to.C) == 0)
              stop("Do you really want to sweep out all the variables?")
      } else {
          X.to.C <- X
      }
  } else {
      X.to.C <- X
  }

  stopifnot(nobs == nrow(X.to.C))
  p <- ncol(X.to.C)

  if(length(guess.constrained.coefficients) != p)
    stop("invalid argument for 'guess.constrained.coefficients'")

  keep <- c("coefficients", "fitted.values", "residuals", "success",
            "Lagrangian", "bound")
  check <- 1
  Cvector <- double(if(something.to.sweep.out) ncol(X) else ncol(X.to.C))
  continue <- TRUE
  j <- 0

  Y.to.C <- y

  if(do.trace) {
      cat("glqce() -- tracing -- before entering IRLS loop:",
        "\n------- nobs =",nobs,", p =",p,";  summary(weights):\n")
      print(summary(weights))
  }
  while(continue) {

    if(is.R()) { ## mu.eta() instead of deriv()
        iM <- mu.eta(eta)
        Y.to.C <- eta + (y - mu) / iM
        W <- sqrt(weights*iM^2 / variance.function(mu))
    } else { ## S+
        M <- deriv.link(mu)
        Y.to.C <- eta + (y - mu) * M
        W <- sqrt(weights/(M^2 * variance.function(mu)))
    }

    X.to.C.w <- X.to.C * W
    Y.to.C   <- Y.to.C * W

    if(something.to.sweep.out) {
      X.sweep.out.w <- X.sweep.out * W

      X.so.qr <- qr(X.sweep.out.w)
      X.so.coefficients <- qr.coef(X.so.qr, Y.to.C)   #beta
      X.to.C.w <- qr.resid(X.so.qr, X.to.C.w) #Z*
      Y.to.C <- qr.resid(X.so.qr, Y.to.C)     #Y*
    }

    if(standardize) {
      X.to.C.stds <- sqrt(apply(X.to.C.w, 2, var))
      if(any(X.to.C.stds < sqrt(.Machine$double.eps)))
        stop("Matrix build from transformed variables has a constant column")
      X.to.C.w <- sweep(X.to.C.w, 2, X.to.C.stds, "/")
    }

    Cvector.old <- Cvector

    fit <- .C("lasso",
              X = as.double(X.to.C.w),
              n = as.integer(nobs),
              p = as.integer(p),
              bound = as.double(bound),
              coefficients = as.double(guess.constrained.coefficients),
              Y = as.double(Y.to.C),
              fitted.values = double(nobs),
              residuals = double(nobs),
              Lagrangian = double(1),
              success = integer(1),
              trace = do.trace,
              assub = as.logical(FALSE),
              PACKAGE = "lasso2")[keep]

    if(fit$success < 0)
      stop("Uups, something went wrong in the C-routine")

    ## else
    fit$success <- NULL

    guess.constrained.coefficients <- gamma <- fit$coefficients

    if(something.to.sweep.out){
      Cvector <- c(X.so.coefficients, gamma)
      eta <- as.vector((X.to.C.w %*% gamma + X.sweep.out.w %*%
                        X.so.coefficients)/W)
    }else{
      Cvector <- c(gamma)
      eta <- as.vector((X.to.C.w %*% gamma)/W)
    }

    fit$fitted.values <- mu <- inv.link(eta + offset) ## Our new values of mu and eta

    j <- j + 1

    continue <-
        j < control$maxit &&
        sum((Cvector - Cvector.old)^2) > control$epsilon *sum(Cvector^2)

    if(do.trace) {
        cat(">>> while(continue = ", continue,") { .. } : j = ",j,
            "  summary(W) weights:\n")
        print(summary(W)) ; cat("\n")
    }

  } ## while(continue)

  if(j == control$maxit)
      warning("Maximal number of iterations reached.")

  fit$iter <- j
  deviance <- family$deviance(mu, y, weights)
  fit$linear.predictors <- eta

  fit$xtx <- crossprod(X.to.C)
  fit$xtr <- crossprod(X.to.C,fit$residuals)
  fit$constrained.coefficients <- fit$coefficients
  names(fit$coefficients) <- names(fit$constrained.coefficients) <-
    dimnames(X.to.C)[[2]]
  if(standardize) {
    fit$X.stds <- X.to.C.stds
    fit$coefficients <- fit$coefficients/X.to.C.stds
  }
  if(something.to.sweep.out) {
    X.so.X <- qr.coef(X.so.qr, X.to.C*W)
    X.so.coefficients <- X.so.coefficients - X.so.X %*% fit$coefficients
    tmp <- fit$coefficients
    fit$coefficients <- rep(0, ncol(X))
    names(fit$coefficients) <- dimnames(X)[[2]]
    fit$coefficients[X.names] <- tmp
    names(X.so.coefficients) <- X.sweep.out.names
    ind <- names(fit$coefficients[name.matches])
    fit$coefficients[name.matches] <- X.so.coefficients[ind]#fixed "coeff" !
    fit$sweep.out <-
        list(sweep.out = sweep.out, X.sweep.out.names = X.sweep.out.names,
             name.matches = name.matches,
             all.matched = all.matched, some.matched = some.matched,
             X.so.rtr.inv = qr.rtr.inv(X.so.qr), X.so.X = X.so.X)
  }
  fit$weights <- W
  fit$prior.weights <- weights
  fit$terms <- Terms
  fit$call <- call
  fit$contrasts <- attr(X, "contrasts")
  fit$assign <- attr(X, "assign")
  fit$family <- family
  fit$deviance <- deviance
  if(ret.x) fit$x <- X
  if(ret.y) fit$y <- y
  fit$df.res <- nobs - p # must always have full rank here
  ##    if(!absolute.t)
  ##      fit$relative.bound <- relative.bound
  structure(fit, class = c("gl1ce", "l1ce"),
            na.message = attr(mf, "na.message"))
}

## trivial accessor function:
family.gl1ce <- function (object, ...) object$family
