###  Copyright (C) 1998
###  Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
###  Bill Venables <wvenable@stats.adelaide.edu.au>
### --> ../COPYRIGHT for more details

###--  vcov() is in MASS as well -- this generic is compatible :
vcov <- function(object, ...) UseMethod("vcov")

vcov.l1ce <- function(object, type = c("OPT", "Tibshirani"),
                      gen.inverse.diag = 0, ...)
{
  type <- match.arg(type)

  xtx <- object$xtx
  cov.mat <- switch(type,
                    OPT = opt.covmat(xtx, object$xtr, object$bound,
                      object$Lagrangian),
                    Tibshirani = tib.covmat(xtx, gen.inverse.diag,
                      object$constrained.coefficients, object$Lagrangian)
                    )
  cov.mat <- solve(cov.mat)
  df <- tr(junk <- cov.mat %*% xtx)
  cov.mat <- junk %*% cov.mat

  if(!is.null(object$X.stds)) {
    X.stds <- as.vector(object$X.stds)
    p <- length(X.stds)
    cov.mat <- cov.mat/X.stds
    cov.mat <- cov.mat/X.stds[rep(1:p,rep(p,p))]
  }

  if(!is.null(object$sweep.out)) {
    if(!object$sweep.out$all.matched)
      warning("Variables in sweep.out were not a subset of variables in model.  Results could be meaningless.")

    psi1 <- object$sweep.out$X.so.rtr.inv
    psi2 <- cov.mat
    psi12 <- object$sweep.out$X.so.X %*% psi2
    psi1 <- psi1 + psi12 %*% t(object$sweep.out$X.so.X)

    p1 <- nrow(psi1)
    p2 <- nrow(psi2)
    ## 2 x 2 block matrix -- wouldn't rbind(cbind(),cbind()) be faster? FIXME
    cov.mat <- matrix(0,p1+p2,p1+p2)
    cov.mat[1:p1,1:p1] <- psi1
    cov.mat[(p1+1):(p1+p2),(p1+1):(p1+p2)] <- psi2
    cov.mat[(p1+1):(p1+p2),1:p1] <- -t(psi12)
    cov.mat[1:p1,(p1+1):(p1+p2)] <- -psi12

    names2 <- names(object$coefficients)
    if(paste(R.version$major, R.version$minor, sep=".") < "1.4.1") {
        ## in R versions <= 1.4.0, most qr.* didn't have dimnames
        names1 <- object$sweep.out$X.sweep.out.names
        names1 <- c(names1, names2[! names2 %in% names1])
    }
    else ## X.so.X has dimnames in R version >= 1.4.1
        names1 <- do.call("c",dimnames(object$sweep.out$X.so.X))
    dimnames(cov.mat) <- list(names1, names1)
    cov.mat <- cov.mat[names2,names2]
    df <- df+p1
  }
  list(cov.unscaled=cov.mat,
       df = c(df,length(object$fitted)-df))
}

vcov.l1celist <- function(object, type = c("OPT", "Tibshirani"),
                          gen.inverse.diag = 0, ...)
{
    type <- match.arg(type)

    res <- list()
    for(i in 1:length(object))
        res[[i]] <- vcov(object[[i]],type,gen.inverse.diag)
    res
}

opt.covmat <- function(xtx, xtr, bound, Lagrangian)
{
  if( bound == 0 )
    stop("\"bound\" must be unequal zero when calculating covariance matrix")
  if( abs(Lagrangian) < sqrt(.Machine$double.eps) ){
    xtx
  }else{
    xtx + xtr%*%t(xtr)/(bound*Lagrangian)
  }
}

tib.covmat <- function(xtx, gen.inverse.diag, beta, Lagrangian)
{
  ind <- abs(beta) > sqrt(.Machine$double.eps)
  beta[ind]  <- 1/abs(beta[ind])
  beta[!ind] <- gen.inverse.diag

  xtx + Lagrangian*diag(beta)
}
