% Copyright (C) 1998
% Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
% Bill Venables <wvenable@stats.adelaide.edu.au>
% $Id: l1ce.object,v 1.2 1998/09/20 09:28:27 bturlach Exp $
% --> ../COPYRIGHT for more details
\name{l1ce.object}
\alias{l1ce.object}
\title{L1 Constrainted Estimation Model Object}
\description{
  These are objects of class \code{"l1ce"}.
  They represent the fit of a regression model under an L1 constraint on
  (some of) the parameters.
}
\section{GENERATION}{
  This class of objects is returned from the \code{\link{l1ce}} function
  to represent a fitted model.
}
\section{METHODS}{
  The \code{"l1ce"} class of objects has methods for the following generic
  functions:\cr
  \code{coef}, \code{deviance}, \code{fitted}, \code{formula},
  \code{gcv}, \code{labels}, \code{predict}, \code{print},
  \code{residuals}, \code{summary}, \code{vcov}.
}
\section{STRUCTURE}{
  The following components must be included in a legitimate \code{l1ce} object.
  \describe{
    \item{coefficients}{
      the coefficients of the fit of the response to the columns of the model
      matrix.  The names of the coefficients are the names of the columns
      of the model matrix.
    }
    \item{residuals}{
      the residuals from the fit.
      If weights were used, then the residuals are the raw residuals - the
      weights are not taken into account.
      If you need residuals that all have the same variance, then use the
      \code{residuals} function with \code{type="pearson"}.
    }
    \item{fitted.values}{
      the fitted values from the fit.
      If weights were used, the fitted values are not adjusted for the weights.
    }
    \item{bound}{the (absolute) L1 constraint imposed on the parameters.    }
    \item{relative.bound}{
      optional, the (relative) L1 constraint imposed on the parameters.
      Present if \code{absolute.t=F}.
    }
    \item{Lagrangian}{
      the value of the Lagrangian that enforces the constraint at the
      solution.
    }
    \item{xtx}{
      the moment matrix of the variables that are under the constraint.
      (After taking weights, sweep-out variables and standardization into
      account).
    }
    \item{xtr}{
      the product of the design matrix of the variables that are under the
      constraint (after taking weights, sweep-out variables and
      standardization into account) with the residual vector.
    }
    \item{constrained.coefficients}{
      the coefficients on the scale on which they are constrained.  Useful
      as initial value for further fits.
    }
    \item{sweep.out}{
      Optional information on the variables that are not under the
      constraint and on which the other variables and the response is
      projected first.%% FIXME?? isn't this needed by summary?
      Optional, not present if \code{sweep.out = NULL}.
    }
    \item{assign}{
      the list of assignments of coefficients (and effects) to the terms in
      the model.
      The names of this list are the names of the terms.
      The \code{i}th element of the list is the vector saying which coefficients
      correspond to the \code{i}th term.
      It may be of length 0 if there were no estimable effects for the term.
      See also \code{R.assign} below.
    }
    \item{terms}{%% FIXME : differs between R and S+
      an object of mode \code{expression} and class \code{term}
      summarizing the formula.  Used by various methods, but typically
      not of direct relevance to users.
    }
    \item{call}{
      an image of the call that produced the object, but with the arguments
      all named and with the actual formula included as the formula argument.
    }
    \item{contrasts}{
      a list containing sufficient information to construct the contrasts
      used to fit any factors occurring in the model.
      The list contains entries that are either matrices or character vectors.
      When a factor is coded by contrasts, the corresponding contrast matrix
      is stored in this list.
      Factors that appear only as dummy variables and variables in the model
      that are matrices correspond to character vectors in the list.
      The character vector has the level names for a factor or the column
      labels for a matrix.
    }
    \item{x}{optionally the model matrix, if \code{x=TRUE}.}
    \item{y}{optionally the response, if \code{y=TRUE}.}
  }%describe
}
\details{
  The residuals, fitted values, coefficients, and effects should be extracted
  by the generic functions of the same name, rather than
  by the \code{$} operator.
}
\seealso{
  \code{\link{l1ce}}, \code{\link{coefficients}}.
}
\keyword{classes}
\keyword{regression}
\keyword{methods}
% Converted by Sd2Rd version 1.21.
