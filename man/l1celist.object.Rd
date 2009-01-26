% Copyright (C) 1998
% Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
% Bill Venables <wvenable@stats.adelaide.edu.au>
% $Id: l1celist.object.d,v 1.2 1998/09/23 02:55:04 bturlach Exp $
% --> ../COPYRIGHT for more details
\name{l1celist.object}
\alias{l1celist.object}
\title{Object of Several L1 Constrained Estimation Models}
\description{
These are objects of class \code{"l1celist"}
They represent the fits of several regression models under an L1
constraint on (some of the parameters).
}
\section{GENERATION}{
This class of objects is returned from the \code{l1ce} function
to represent a fitted model.
}
\section{METHODS}{
The \code{"l1celist"} class of objects has methods for the following generic
functions:
\code{[}, \code{[[}, \code{aux}, \code{coef}, \code{deviance}, \code{fitted}, \code{formula}, \code{gcv},
\code{labels}, \code{plot}, \code{print}, \code{residuals}, \code{vcov}.
}
\section{STRUCTURE}{
  An object of class \code{l1celist} is a list of lists.  Each component of
  the list is a list with the information of an object of class \code{l1ce}
  that is unique for that information.  All shared information is stored
  as attributes.
  
  
  Each component of the list must include the following components
  for it to be a legitimate \code{l1celist} object.
  \describe{
    \item{coefficients}{
      the coefficients of the fit of the response to the
      columns of the model
      matrix.
      The names of the coefficients are the names of the
      columns of the model matrix.
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
    \item{bound}{
      the (absolute) L1 constraint imposed on the parameters.
    }
    \item{relative.bound}{
      optional, the (relative) L1 constraint imposed on the parameters.
      Present if \code{absolute.t=F}.
    }
    \item{Lagrangian}{
      the value of the Lagrangian that enforces the constraint at the
      solution.
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
  }
}
\details{
The residuals, fitted values, coefficients, and effects should be extracted
by the generic functions of the same name, rather than
by using the \code{[[]]} and the \code{$} operator.
}
\seealso{
\code{\link{l1ce}}, \code{\link{coefficients}}.
}
\keyword{classes}
\keyword{regression}
\keyword{methods}
% Converted by Sd2Rd version 1.21.
