% Copyright (C) 1998
% Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
% Bill Venables <wvenable@stats.adelaide.edu.au>
% $Id: Prostate.d,v 1.4 1999/12/01 05:07:34 bturlach Exp $
% --> ../COPYRIGHT for more details
\name{Prostate}
\alias{Prostate}
\title{Prostate Cancer Data}
\description{
  These data come from a study that examined the correlation between the
  level of prostate specific antigen and a number of clinical measures
  in men who were about to receive a radical prostatectomy.  It is data
  frame with 97 rows and 9 columns.
}
\usage{
data(Prostate)
}
\format{
  The data frame has the following components:
  \describe{
    \item{\code{lcavol}}{log(cancer volume)}
    \item{\code{lweight}}{log(prostate weight)}
    \item{\code{age}}{age}
    \item{\code{lbph}}{log(benign prostatic hyperplasia amount)}
    \item{\code{svi}}{seminal vesicle invasion}
    \item{\code{lcp}}{log(capsular penetration)}
    \item{\code{gleason}}{Gleason score}
    \item{\code{pgg45}}{percentage Gleason scores 4 or 5}
    \item{\code{lpsa}}{log(prostate specific antigen)}
  }
}
\source{
  Stamey, T.A., Kabalin, J.N., McNeal, J.E., Johnstone, I.M., Freiha,
  F., Redwine, E.A. and Yang, N. (1989)\cr
  Prostate specific antigen in the diagnosis and treatment of
  adenocarcinoma of the prostate: II. radical prostatectomy treated
  patients,
  \emph{Journal of Urology} \bold{141}(5), 1076--1083.
}
\examples{
data(Prostate)
attach(Prostate)
pairs(Prostate, col = 1+svi, pch = gleason - 5,
      main = paste("Prostate data, n = ", nrow(Prostate)))
detach()

l1c.P <- l1ce(lcavol ~ ., data = Prostate)
coef(l1c.P)[coef(l1c.P) != 0] ## only age, lcp, lpsa (+ intercept)
summary(l1c.P)
}
\keyword{datasets}
% Converted by Sd2Rd version 1.21.
