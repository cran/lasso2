\name{qr.rtr.inv}
\alias{qr.rtr.inv}
\title{Reconstruct the Inverse of R'R from a QR Object}
\description{
  From a QR object, compute the inverse matrix which is implicitely (but
  not explicitly!) used to solve the underlying least squares problem.
}
\usage{
qr.rtr.inv(qr)
}
\arguments{
  \item{qr}{\"qr\" object, typically resulting from \code{\link{qr}(.)}.}
}
\value{
  The \eqn{p \times p}{p * p} matrix \eqn{(R'R)^{-1}}{(R'R)^(-1)}
  or equivalently, the inverse of \eqn{X'X} (i.e. \code{t(X) \%*\% X} in \R).
}
\seealso{\code{\link{qr}}, \code{\link{qr.R}}, \code{\link{backsolve}}.}
\examples{
(h3 <- 1/outer(0:5, 1:3, "+"))
rtr <- qr.rtr.inv(qr(h3))
all.equal(c(rtr \%*\% 1:3), solve(crossprod(h3), 1:3))
}
\keyword{algebra}
\keyword{array}
