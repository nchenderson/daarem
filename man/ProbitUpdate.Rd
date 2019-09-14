\name{ProbitUpdate}
\alias{ProbitUpdate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EM Algorithm Update for Probit Regression
}
\description{
Function performs an EM update (both the E and M steps) of the parameters for a Probit regression model.
}
\usage{
ProbitUpdate(beta.hat, X, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta.hat}{
A vector of length p. The current estimates of the regression parameters.
}
\item{X}{
The n x p design matrix for the Probit regression model.
}
\item{y}{
Vector of length n containing binary outcomes (either 0 or 1).
}
}
\value{
A vector of length p - the updated parameter values.
}
\author{
Nicholas Henderson
}
\seealso{
\code{\link{ProbitSimulate}}, \code{\link{ProbitLogLik}}
}
\examples{
n <- 200
npars <- 5
true.beta <- .5*rt(npars, df=2) + 2
XX <- matrix(rnorm(n*npars), nrow=n, ncol=npars)
yy <- ProbitSimulate(true.beta, XX)

initial.beta <- rep(0.0, npars)
new.beta <- ProbitUpdate(initial.beta, XX, yy)
}
\keyword{ probit }% use one of  RShowDoc("KEYWORDS")
\keyword{ regression }% __ONLY ONE__ keyword per line
