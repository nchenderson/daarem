\name{ProbitSimulate}
\alias{ProbitSimulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Data from a Probit Regression Model
}
\description{
Function to simulate data from a Probit regression model. User provides a design matrix and a vector of regression coefficients.
Output is a vector of 0/1 responses.
}
\usage{
ProbitSimulate(beta.vec, X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta.vec}{
A vector of length p containing the true regression coefficients of the Probit regression model to be simulated from.
}
\item{X}{
An n x p design matrix for the Probit regression model to be simulated from.
}
}
\value{
A vector of length n containing binary outcomes (i.e., 0 or 1).
}
\author{
Nicholas Henderson
}

\seealso{
\code{\link{ProbitUpdate}}, \code{\link{ProbitLogLik}}
}
\examples{
n <- 200
npars <- 5
true.beta <- .5*rt(npars, df=2) + 2
XX <- matrix(rnorm(n*npars), nrow=n, ncol=npars)
yy <- ProbitSimulate(true.beta, XX)
}
\keyword{probit}% use one of  RShowDoc("KEYWORDS")
\keyword{regression}% __ONLY ONE__ keyword per line
