\name{fpiter}
\alias{fpiter}
\title{Fixed-Point Iteration Scheme}
\description{A function to implement the fixed-point iteration algorithm.  This includes monotone, contraction mappings including EM and MM algorithms}
\usage{
  fpiter(par, fixptfn, objfn=NULL, control=list( ), ...)
}
\arguments{
    \item{par}{A vector of parameters denoting the initial guess for the
       fixed-point iteration.}

    \item{fixptfn}{A vector function, \eqn{F}{F} that denotes the fixed-point
       mapping. This function is the most essential input in the package.  It
       should accept a parameter vector as input and should return a parameter
       vector of same length. This function defines the fixed-point iteration:
       \eqn{x_{k+1} = F(x_k)}{x[k+1] = F(x[k]}.
       In the case of EM algorithm, \eqn{F}{F} defines a single E and M step.}

    \item{objfn}{This is a scalar function, $L$, that denotes a ''merit''
       function which attains its local minimum at the fixed-point of $F$.
       This function should accept a parameter vector as input and should
       return a scalar value. In the EM algorithm, the merit function \eqn{L}{L}
       is the log-likelihood. In some problems, a natural merit function may
       not exist, in which case the algorithm works with only \code{fixptfn}.
       The merit function function \code{objfn} does not have to be specified,
       even when a natural merit function is available, especially when its
       computation is expensive.}

    \item{control}{ A list of control parameters to pass on to the algorithm.  Full names of control list elements must be specified, otherwise, user-specifications are ignored.  See *Details* below.}

    \item{...}{Arguments passed to \code{fixptfn} and \code{objfn}.}
}

\value{A list with the following components:
\item{par}{Parameter,\eqn{x*}{x*} that are the fixed-point of \eqn{F}{F} such
   that \eqn{x* = F(x*)}{x* = F(x*)}, if convergence is successful.}

\item{value.objfn}{The value of the objective function \eqn{L}{L} at termination.}

\item{fpevals}{Number of times the fixed-point function \code{fixptfn} was evaluated.}

\item{objfevals}{Number of times the objective function \code{objfn} was evaluated.}

\item{convergence}{An integer code indicating type of convergence.
  \code{0} indicates successful convergence,
  whereas \code{1} denotes failure to converge.}
}

\details{
\code{control} is list of control parameters for the algorithm.

\describe{
\code{control = list(tol = 1.e-07, maxiter = 1500, trace = FALSE)}

    \code{tol}{  A small, positive scalar that determines when iterations
       should be terminated.  Iteration is terminated when
       \eqn{||x_k - F(x_k)|| \leq tol}{abs(x[k] - F(x[k]) <= tol}.
       Default is \code{1.e-07}.}

    \code{maxiter}{  An integer denoting the maximum limit on the number of
       evaluations of  \code{fixptfn}, \eqn{F}{F}.  Default is \code{1500}.}

    \code{trace}{  A logical variable denoting whether some of the intermediate
       results of iterations should be displayed to the user.
       Default is \code{FALSE}.}
  }
}

\seealso{
  \code{\link{daarem}}
   }

\examples{

### Generate outcomes from a probit regression model
n <- 1000
npars <- 5
true.beta <- .5*rt(npars, df=2) + 1
XX <- matrix(rnorm(n*npars), nrow=n, ncol=npars)
yy <- ProbitSimulate(true.beta, XX)
max.iter <- 1000
beta.init <- rep(0.0, npars)

### EM algorithm for estimating parameters from probit regression

em.probit <- fpiter(par=beta.init, fixptfn = ProbitUpdate, X=XX, y=yy,
                    control=list(maxiter=max.iter))
}
\keyword{optimization}
\keyword{EM algorithm}

