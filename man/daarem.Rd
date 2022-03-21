\name{daarem}
\alias{daarem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Damped Anderson Acceleration with Restarts and Epsilon-Montonicity for Accelerating Slowly-Convergent, Monotone
Fixed-Point Iterations
}
\description{
An `off-the-shelf' acceleration scheme for accelerating the convergence of \emph{any} smooth, monotone, slowly-converging
fixed-point iteration. It can be used to accelerate the convergence of a wide
variety of montone iterations including, for example, expectation-maximization (EM) algorithms
and majorization-minimization (MM) algorithms.
}
\usage{
daarem(par, fixptfn, objfn, ..., control=list())
}
\arguments{
     \item{par}{A vector of starting values of the parameters.}

    \item{fixptfn}{A vector function, \eqn{G}{G} that denotes the fixed-point
      mapping.  This function is the most essential input in the package.
      It should accept a parameter vector as input and should return a
      parameter vector of the same length. This function defines the fixed-point
      iteration: \eqn{x_{k+1} = G(x_k)}{x[k+1] = G(x[k])}.
      In the case of an EM algorithm, \eqn{G}{G} defines a single E and M step.}

    \item{objfn}{This is a scalar function, \eqn{L}{L}, that denotes
      a ''merit'' function which attains its local maximum at the fixed-point
      of \eqn{G}{G}. The function \eqn{L}{L} should accept a parameter vector as
      input and should return a scalar value.  In the EM algorithm, the merit
      function \eqn{L}{L} is the log-likelihood function. It is not necessary
      for the user to provide this argument though it is preferable.}

    \item{control}{A list of control parameters specifying any changes to
      default values of algorithm control parameters. Full names of control list elements must be specified, otherwise, user-specifications are ignored.
	See *Details*.}
    \item{...}{Arguments passed to \code{fixptfn} and  \code{objfn}.}
}
\value{
A list with the following components:
\item{par}{Parameter, \eqn{x*}{x*} that are the fixed-point of \eqn{G}{G}
  such that \eqn{x* = G(x*)}{x* = G(x*)}, if convergence is successful.}

\item{value.objfn}{The value of the objective function \eqn{L}{L} at termination.}

\item{fpevals}{Number of times the fixed-point function \code{fixptfn} was evaluated.}

\item{objfevals}{Number of times the objective function \code{objfn} was evaluated.}

\item{convergence}{An integer code indicating type of convergence.  \code{0}
  indicates successful convergence, whereas \code{1} denotes failure to
  converge.}
\item{objfn.track}{A vector containing the value of the objective function at each iteration.}
\item{p.intermed}{A matrix where each row corresponds to parameters at each iteration, 
      along with the corresponding value of the objective function (in the first column). 
      This object is returned only when the control parameter \code{intermed} is set to \code{TRUE}. It is not returned when 
      \code{objfn} is not specified.}
}
\details{
Default values of \code{control} are:
    \code{maxiter=2000},
    \code{order=10},
    \code{tol=1e-08},
    \code{mon.tol=0.01},
    \code{cycl.mon.tol=0.0},
    \code{alpha=1.2},
    \code{kappa=25},
    \code{resid.tol=0.95},
    \code{convtype="param"}

\describe{
	\item{\code{maxiter}}{An integer denoting the maximum limit on the number
	  of evaluations of  \code{fixptfn}, \eqn{G}{G}.  Default value is 2000.}

\item{\code{order}}{An integer \eqn{1}{>= 1} denoting the order of the DAAREM acceleration scheme.}
\describe{
	\item{\code{tol}}{A small, positive scalar that determines when iterations
	  should be terminated.  When \code{convtype} is set to "param", iteration is terminated when
	  \eqn{||x_k - G(x_k)|| < tol}{|| x[k] - G(x[k]) || < tol}.
	  Default is \code{1.e-08}.}
	\item{\code{mon.tol}}{A nonnegative scalar that determines whether the montonicity condition
	is violated. The monotonicity condition is violated whenver \eqn{L(x[k+1]) < L(x[k]) - mon.tol }{L(x[k+1]) < L(x[k]) - mon.tol }.
	Such violations determine how much damping is to be applied on subsequent steps of the algorithm. Default
	value of mon.tol is \code{1.e-02}.}
	\item{\code{cycl.mon.tol}}{A nonegative scalar that determines whether a montonicity condition is violated
	after the end of the cycle. This cycle-level monotonicity condition is violated whenver
	\eqn{L(x[end cycle]) < L(x[start cycle]) - cycl.mon.tol }{L(x[end cycle]) < L(x[start cycle]) - cycl.mon.tol }.
	Here, \eqn{x[start cycle]}{x[start cycle]} refers to the value of \eqn{x}{x} at the beginning of the current cycle
	while \eqn{x[end cycle]}{x[end cycle]} refers to the value of \eqn{x}{x} at the end of the current cycle.
	Such violations also determine how much damping is to be applied on subsequent steps of the algorithm.}
	\item{\code{kappa}}{A nonnegative parameter which determines the ``half-life'' of relative damping and how quickly relative
	damping tends to one. In the absence of monotonicity
	violations, the relative damping factor is \eqn{<= 1/2}{<= 1/2} for the first \code{kappa} iterations, and it is
	then greater than \eqn{1/2}{1/2} for all subsequent iterations. The relative damping factor is the ratio between
	the norm of the unconstrained coefficients in Anderson acceleration and the norm of the damped coefficients.
	In the absence of any monotonicity violations, the relative damping factor in iteration \eqn{k}{k} is
	\eqn{1/(1 + \alpha^(\kappa - k))}{1/(1 + \alpha^(\kappa - k))}. 
	}
	\item{\code{alpha}}{A parameter \eqn{> 1}{> 1} that determines the initial relative damping factor and how quickly
	the relative damping factor tends to one. The initial relative damping factor is \eqn{1/(1 + \alpha^\kappa)}{1/(1 + \alpha^\kappa)}.
	In the absence of any monotonicity violations, the relative damping factor in iteration \eqn{k}{k} is
	\eqn{1/(1 + \alpha^(\kappa - k))}{1/(1 + \alpha^(\kappa - k))}.}
	}
	\item{\code{resid.tol}}{A nonnegative scalar \eqn{< 1}{< 1} that determines whether a residual change condition is violated.
	The residual change condition is violated whenever \eqn{||x_k+1 - G(x_k+1)|| > ||x_k - G(x_k)|| (1 + resid.tol^k)}{|| G(x[k+1]) - x[k+1] || > || G(x[k]) - x[k] ||*(1 + resid.tol^k)}. Default value of resid.tol is \code{0.95}.
	}
	\item{\code{convtype}}{This can equal either "param" or "objfn". When set to "param", convergence is determined by the criterion: \eqn{||x_k - G(x_k)|| \leq tol}{|| x[k] - G(x[k]) || < tol}. When set to "objfn", convergence is determined by the objective function-based criterion: \eqn{| L(x[k+1]) - L(x[k])| < tol }{| L(x[k+1]) - L(x[k])| < tol }.
	}
	\item{\code{intermed}}{
	A logical variable indicating whether or not the intermediate results of iterations should be returned. If set to
	\code{TRUE}, the function will return a matrix where each row corresponds to parameters at each iteration, along with the corresponding value of the objective function in the first column. This option is inactive when objfn is not specified. Default is \code{FALSE}.
	}
}
}
\references{
Henderson, N.C. and Varadhan, R. (2019) \emph{Damped Anderson acceleration with restarts and monotonicity
control for accelerating EM and EM-like algorithms}, \emph{Journal of Computational and Graphical Statistics, Vol. 28(4), 834-846.} \url{https://doi.org/10.1080/10618600.2019.1594835}
}
\author{
Nicholas Henderson and Ravi Varadhan
}
\seealso{
\code{\link{fpiter}}
}
\examples{

n <- 2000
npars <- 25
true.beta <- .5*rt(npars, df=2) + 2
XX <- matrix(rnorm(n*npars), nrow=n, ncol=npars)
yy <- ProbitSimulate(true.beta, XX)
max.iter <- 1000
beta.init <- rep(0.0, npars)

# Estimating Probit model with DAAREM acceleration
aa.probit <- daarem(par=beta.init, fixptfn = ProbitUpdate, objfn = ProbitLogLik,
                    X=XX, y=yy, control=list(maxiter=max.iter))

plot(aa.probit$objfn, type="b", xlab="Iterations", ylab="log-likelihood")

\donttest{
# Compare with estimating Probit model using the EM algorithm

max.iter <- 25000  # need more iterations for EM convergence
beta.init <- rep(0.0, npars)

em.probit <- fpiter(par=beta.init, fixptfn = ProbitUpdate, objfn = ProbitLogLik,
                    X=XX, y=yy, control=list(maxiter=max.iter))
c(aa.probit$fpevals, em.probit$fpevals)
c(aa.probit$value, em.probit$value)
}

# Accelerating using SQUAREM if the SQUAREM package is loaded
# library(SQUAREM)
# max.iter <- 5000
# sq.probit <- squarem(par=beta.init, fixptfn = ProbitUpdate, objfn = ProbitLogLik,
#                     X=XX, y=yy, control=list(maxiter=max.iter))

# print( c(aa.probit$fpevals, em.probit$fpevals, sq.probit$fpevals) )
# print( c(aa.probit$value, em.probit$value, sq.probit$value) )
# print( c(aa.probit$objfeval, em.probit$objfeval, sq.probit$objfeval) )

}
\keyword{ optimization }
\keyword{ EM algorithm }
