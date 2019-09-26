daarem <- function(par, fixptfn, objfn, ..., control=list()) {

  control.default <- list(maxiter=2000, order=5, tol=1.e-08, mon.tol=0.01, cycl.mon.tol=0.0, kappa=25, alpha=1.2, resid.tol=0.95)
  namc <- names(control)
  if (!all(namc %in% names(control.default))) {
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  }
  control <- modifyList(control.default, control)

  maxiter <- control$maxiter
  tol <- control$tol
  mon.tol <- control$mon.tol  ## monotonicity tolerance
  cycl.mon.tol <- control$cycl.mon.tol
  a1 <- control$alpha
  kappa <- control$kappa
  resid.tol <- control$resid.tol

  num.params <- length(par)
  nlag <- min(control$order, ceiling(num.params/2))

  if(!missing(objfn)) {
      ans <- daarem_base_objfn(par, fixptfn, objfn, maxiter, tol, mon.tol, 
                               cycl.mon.tol, a1, kappa, num.params, nlag, ...) 
  } else {
      ans <- daarem_base_noobjfn(par, fixptfn, maxiter, tol, resid.tol, 
                                 a1, kappa, num.params, nlag, ...) 
  }
  if(!ans$convergence) {
     warning("Algorithm did not converge")
  }
  return(ans)
}
