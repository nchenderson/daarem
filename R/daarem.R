daarem <- function(par, fixptfn, objfn, ..., control=list()) {

  control.default <- list(maxiter=2000, order=5, tol=1.e-08, mon.tol=0.01, cycl.mon.tol=0.0, kappa=25, 
                          alpha=1.2, resid.tol=0.95, convtype="param", intermed=FALSE)
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
  if(control$convtype=="param") {
      check.par.resid <- TRUE
  } else if(control$convtype=="objfn") {
      check.par.resid <- FALSE
  }
  intermed <- control$intermed
  
  num.params <- length(par)
  nlag <- min(control$order, ceiling(num.params/2))

 
  if(!missing(objfn)) {
      ans <- daarem_base_objfn(par=par, fixptfn=fixptfn, objfn=objfn, maxiter=maxiter, 
                               tol=tol, mon.tol=mon.tol, cycl.mon.tol=cycl.mon.tol, a1=a1, 
                               kappa=kappa, num.params=num.params, nlag=nlag, 
                               check.par.resid=check.par.resid, intermed=intermed, ...) 
  } else {
      ans <- daarem_base_noobjfn(par=par, fixptfn=fixptfn, maxiter=maxiter, tol=tol, 
                                 resid.tol=resid.tol, a1=a1, kappa=kappa, 
                                 num.params=num.params, nlag=nlag, intermed=intermed, ...) 
  }
  if(!ans$convergence) {
     warning("Algorithm did not converge")
  }
  return(ans)
}
