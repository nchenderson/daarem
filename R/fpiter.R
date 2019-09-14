fpiter <- function(par, fixptfn, objfn=NULL, control=list( ), ...){

  control.default <- list(tol=1.e-07, maxiter=5000, trace=FALSE)
  namc <- names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  ctrl <- modifyList(control.default, control)

  #
  # method = reduced-rank ("rre") or minimal-polynomial ("mpe") extrapolation
  # K = order of extrapolation scheme; K=2,3,4 are typical choices
  # square = a logical variable indicating whether or not "squaring" is used

  tol <- ctrl$tol
  maxiter <- ctrl$maxiter
  trace <- ctrl$trace

  if (trace) cat("fpiter \n")

  iter <- 1
  resid <- rep(NA,1)
  objeval <- 0
  conv <- FALSE

  while (iter < maxiter) {

    p.new <- fixptfn(par, ...)
    res <- sqrt(crossprod(p.new - par))

    if ( res < tol) {conv <- TRUE; break}

    if (trace) {
      if (!is.null(objfn)) {cat("Iter: ", iter, "Objective fn: ",objfn(par, ...), "\n"); objeval <- objeval + 1}
      else cat("Iter: ", iter, "Residual: ",res, "\n")
    }
    par <- p.new
    iter <- iter+1
  }

  loglik.best <-  if (!is.null(objfn)) objfn(par, ...) else NA

  return(list(par=par, value.objfn=loglik.best, fpevals=iter, objfevals = objeval, convergence=conv))
}
