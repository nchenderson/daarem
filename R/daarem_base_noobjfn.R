daarem_base_noobjfn <- function(par, fixptfn, maxiter, tol, mon.tol, 
                               cycl.mon.tol, a1, kappa, num.params, nlag, ...) {
  
  
  Fdiff <- Xdiff <- matrix(0.0, nrow=num.params, ncol=nlag)
  obj_funvals <- rep(NA, maxiter + 2)
  
  xold <- par
  xnew <- fixptfn(xold, ...)
  #obj_funvals[1] <- objfn(xold, ...)
 # obj_funvals[2] <- objfn(xnew, ...)
  
  fold <- xnew - xold
  k <- 1
  count <- 0
  shrink.count <- 0
  shrink.target <- 1/(1 + a1^kappa)
  lambda.ridge <- 100000
  r.penalty <- 0
  conv <- TRUE
  new.objective.val <- 0
  while(k < maxiter) {
    count <- count + 1
    
    fnew <- fixptfn(xnew, ...) - xnew
    ss.resids <- sqrt(crossprod(fnew))
    if(ss.resids < tol & count==nlag) break
    
    Fdiff[,count] <- fnew - fold
    Xdiff[,count] <- xnew - xold
    
    np <- count
    if(np==1) {
      Ftmp <- matrix(Fdiff[,1], nrow=length(fnew), ncol=np)
      Xtmp <- matrix(Xdiff[,1], nrow=length(fnew), ncol=np)  
    } else {
      Ftmp <- Fdiff[,1:np]
      Xtmp <- Xdiff[,1:np]  
    }
    tmp <- La.svd(Ftmp)
    dvec <- tmp$d
    dvec.sq <- dvec*dvec
    uy <- crossprod(tmp$u, fnew)
    uy.sq <- uy*uy
    
    ### Still need to compute Ftf
    Ftf <- sqrt(sum(uy.sq*dvec.sq))
    tmp_lam <- DampingFind(uy.sq, dvec, a1, kappa, shrink.count, Ftf, lambda.start=lambda.ridge, r.start=r.penalty)
    lambda.ridge <- tmp_lam$lambda
    r.penalty <- tmp_lam$rr
    dd <- (dvec*uy)/(dvec.sq + lambda.ridge)
    gamma_vec <- crossprod(tmp$vt, dd)
    
    if(class(gamma_vec) != "try-error"){
      
      xbar <- xnew - drop(Xtmp%*%gamma_vec)
      fbar <- fnew - drop(Ftmp%*%gamma_vec)
      
      x.propose <- xbar + fbar
      #new.objective.val <- try(objfn(x.propose, ...), silent=TRUE)

      if(class(new.objective.val) != "try-error") {
        if(TRUE) {  ## just change this line in daarem_base_noobjfn
          ## Increase delta
          fold <- fnew
          xold <- xnew
          
          xnew <- x.propose
          shrink.count <- shrink.count + 1
        } else {
          ## Keep delta the same
          fold <- fnew
          xold <- xnew
          
          xnew <- fold + xold
        }
      } else {
        ## Keep delta the same
        fold <- fnew
        xold <- xnew
        
        xnew <- fold + xold
        count <- 0
      }
    } else {
      ## Keep delta the same
      fold <- fnew
      xold <- xnew
      
      xnew <- fold + xold
      count <- 0
    }
    if(count==nlag) {
      count <- 0
      ## restart count
      ## make comparison here l.star vs. obj_funvals[k+2]
      if(FALSE) {
        ## Decrease delta
        shrink.count <- max(shrink.count - nlag, -2*kappa)
      }
      #ell.star <- obj_funvals[k+2]
    }
    
    shrink.target <-  1/(1 + a1^(kappa - shrink.count))
    k <- k+1
  }
  if(k >= maxiter) {
    conv <- FALSE
  }
  return(list(par=c(xnew), fpevals = k, value.objfn=NULL, objfevals=NULL, convergence=conv, objfn.track=NULL))
}
