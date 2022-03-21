daarem_base_noobjfn <- function(par, fixptfn, maxiter, tol, resid.tol, 
                                a1, kappa, num.params, nlag, intermed, ...) {
  
  Fdiff <- Xdiff <- matrix(0.0, nrow=num.params, ncol=nlag)
  rho <- resid.tol ## should this be user specified?
  
  resid_vals <- rep(NA, maxiter + 2)
  xold <- par
  xnew <- fixptfn(xold, ...)
  fold <- xnew - xold
  resid_vals[1] <- sqrt(crossprod(fold))
  fnew <- fixptfn(xnew, ...) - xnew
  resid_vals[2] <- sqrt(crossprod(fnew))
  ss.resids <- resid_vals[2]
  
  if(intermed) {
    p.inter <- rbind(par, xnew)
  } else {
    p.inter <- NULL
  }
  
  fp.evals <- 2
  k <- 1
  n.aa <- 0
  count <- 0
  shrink.count <- 0
  shrink.target <- 1/(1 + a1^kappa)
  lambda.ridge <- 100000
  r.penalty <- 0
  conv <- TRUE
  while(fp.evals < maxiter) {
     count <- count + 1
    
     Fdiff[,count] <- fnew - fold
     Xdiff[,count] <- xnew - xold
    
     np <- count
     if(np==1) {
        Ftmp <- matrix(Fdiff[,1], nrow=num.params, ncol=np)
        Xtmp <- matrix(Xdiff[,1], nrow=num.params, ncol=np)  
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
    
     xbar <- xnew - drop(Xtmp%*%gamma_vec)
     fbar <- fnew - drop(Ftmp%*%gamma_vec)
     x.propose <- xbar + fbar
    
     ftmp <- try(fixptfn(x.propose, ...)) 
     if(class(ftmp)[1] == "try-error") {
       break
     } else {
       f.propose <- ftmp - x.propose
     }
     #f.propose <- fixptfn(x.propose, ...) - x.propose
     fp.evals <- fp.evals + 1
     ss.propose <- sqrt(crossprod(f.propose))
     if(ss.propose <= ss.resids*(1.00 + rho^k)) {  
         ## Increase delta
         fold <- fnew
         xold <- xnew
         xnew <- x.propose
         fnew <- f.propose
          
         shrink.count <- shrink.count + 1
         ss.resids <- ss.propose
         n.aa <- n.aa + 1
     } else {
         ## Keep delta the same
         fold <- fnew
         xold <- xnew
         xnew <- fold + xold
          
         fnew <- fixptfn(xnew, ...) - xnew
         ss.resids <- sqrt(crossprod(fnew))
         fp.evals <- fp.evals + 1
     }

     resid_vals[k + 2] <- ss.resids
     if(ss.resids < tol & count==nlag) break
    
     if(count==nlag) {
        count <- 0
        ## restart count
     }
     shrink.target <-  1/(1 + a1^(kappa - shrink.count))
     if(intermed) {
       p.inter <- rbind(p.inter, xnew)
     }
     k <- k+1
  }
  if(fp.evals >= maxiter) {
    conv <- FALSE
  }
  return(list(par=c(xnew), fpevals = fp.evals, value.objfn=NULL, objfevals=NULL, convergence=conv, objfn.track=NULL,
              residuals=resid_vals[!is.na(resid_vals)], p.intermed=p.inter))
}
