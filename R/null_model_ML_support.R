#' @export
fit_generic_null<- function(data.vec,
                            family = c("htlnorm","zoibeta","cb"),
                            init = 0,
                            lower = -Inf,
                            upper = Inf,
                            method=c("Nelder-Mead","BFGS","L-BFGS-B","nlminb"),
                            id = NULL,
                            ...){
  method <- match.arg(method)
  family <- match.arg(family)

  data.vec <- .herb_data_check(data.vec, min.phi = 0, allow.zero = TRUE)

  if(family == "htlnorm"){
    optim.vars <- c("theta","meanlog","sdlog")
    param.val.trans = c(
      "theta" = function(x) {plogis(x)},
      "meanlog" = function(x) {x},
      "sdlog" = function(x) {exp(x)})
    nll <- function(theta){
      -sum(dhtlnorm(x = data.vec,
                    theta = plogis(theta[1]),
                    meanlog = theta[2],
                    sdlog = exp(theta[3]),
                    endpoint = 1,
                    log = TRUE))

    }
  } else if(family == "zoibeta"){
    optim.vars <- c("p0","p1","alpha","beta")
    param.val.trans = c(
      "p0" = function(x) {plogis(x)},
      "p1" = function(x) {plogis(x)},
      "alpha" = function(x) {exp(x)},
      "beta" = function(x) {exp(x)})
    nll <- function(theta){
      -sum(dzoibeta(x = data.vec,
                    p0 = plogis(theta[1]),
                    p1 = plogis(theta[2]),
                    alpha = exp(theta[3]),
                    beta = exp(theta[4]),
                    log = TRUE))

    }
  } else if(family == "cb"){
    optim.vars <- c("lambda")
    param.val.trans = c("lambda" = function(x) {plogis(x)})
    nll <- function(theta){
      -sum(dcb(x = data.vec,
                    lambda = plogis(theta[1]),
                    log = TRUE))

    }
  }


  if(length(lower) != length(optim.vars) ){
    lower <- rep(lower, length(optim.vars))
  }
  if(length(upper) != length(optim.vars)){
    upper <- rep(upper, length(optim.vars))
  }
  if(length(init) != length(optim.vars)){
    init <- rep(init, length(optim.vars))
  }

  if(method == "nlminb"){
    ML.fit<-nlminb(start = init,
                   objective = nll,
                   lower = lower,
                   upper = upper,
                   ...)
    hessian <- hessian(nll, ML.fit$par)
    loglik <- -ML.fit$objective
    iters <- ML.fit$evaluations

  } else {
    ML.fit<-optim2(par = init,
                  fn = nll,
                  hessian = T,
                  method = method,
                  lower=lower,
                  upper=upper,
                  ...)
    hessian <- ML.fit$hessian
    loglik <- -ML.fit$value
    iters <- ML.fit$counts
  }
  generic_null_fit.out<-list("theta.names" = optim.vars,
                          "par" = ML.fit$par,
                          "se" = tryCatch(sqrt(diag(solve(hessian))),
                                          error = function(e) {
                                            rep(NA,length(diag(hessian)))
                                          }),
                          "family" = family,
                          "loglik" = loglik,
                          "hessian" = hessian,
                          "iters" = iters,
                          "convergence" = ML.fit$convergence, # > 0 indicates issues
                          "method" = method,
                          "init" = init,
                          "data" = data.vec,
                          "param.val.trans" = param.val.trans,
                          "id" = id,
                          "n" = length(data.vec),
                          "df" = length(optim.vars)

  )
  if(any(is.na(generic_null_fit.out$se))){
    warning("Non-positive definite hessian matrix")
    generic_null_fit.out$convergence <- 2
  } else if(generic_null_fit.out$convergence!=0){
    warning("Model did not converge","  ",generic_null_fit.out$convergence)
  }

  generic_null_fit.out<-structure(generic_null_fit.out,
                               class=c("generic_null_fit","list"))
  return(generic_null_fit.out)
}

#' @title Print Values
#' @description Print 'generic_null_fit' objects
#' @param x a 'generic_null_fit' object
#' @param ... additional arguments
#' @param trans a function that is used to transform the estimated model coefficients. No transformation is applied when set to \code{NULL} (default).
#' @param digits number of digits to display. Default is 3.
#' @return a matrix
#' @export
print.generic_null_fit<-function(x, ..., trans = NULL, digits = 3){
  if(!inherits(x,"generic_null_fit")){
    stop("Needs to by of object type 'generic_null_fit'")
  }
  if(is.null(trans)){
    FUN <- function(x) x
  } else {
    FUN <- match.fun(trans)
  }

  cat("Fitted result on", x$id,"\n")
  cat("method =",x$method,"\n")
  cat("family =",x$family,"\n")
  cat("n =",x$n,"\n")
  cat("\n")
  print(c("loglik" = x$loglik, "AIC" = AIC(x), "AICc" = AICc(logLik(x)), "BIC" = BIC(x)))
  cat("\n")
  print(c(x$iters,"convergence"=x$convergence))
  cat("\n")
  print(noquote(vapply(
    x$param.val.trans[x$theta.names],
    function(z){
      gsub(".*\\{| |\\}","",paste0(deparse(z),collapse = ""))
    },
    character(1)
  )))
  cat("\n")
  out<-round(
    apply(as.matrix(
      data.frame(row.names = x$theta.names,
                 "Estimate" = x$par,
                 "Std." = x$se,
                 "lower" = x$par-x$se*1.96,
                 "upper" = x$par+x$se*1.96)
    ), 2, FUN = FUN),
    digits = digits)
  if(is.vector(out)){
    out<-t(as.matrix(out))
    rownames(out)<-x$theta.names
  }
  print(out)
}

#' @title Extract Log-Likelihood
#' @rdname logLik
#' @export
logLik.generic_null_fit <- function(object,...){
  out <- object$loglik
  out <- structure(out,
                   class= "logLik",
                   nall = object$n,
                   nobs= object$n,
                   df = object$df)
  return(out)
}

#' @title Akaike's An Information Criterion
#' @rdname AIC
#' @export
AIC.generic_null_fit <- function(object,...,k=2){
  dots <- as.list(match.call())[-1]
  dots <- dots[!names(dots) == "k"]

  if(length(dots) == 1){
    return(AIC(logLik(object),..., k = k))
  } else {
    out <- do.call("rbind",
                   lapply(unname(dots), function(x){
                     loglik <- logLik(eval(x))
                     df <- attributes(loglik)$df
                     aic <- AIC(loglik)
                     return(data.frame("model"= deparse(x), "AIC" = aic, "df"= df))
                   }))

    out$dAIC <- out$AIC - min(out$AIC)
    out <- out[order(out$AIC),]
    return(out)
  }
}

#' @title Bayesian information criterion
#' @rdname AIC
#' @export
BIC.generic_null_fit <- function(object,...){
  dots <- as.list(match.call())[-1]

  if(length(dots) == 1){
    return(BIC(logLik(object),...))
  } else {
    out <- do.call("rbind",
                   lapply(unname(dots), function(x){
                     loglik <- logLik(eval(x))
                     df <- attributes(loglik)$df
                     bic <- BIC(loglik)
                     return(data.frame("model"= deparse(x), "BIC" = bic, "df"= df))
                   }))

    out$dBIC <- out$BIC - min(out$BIC)
    out <- out[order(out$BIC),]
    return(out)
  }
}

#' @title Extract Model Coefficients
#' @description extract model coefficient from 'generic_null_fit' objects
#' @param object A 'generic_null_fit' object
#' @param ... additional arguments
#' @param backtransform A logic value indicating whether to back transform coefficients from the scale the coefficient was estimated at to the scale the coefficient is parameterized as for the neutral model. Default is \code{TRUE}.
#' @return A named vector
coef.generic_null_fit<-function(object, ..., backtransform = T){
  if(backtransform){
    out<-vapply(seq_along(object$par), function(i){
      vapply(object$par[i],object$param.val.trans[[object$theta.names[i]]],numeric(1))
    }, numeric(1))
  } else {
    out <- object$par
  }
  names(out) <- object$theta.name
  return(out)
}



#' @export
fit_bite_size<-function(object,
                        family = c("allo","allo_a","allo_M",
                                   "tlnorm","beta","kumar", "tpareto"),
                        min.phi = NA,
                        method = c("Nelder-Mead","BFGS"),
                        IC = "AICc"){

  if(inherits(object,"split_herb")){
    data.vec <- object$hole_prop
    if(is.na(min.phi)){
      min.phi <- object$min_prop
    }
  } else {
    data.vec <- object
  }

  if(is.na(min.phi)){
    min.phi <- min(data.vec[data.vec > 0])
    warning("No min.phi supplied; defaulting to minimum value in data. Might be a bad idea.")
  }

  family <- match.arg(family, several.ok = TRUE)
  method <- match.arg(method)

  if(min.phi <= 0 || min.phi >=1){
    stop("min.phi must be between zero and one, not inclusive.")
  }

  data.vec <- .herb_data_check(data.vec, min.phi = min.phi, allow.zero = FALSE)

  n <- length(data.vec)
  if(n < 3){
    warning("Sample size too low. Some models not identifiable.")
  }

  if("allo" %in% family){
    allo.loglik<- structure(sum(dallo(x = data.vec,
                                      min.phi = min.phi,
                                      max.phi = 1,
                                      a = 14/9,
                                      log = T)),
                            class= "logLik",
                            nall = n,
                            nobs= n,
                            df = 0)

    allo_out<-list("param" = c(NULL),
                   "se" = c(NULL),
                   "logLik" = allo.loglik,
                   "convergence" = ifelse(is.finite(allo.loglik),0,1)
    )
  } else {
   allo_out <- NULL
  }

  if("allo_a" %in% family){

    allo_a.fit<-optim2(par = 1.5,
                      fn = function(theta){
                        -sum(dallo(x = data.vec,
                                   min.phi = min.phi,
                                   max.phi = 1,
                                   a = theta,
                                   log = T))
                      }, method = "Brent",
                      lower = 0,
                      upper = 10,
                      hessian = T)

    allo_a.loglik<- structure(-allo_a.fit$value,
                              class= "logLik",
                              nall = n,
                              nobs= n,
                              df = 1)

    allo_a_out<-list("param" = c("a"=allo_a.fit$par),
                     "se" = tryCatch(sqrt(diag(solve(allo_a.fit$hessian))),
                                     error = function(e) {
                                       rep(NA,length(diag(allo_a.fit$hessian)))
                                     }),
                     "logLik" = allo_a.loglik)
    allo_a_out$convergence <- ifelse(
      is.na((allo_a_out$se)) ||
        allo_a.fit$convergence > 0 ||
        !is.finite(allo_a.fit$value),
      1,
      0
    )
  } else {
    allo_a_out <- NULL
  }

  if("allo_M" %in% family){
    # allo_M.fit<-optim2(par = 1,
    #                   fn = function(theta){
    #                     -sum(dallo(x = data.vec,
    #                                min.phi = min.phi,
    #                                max.phi = theta,
    #                                a = 14/9,
    #                                log = T))
    #                   }, method = "Brent",
    #                   lower = max(data.vec)+0.01,
    #                   upper = 10,
    #                   hessian = T)

    allo_M.loglik<- structure(sum(dallo(x = data.vec,
                                        min.phi = min.phi,
                                        max.phi = max(data.vec),
                                        a = 14/9,
                                        log = T)),
                              class= "logLik",
                              nall = n,
                              nobs= n,
                              df = 1)

    allo_M_out<-list("param" = c("max.phi"=max(data.vec)),
                     "se" = NA,
                     "logLik" = allo_M.loglik)
    allo_M_out$convergence <- ifelse(is.finite(allo_M.loglik),2,1)
  } else {
    allo_M_out <- NULL
  }

  if("tlnorm" %in% family){
    tlnorm.fit<-optim2(par = c(1,0),
                      fn = function(theta){
                        -sum(dhtlnorm(x = data.vec,
                                      theta = 0,
                                      meanlog = theta[1],
                                      sdlog = exp(theta[2]),
                                      endpoint = 1,
                                      log = T))
                      }, method = method,
                      hessian = T)

    tlnorm.loglik<- structure(-tlnorm.fit$value,
                              class= "logLik",
                              nall = n,
                              nobs= n,
                              df = 2)

    tlnorm_out<-list("param" = c("meanlog"=tlnorm.fit$par[1],
                                 "sdlog_log"=tlnorm.fit$par[2]),
                     "se" = tryCatch(sqrt(diag(solve(tlnorm.fit$hessian))),
                                     error = function(e) {
                                       rep(NA,length(diag(tlnorm.fit$hessian)))
                                     }) ,
                     "logLik" = tlnorm.loglik)
    tlnorm_out$convergence <- ifelse(
      any(is.na((tlnorm_out$se))) ||
        tlnorm.fit$convergence > 0 ||
        any(!is.finite(tlnorm.fit$value)),
      1,
      0
    )
  } else {
    tlnorm_out <- NULL
  }

  if("beta" %in% family){
    beta.fit<-optim2(par = c(0,0),
                    fn = function(theta){
                      -sum(dbeta(x = data.vec,
                                 shape1 = exp(theta[1]),
                                 shape2 = exp(theta[2]),
                                 log = T))
                    }, method = method,
                    hessian = T)

    beta.loglik<- structure(-beta.fit$value,
                            class= "logLik",
                            nall = n,
                            nobs= n,
                            df = 2)

    beta_out<-list("param" = c("shape1_log"=beta.fit$par[1],
                               "shape2_log"=beta.fit$par[2]),
                   "se" = tryCatch(sqrt(diag(solve(beta.fit$hessian))),
                                   error = function(e) {
                                     rep(NA,length(diag(beta.fit$hessian)))
                                   }) ,
                   "logLik" = beta.loglik)
    beta_out$convergence <- ifelse(
      any(is.na((beta_out$se))) ||
        beta.fit$convergence > 0 ||
        any(!is.finite(beta.fit$value)),
      1,
      0
    )
  } else {
    beta_out <- NULL
  }

  if("kumar" %in% family){
    kumar.fit<-optim2(par = c(0,0),
                    fn = function(theta){
                      -sum(extraDistr::dkumar(x = data.vec,
                                 a = exp(theta[1]),
                                 b = exp(theta[2]),
                                 log = T))
                    }, method = method,
                    hessian = T)

    kumar.loglik<- structure(-kumar.fit$value,
                            class= "logLik",
                            nall = n,
                            nobs= n,
                            df = 2)

    kumar_out<-list("param" = c("a_log"=kumar.fit$par[1],
                               "b_log"=kumar.fit$par[2]),
                   "se" = tryCatch(sqrt(diag(solve(kumar.fit$hessian))),
                                   error = function(e) {
                                     rep(NA,length(diag(kumar.fit$hessian)))
                                   }) ,
                   "logLik" = kumar.loglik)
    kumar_out$convergence <- ifelse(
      any(is.na((kumar_out$se))) ||
        kumar.fit$convergence > 0 ||
        any(!is.finite(kumar.fit$value)),
      1,
      0
    )
  } else {
    kumar_out <- NULL
  }

  if("tpareto" %in% family){
    tpareto.fit<-optim2(par = 1,
                     fn = function(theta){
                       -sum(dtpareto(x = data.vec,
                                               a = exp(theta),
                                               b = min.phi,
                                               log = T))
                     }, method = "Brent",
                     lower = -10,
                     upper = 4,
                     hessian = T)

    tpareto.loglik<- structure(-tpareto.fit$value,
                             class= "logLik",
                             nall = n,
                             nobs= n,
                             df = 1)

    tpareto_out<-list("param" = c("a_log"=tpareto.fit$par),
                    "se" = tryCatch(sqrt(diag(solve(tpareto.fit$hessian))),
                                    error = function(e) {
                                      rep(NA,length(diag(tpareto.fit$hessian)))
                                    }) ,
                    "logLik" = tpareto.loglik)
    tpareto_out$convergence <- ifelse(
      is.na((tpareto_out$se)) ||
        tpareto.fit$convergence > 0 ||
        !is.finite(tpareto.fit$value),
      1,
      0
    )
  } else {
    tpareto_out <- NULL
  }

  out<-list(allo_out, allo_a_out, allo_M_out,
            tlnorm_out, beta_out, kumar_out, tpareto_out)
  out<-out[!simplify2array(lapply(out,is.null))]
  names(out) <- supported.families[supported.families %in% family]

  out.tab<-simplify2array(lapply(out,function(x){
    c(round(do.call(match.fun(IC),
                    list(x$logLik)),
                   2),
      attr(x$logLik,"df"),x$convergence)
  }))
  out.tab<-out.tab[,order(out.tab[1,]),drop=FALSE]
  out.tab<-rbind(out.tab[1,]-out.tab[1,1],out.tab)[c(2,1,3,4),,drop=FALSE]
  rownames(out.tab) <- c(IC,paste0("d",IC),"df","convergence")
  print(t(out.tab))

  out.final<-list("models" = out,
                  "settings" = list("data" = data.vec,
                         "min.phi" = min.phi,
                         "method" = method,
                         "family" = family),
                  "summary" = out.tab)
  if(any(out.tab[4,] == 1)){
    warning("Model did not converge (1).")
  }

  invisible(out.final)
}


