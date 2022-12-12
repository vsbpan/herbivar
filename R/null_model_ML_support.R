#' @title Fit Non-Neutral Generic Null Models Using Maximum Likelihood Estimation
#' @description Estimate unknown parameters in non-neutral generic null herbivory models from a vector of observed herbivory data using Maximum Likelihood Estimation (MLE).
#' @param data.vec A vector of numeric data. \code{NA}s are ignored. If the whole vector does not contain any non-zero values or contains values outside of \eqn{[0,1]}, the function would throw an error. Values less than "min.phi" are coerced to zero by default.
#' @param family The generic null distribution. Valid options are "htlnorm" for hurdle truncated log-normal distribution, and "zoibeta" for zero-one-inflated beta distribution. see \code{?dzoibeta()} and \code{?dhtlnorm()}
#' @param init Initial value for optimizer to start. Defaults to 0.
#' @param lower,upper A numeric vector indicating the bounds of each estimated variable for the "L-BFGS-B" and "Brent" method. Defaults to \code{-Inf} and \code{Inf}.
#' @param method The method of maximum likelihood estimation algorithm. Acceptable methods are "Nelder-Mead" (Default), "BFGS", "L-BFGS-B", "Brent", "nlminb", and "nlm". Defaults to "Brent" for one dimensional optimization. For details see \code{?stats::optim()} and \code{?stats::nlminb()}.
#' @param id A character string supplied for book keeping purposes. Default is \code{NULL}. Useful for when storing a large number of "generic_null_fit" objects in a list.
#' @param ... Additional arguments that are passed to the optimizer functions \code{optim()} or \code{nlminb()}.
#' @return
#' An 'generic_null_fit' object with the following slots:
#'
#' \code{theta.names}: A vector of character string naming the fitted parameters
#'
#' \code{par}: the value of fitted parameters on the transformed scale in the order of theta.names
#'
#' \code{se}: standard error of fitted parameters on the transformed scale.
#'
#' \code{loglik}: log likelihood of the data given the parameter combinations
#'
#' \code{hessian}: The Hessian matrix
#'
#' \code{iters}: The number of times the likelihood function or its gradient is evaluated.
#'
#' \code{convergence}: 0 indicates successful convergence. Error codes are passed from the optimizer function. 1 indicate that the iteration limit has been reached and convergence has failed. 2 indicate that the Hessian matrix is non-positive definite.
#'
#' \code{optim.param}: A list of parameters passed to \code{optim2()}, including the method of optimization, a vector of the used initial values for optimization, the lower and upper bounds.
#'
#' \code{data}: The data to which the model is fitted to
#'
#' \code{param.val.trans}: A vector of functions that transform each estimated parameter value in the optimization process (the optimizer find the values on the transformed scale).
#'
#' \code{id}: name of model supplied via the id argument
#'
#' \code{n}: the sample size
#'
#' \code{df}: the degrees of freedom of the model
#'
#' \code{herbivar.version}: version of herbivar used to fit the model

#' @export
fit_generic_null<- function(data.vec,
                            family = c("htlnorm","zoibeta"),
                            init = 0,
                            lower = -Inf,
                            upper = Inf,
                            method=c("Nelder-Mead","BFGS","L-BFGS-B","nlminb","nlm"),
                            id = NULL,
                            ...){
  method <- match.arg(method)
  family <- match.arg(family)

  data.vec <- .herb_data_check(data.vec, min.phi = 0, allow.zero = TRUE, coerce.to.zero = TRUE)

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

 ML.fit<-optim2(init = init,
                  fn = nll,
                  method = method,
                  lower = lower,
                  upper = upper,
                  ...)

  generic_null_fit.out<-list("theta.names" = optim.vars,
                          "par" = ML.fit$par,
                          "se" = ML.fit$se,
                          "family" = family,
                          "loglik" = -ML.fit$value,
                          "hessian" = ML.fit$hessian,
                          "iters" = ML.fit$counts,
                          "convergence" = ML.fit$convergence, # > 0 indicates issues
                          "optim.param" = list("init" = init,
                                               "method" = method,
                                               "lower" = lower,
                                               "upper" = upper),
                          "data" = data.vec,
                          "param.val.trans" = param.val.trans,
                          "id" = id,
                          "n" = length(data.vec),
                          "df" = length(optim.vars),
                          "herbivar.version" = herbivar.version(silent = TRUE)

  )
  if(generic_null_fit.out$convergence == 2){
    warning("Non-positive definite hessian matrix")
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
#' @param trans a function that is used to transform the estimated model coefficients. No transformation is applied when set to \code{FALSE}. If set to "backtransform" (default), then the values are back transformed to the original scale.
#' @param digits number of digits to display. Default is 3.
#' @return a matrix
#' @export
print.generic_null_fit<-function(x, ..., trans = "backtransform", digits = 3){
  if(!inherits(x,"generic_null_fit")){
    stop("Needs to by of object type 'generic_null_fit'")
  }

  cat("Fitted result on", x$id,"\n")
  cat("method =",x$optim.param$method,"\n")
  cat("family =",x$family,"\n")
  cat("n =",x$n,"\n")
  cat("\n")
  print(c("loglik" = x$loglik, "AIC" = AIC(x), "AICc" = AICc(logLik(x)), "BIC" = BIC(x)))
  cat("\n")
  print(c(x$iters,"convergence"=x$convergence))
  cat("\n")
  if(trans != "backtransform"){
    cat("Estimate needs backtransformation")
    print(noquote(vapply(
      x$param.val.trans[x$theta.names],
      function(z){
        gsub(".*\\{| |\\}","",paste0(deparse(z),collapse = ""))
      },
      character(1)
    )))
    cat("\n")
  }

  if(!is.null(trans) && trans == "backtransform"){
    back.trans.coef <-coef(object = x, backtransform = TRUE, se = TRUE)
    x$par <- back.trans.coef$Estimate
    x$se <- back.trans.coef$Std.

    FUN <- function(x) x
  } else {
    if(is.null(trans) || isFALSE(trans)){
      FUN <- function(x) x
    } else {
      FUN <- match.fun(trans)
    }

  }

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
#' @return A named vector of estimates. If \code{se = TRUE}, a named list of vectors will be returned instead.
coef.generic_null_fit<-function(object, ..., backtransform = TRUE, se = FALSE){
  if(backtransform){
    out<-vapply(seq_along(object$par), function(i){
      object$param.val.trans[[object$theta.names[i]]](object$par[i])
    }, numeric(1))

    if(se){
      se.out <- vapply(seq_along(object$par), function(i){
        sqrt(FOSM(object$par[i],
                  object$se[i]^2,
                  object$param.val.trans[[object$theta.names[i]]]))
      }, numeric(1))
      names(se.out) <- object$theta.name
    }
  } else {
    out <- object$par
    if(se){
      se.out <- object$se
      names(se.out) <- object$theta.name
    }
  }
  names(out) <- object$theta.name

  if(se){
    return(list("Estimate" = out,
                "Std." = se.out))
  } else {
    return(out)
  }
}


#' @title Fit A 'Bite Size' Distribution To Data
#' @description Fit multiple distributions to a vector of data using maximum likelihood.
#' @param object A vector of proportion 'bite sizes' between zero and one, not inclusive. If an object of class 'split_herb' is supplied, the \code{hole_prop} vector will be extracted.
#' @param family A vector of distributions to fit to a vector of 'bite sizes'. Valid options are "allo", "allo_a", "allo_M", "tlnorm", "beta", "kumar", "tpareto", and "cb". See details. If "all", all available families will be selected.
#' @param min.phi The minimum proportion 'bite size'. If \code{NA} (default), the minimum non-zero value in the data is used. If \code{NA} and \code{object} supplied is of class 'split_herb', the \code{min_prop} is extracted from the \code{object}.
#' @param method Optimizer used to find the maximum likelihood estimates. Valid options are "Nelder-Mead", "BFGS" (default), and "nlminb".
#' @param IC An information criteria used to rank the fitted distributions. Default is "AICc".
#'
#' @details
#'
#' ## allo
#' Neutral 'bite size' distribution with zero degree of freedom.
#' \deqn{P(\phi| \phi_M = 1, \phi_m = min.phi, \alpha = \frac{14}{9}) = \frac{1-\alpha}{\phi_M^{1-\alpha} - \phi_m^{1-\alpha}} \phi^{-\alpha}}.
#'
#'
#' ## allo_a
#' Neutral 'bite size' distribution with one degree of freedom, allowing \eqn{\alpha} to be fitted.
#' \deqn{P(\phi| \phi_M = 1, \phi_m = min.phi, \alpha) = \frac{1-\alpha}{\phi_M^{1-\alpha} - \phi_m^{1-\alpha}} \phi^{-\alpha}}.
#'
#'
#' ## allo_M
#' Neutral 'bite size' distribution with one degree of freedom, allowing \eqn{\phi_M} to be fitted.
#' \deqn{P(\phi| \phi_M, \phi_m = min.phi, \alpha = \frac{14}{9}) = \frac{1-\alpha}{\phi_M^{1-\alpha} - \phi_m^{1-\alpha}} \phi^{-\alpha}}.
#'
#'
#' ## tlnorm
#' Truncated log-normal distribution with two degrees of freedom.
#' \deqn{P(\phi) = \frac{g(\phi, \mu, \sigma)}{G(\phi = 1; \mu ,\sigma)} }
#' The log-normal PDF is defined as
#' \deqn{ g(\phi; \mu , \sigma) = \frac{1}{\phi \sigma \sqrt{2 \pi}} \exp{(-\frac{(\ln{\phi}-\mu)^2}{2\sigma^2})} }
#' and the log-normal CDF is \eqn{G(\phi = 1; \mu ,\sigma)}, \eqn{\mu} is the mean, \eqn{\sigma} is the standard deviation.
#'
#'
#' ## beta
#' Beta distribution with two degrees of freedom.
#' \deqn{P(\phi) = \frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)}\phi^{a-1}(1-\phi)^{b-1}}
#' \eqn{a} and \eqn{b} are the two shape parameters.
#'
#'
#' ## kumar
#' Kumaraswamy distribution with two degrees of freedom.
#' \deqn{P(\phi) = ab\phi^{a-1}(1-\phi^{a-1})(1-\phi^{a})^{b-1}}
#' \eqn{a} and \eqn{b} are the two shape parameters.
#'
#'
#' ## tpareto
#' Truncated Pareto distribution with two degrees of freedom.
#' \deqn{P(x) = \frac{a b^a}{\phi^{a+1}}\frac{1}{1 - (\frac{b}{1})^a},}.
#' \eqn{a} is the exponent of the Pareto distribution and \eqn{b} is the minimum value of \eqn{\phi}.
#'
#'
#' ## cb
#' Continuous Bernoulli distribution with one degree of freedom.
#' \deqn{P(\phi) = C(\lambda) \lambda^\phi (1 - \lambda)^{1-\phi},}
#' \eqn{C(\lambda)} is defined as
#' \deqn{C(\lambda) = 2}
#' if \eqn{\lambda = \frac{1}{2}}, otherwise,
#' \deqn{C(\lambda) = \frac{2 \tanh^{-1}(1-2\lambda)}{1-2\lambda}.}
#'
#' \eqn{\lambda} is the rate parameter.
#'
#' @return
#' A object of class 'bite_size_fit' and 'list'
#'
#'
#' \code{models}:
#' A named list of model fit details, including MLE values, SE, log likelihood, and convergence code.
#'
#' \code{settings}:
#' A named list of function inputs, including the original vector of data, \code{min.phi}, \code{method}, \code{family}, and \code{herbivar.version}.
#'
#' \code{summary}:
#' A matrix array of ordered IC values of fitted models.
#'
#' @export
fit_bite_size<-function(object,
                        family = c("allo","allo_a","allo_M",
                                   "tlnorm","beta","kumar", "tpareto", "cb"),
                        min.phi = NA,
                        method = c("BFGS","nlminb", "Nelder-Mead"),
                        IC = "AICc"){

  supported.families <- c("allo","allo_a","allo_M",
                          "tlnorm","beta","kumar", "tpareto", "cb")
  if(any("all" %in% family)){
    family <- supported.families
  }

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

  data.vec <- .herb_data_check(data.vec, min.phi = min.phi, allow.zero = FALSE, allow.one = FALSE)

  n <- length(data.vec)
  if(n < 3){
    warning("Sample size too low. Some models not identifiable.")
  }

  if("allo" %in% family){
    allo.loglik<- structure(sum(dallo(x = data.vec,
                                      min.phi = min.phi,
                                      max.phi = 1,
                                      a = 14/9,
                                      log = TRUE)),
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

    allo_a.fit<-optim2(init = 0.01,
                      fn = function(theta){
                        -sum(dallo(x = data.vec,
                                   min.phi = min.phi,
                                   max.phi = 1,
                                   a = theta,
                                   log = TRUE))
                      }, method = "Brent",
                      lower = 0,
                      upper = 1)

    allo_a.loglik<- structure(-allo_a.fit$value,
                              class= "logLik",
                              nall = n,
                              nobs= n,
                              df = 1)

    allo_a_out<-list("param" = c("a"=allo_a.fit$par),
                     "se" = allo_a.fit$se,
                     "logLik" = allo_a.loglik,
                     "convergence" = allo_a.fit$convergence)
  } else {
    allo_a_out <- NULL
  }

  if("allo_M" %in% family){
    # allo_M.fit<-optim2(init = 0.9,
    #                   fn = function(theta){
    #                     -sum(dallo(x = data.vec,
    #                                min.phi = min.phi,
    #                                max.phi = exp(theta),
    #                                a = 14/9,
    #                                log = TRUE))
    #                   }, method = "Brent",
    #                   lower = exp(min.phi+0.0001),
    #                   upper = exp(30))
    allo_M.fit <- list("par" = max(data.vec),
                       "se" = NA,
                       "value" = -sum(dallo(x = data.vec,
                                            min.phi = min.phi,
                                            max.phi = max(data.vec),
                                            a = 14/9,
                                            log = TRUE)),
                       "convergence" = 2)

    allo_M.loglik<- structure(-allo_M.fit$value,
                              class= "logLik",
                              nall = n,
                              nobs= n,
                              df = 1)

    # Calculation of the variance of the MLE of the truncation point is nonregular -- cannot be found from the fisher infromation matrix

    allo_M_out<-list("param" = c("max.phi_log" = allo_M.fit$par),
                     "se" = allo_M.fit$se,
                     "logLik" = allo_M.loglik,
                     "convergence" = allo_M.fit$convergence)

  } else {
    allo_M_out <- NULL
  }

  if("tlnorm" %in% family){
    tlnorm.fit<-optim2(init = c(1,0),
                      fn = function(theta){
                        -sum(dhtlnorm(x = data.vec,
                                      theta = 0,
                                      meanlog = theta[1],
                                      sdlog = exp(theta[2]),
                                      endpoint = 1,
                                      log = TRUE))
                      }, method = method,
                      upper = c(Inf, Inf),
                      lower = c(-Inf, -Inf))

    tlnorm.loglik<- structure(-tlnorm.fit$value,
                              class= "logLik",
                              nall = n,
                              nobs= n,
                              df = 2)

    tlnorm_out<-list("param" = c("meanlog"=tlnorm.fit$par[1],
                                 "sdlog_log"=tlnorm.fit$par[2]),
                     "se" = tlnorm.fit$se,
                     "logLik" = tlnorm.loglik,
                     "convergence" = tlnorm.fit$convergence)
  } else {
    tlnorm_out <- NULL
  }

  if("beta" %in% family){
    beta.fit<-optim2(init = c(0,0),
                    fn = function(theta){
                      -sum(dbeta(x = data.vec,
                                 shape1 = exp(theta[1]),
                                 shape2 = exp(theta[2]),
                                 log = TRUE))
                    }, method = method,
                    upper = c(Inf, Inf),
                    lower = c(-Inf, -Inf))

    beta.loglik<- structure(-beta.fit$value,
                            class= "logLik",
                            nall = n,
                            nobs= n,
                            df = 2)

    beta_out<-list("param" = c("shape1_log"=beta.fit$par[1],
                               "shape2_log"=beta.fit$par[2]),
                   "se" = beta.fit$se,
                   "logLik" = beta.loglik,
                   "convergence" = beta.fit$convergence)
  } else {
    beta_out <- NULL
  }

  if("kumar" %in% family){
    kumar.fit<-optim2(init = c(0,0),
                    fn = function(theta){
                      -sum(extraDistr::dkumar(x = data.vec,
                                 a = exp(theta[1]),
                                 b = exp(theta[2]),
                                 log = TRUE))
                    }, method = method,
                    upper = c(Inf, Inf),
                    lower = c(-Inf, -Inf))

    kumar.loglik<- structure(-kumar.fit$value,
                            class= "logLik",
                            nall = n,
                            nobs= n,
                            df = 2)

    kumar_out<-list("param" = c("a_log"=kumar.fit$par[1],
                               "b_log"=kumar.fit$par[2]),
                   "se" = kumar.fit$se,
                   "logLik" = kumar.loglik,
                   "convergence" = kumar.fit$convergence)
  } else {
    kumar_out <- NULL
  }

  if("tpareto" %in% family){
    tpareto.fit<-optim2(init = 1,
                     fn = function(theta){
                       -sum(dtpareto(x = data.vec,
                                               a = exp(theta),
                                               b = min.phi,
                                               log = TRUE))
                     }, method = method,
                     lower = -10,
                     upper = 4)

    tpareto.loglik<- structure(-tpareto.fit$value,
                             class= "logLik",
                             nall = n,
                             nobs= n,
                             df = 1)

    tpareto_out<-list("param" = c("a_log"=tpareto.fit$par),
                    "se" = tpareto.fit$se,
                    "logLik" = tpareto.loglik,
                    "convergence" = tpareto.fit$convergence)
  } else {
    tpareto_out <- NULL
  }

  if("cb" %in% family){
    cb.fit<-optim2(init = 0,
                      fn = function(theta){
                        -sum(dcb(x = data.vec,
                                lambda = plogis(theta[1]),
                                log = TRUE))
                      }, method = method,
                   upper = c(Inf),
                   lower = c(-Inf))

    cb.loglik<- structure(-cb.fit$value,
                             class= "logLik",
                             nall = n,
                             nobs= n,
                             df = 2)

    cb_out<-list("param" = c("lambda_logit"=cb.fit$par[1]),
                    "se" = cb.fit$se,
                    "logLik" = cb.loglik,
                 "convergence" = cb.fit$convergence)
  } else {
    cb_out <- NULL
  }

  out<-list(allo_out, allo_a_out, allo_M_out,
            tlnorm_out, beta_out, kumar_out,
            tpareto_out, cb_out)
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

  out.final<-list("models" = out,
                  "settings" = list("data" = data.vec,
                         "min.phi" = min.phi,
                         "method" = method,
                         "family" = family,
                         "herbivar.version" = herbivar.version(silent = TRUE)),
                  "summary" = t(out.tab))
  if(any(out.tab[4,] == 1)){
    warning("Model did not converge (1).")
  }

  out.final <- structure(out.final,
            class = c("bite_size_fit","list"))

  return(out.final)
}

#' @title Print Values
#' @description Prints fitted object and returns an invisible matrix array of the fitted model coefficients
#' @param x An object of class 'bite_size_fit'
#' @param ... additional arguments
#' @return A matrix array
#' @export
print.bite_size_fit <- function(x,...){
  if(inherits(x,"bite_size_fit")){
    print(x$summary)
    invisible(x$summary)
  } else {
    stop("Object must be of class 'bite_size_fit'.")
  }
}




