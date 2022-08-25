#### Convenience functions ####

#' @export
slide.show<-function(data.list,expression,sleep.time=1,length=1){
  for (i in seq_len(length)){
    eval(parse(text=expression))
    Sys.sleep(sleep.time)
    cat(i,"/",length,      "\r")
  }
}



#' @export
combine.data.lists<-function(data.list,data.list2){
  for (i in seq_along(data.list2)){
    data.list2[[i]]<-data.list2[[i]][names(data.list[[1]])]
  }
  out<-c(data.list,data.list2)
  return(out)
}


.herb_data_check <- function(x, min.phi, allow.zero){
  x <- as.numeric(x)
  if(any(is.na(x))){
    message(sum(is.na(x))," NA detected and removed from data")
    x <- x[!is.na(x)]
  }

  if(sum(x) == 0){
    stop("Data does not contain non-zero herbivory values; model not identifiable.")
  }

  if(any(x < 0) || any(x > 1)){
    stop("Data must be bounded between 0 and 1.")
  }

  if(allow.zero){
    z <- (x < min.phi) & (x != 0)
    if(any(z)){
      message(sum(z)," non-zero values less than min.phi detected and removed from data")
    }
  } else {
    z <- x < min.phi
    if(any(z)){
      message(sum(z)," values less than min.phi detected and removed from data")
    }
  }
  x <- x[!z]

  return(x)
}


.choose_new_theta_val <- function(x, new.param, theta.name){
  ifelse(
    is.na(new.param[theta.name]),
    ifelse(
      is.na(x$param.vals[theta.name]),
      match.fun(x$param.val.trans[[theta.name]])(x$par[(x$theta.names == theta.name)]),
      x$param.vals[theta.name]
    ),
    new.param[theta.name]
  )
}


#' @title General-Purpose Optimization With Error Handling
#' @description A wrapper function of \code{stats::optim()} with error handling.
#' @param silent if \code{TRUE}, no message is returned. Default to \code{TRUE}.
#' @param ... extra arguments passed to \code{optim()}
#' @return a list
#' @export
optim2 <- function(silent=TRUE,...){
  mcall<-match.call()
  tryCatch(
    if(silent){
      suppressWarnings(optim(...))
    } else{
      optim(...)
    },
    error = function(e){
      if(!silent){
        message("Optimization failed: \n", deparse(mcall))
      }
      list("par" = NA,
           "value" = NA,
           "counts" = c("function" = NA, "gradient" = NA),
           "convergence" = 1,
           "message" = "Model did not converge",
           "hessian" = NA)
    })
}



#' @title Adjust Proportion Data
#' @description Adjusts proportion data to get rid of 0's and 1's, and optionally transforms the resulting adjusted proportion.
#' @param x a vector of numeric data bounded between 0 and 1 inclusive.
#' @param trans a character value indicating the transformation to be applied to the nudged data. Valid options are "identity" (default-- no transformation), "logit", "empirical_logit", "probit", "log", and "sqarcsin". See details.
#' @param nudge.method a character value indicating the method of getting rid of 0's and 1's. Valid options are "replace" (default), "smithson", "add", "subtract", and "drop". see details.
#' @param nudge.size a numeric value of the size of a small nudge (\eqn{epsilon}). Valid character values can also be supplied to choose the method for estimating \eqn{epsilon}. Valid methods are "macmillan" (default), "warton_min", and "warton_max". see details.
#' @param bounds a character value indicating bounds for which the \code{nudge.method} is applied to. Valid options are "both" (default), "zero", and "one". Only relevant for \code{nudge.method} set as "replace" or "drop".
#' @param na.action a character value indicating how to handle missing values. Valid options are "ignore", "remove", "as.is", and "fail". See details.
#'
#' @return a vector of numeric value with the same length as \code{x}
#' @details
#'
#' **Add commentary**
#' ## Handling 1's and 0's.
#' ### Nudge method
#'
#' \code{replace}: Replace 0 by \eqn{\epsilon} and 1 by \eqn{1-\epsilon}
#'
#' \code{smithson}: the vector is transformed according to the recommendation of Smithson & Verkuilen (2006): \eqn{\frac{x(n-1) + \epsilon}{n}}
#'
#' \code{add}: the vector is transformed as \eqn{x + \epsilon}
#'
#' \code{subtract}: the vector is transformed as \eqn{x - \epsilon}
#'
#' \code{drop}: 0's and/or 1's are dropped from the vector
#'
#'
#' ### Epsilon
#'
#' \code{macmillan}: estimated as \eqn{\frac{1}{2n}} per the recommendation of Macmillan & Creelman (2004), where \eqn{n} is the sample size.
#'
#' \code{warton_min}: estimated as the smallest non-zero value in the vector x per Warton & Hui (2011).
#' \code{warton_max}: estimated as one minus the largest non-one value in the vector x, useful when \code{warton_min} is too big (Warton & Hui 2011).
#'
#'
#' ## Transformation
#'
#' \code{identity}: no transformation
#'
#' \code{logit}: \eqn{\ln{\frac{x}{1-x}}}. Undefined when \eqn{x = 1} or \eqn{x = 0}. Recommended by Warton & Hui (2011) over the square root arcsine transformation when transforming true proportions for analysis that assumes a Gaussian distribution.
#'
#' \code{empirical_logit}: \eqn{\ln{\frac{x+\epsilon}{1-x+\epsilon}}}. A modified logit transformation recommend by Warton & Hui (2011) when 0's or 1's are present in the data.
#'
#' \code{probit}: \eqn{\int_{-\infty}^{x}\frac{1}{\sqrt{2 \pi}}\exp{-\frac{t^2}{2}}dt}
#'
#' \code{log}: \eqn{\ln(x)}
#'
#' \code{sqarcsin}: \eqn{\arcsin(\sqrt{x})}
#'
#'
#' ## Missing Value Handling
#'
#' \code{ignore}: Function will ignore missing values in the data and return a vector of adjusted values with the original missing values in place.
#'
#' \code{remove}: Function will remove missing values in the data and return a shorter vector or adjusted values.
#'
#' \code{as.is}: Function will perform the computation with the missing values, possibly returning a vector of \code{NA}s.
#'
#' \code{fail}: Function will throw an error if the data contains any missing value.
#'
#' @references
#' Kubinec, R. 2022. Ordered Beta Regression: A Parsimonious, Well-Fitting Model for Continuous Data with Lower and Upper Bounds. Political Analysis:1–18.
#'
#' Macmillan, N. A., and C. D. Creelman. 2004. Detection Theory: A User’s Guide. Second edition. Psychology Press, New York.
#'
#' Smithson, M., and J. Verkuilen. 2006. A better lemon squeezer? Maximum-likelihood regression with beta-distributed dependent variables. Psychological Methods 11:54–71.
#'
#' Warton, D. I., and F. K. C. Hui. 2011. The arcsine is asinine: the analysis of proportions in ecology. Ecology 92:3–10.
#' @export
#'
#' @examples
#' x <- c(0,runif(0,0,1),NA,1,1,1)
#' adjust_prop(x, nudge.method = "replace", nudge.size = 0.01)
#'
#' adjust_prop(x, nudge.method = "replace", nudge.size = "macmillan", trans = "logit")
#'
#'
adjust_prop <- function(x,
                        nudge.method = c("replace","smithson", "add", "subtract", "drop"),
                        nudge.size = c("macmillan", "warton_min", "warton_max"),
                        trans = c("identity","logit","empirical_logit","probit","log","sqarcsin"),
                        bounds = c("both","zero","one"),
                        na.action = c("ignore","remove","as.is","fail")){

  nudge.method <- match.arg(nudge.method)
  bounds <- match.arg(bounds)
  trans <- match.arg(trans)
  na.action <- match.arg(na.action)

  bounds.numeric <- switch(bounds,
                           "both" = c(0,1),
                           "zero" = 0,
                           "one" = 1)

  if(is.character(nudge.size)){
    nudge.size <- match.arg(nudge.size)
  } else {
    nudge.size <- as.numeric(nudge.size)[1]
  }

  na.v <- is.na(x)

  if(any(x[!na.v] > 1)){
    stop("'x' contains values greater than 1")
  }
  if(any(x[!na.v] < 0)){
    stop("'x' contains negative values")
  }

  if(nudge.method == "drop"){
    x[x %in% bounds.numeric] <- NA
    if(na.action != "remove"){
      if(na.action != "ignore"){
        na.action <- "ignore"
        message("na.action defaulting to 'ignore'")
      }
    }
  }

  if(any(na.v)){
    if(na.action == "fail"){
      stop("'x' contains missing values")
    }

    if(na.action == "remove"){
      x <- x[!na.v]
    } else {
      if(na.action == "ignore"){
        x_na_rm <- x[!na.v]
        warning("'x' contains ignored missing values.")
      } else {
        warning("'x' contains missing values. Values computed as is.")
      }
    }
  } else {
    x_na_rm <- x
  }

  if(is.numeric(nudge.size)){
    epsilon <- nudge.size
    nudge.size <- "user_supplied"
  } else {
    if(nudge.method == "smithson"){
      message("epsilon defaulting to 0.5")
      nudge.size <- "smithson_default"
      epsilon <- 0.5
    } else {
      if(na.action == "ignore"){
        epsilon <- switch(nudge.size,
                          "macmillan" = 1/(2*length(x_na_rm)),
                          "warton_min" = min(x_na_rm[x_na_rm > 0]),
                          "warton_max" = min((1-x_na_rm)[(1-x_na_rm) > 0])
        )
      } else {
        epsilon <- switch(nudge.size,
                          "macmillan" = 1/(2*length(x)),
                          "warton_min" = min(x[x > 0]),
                          "warton_max" = min((1-x)[(1-x) > 0])
        )
      }
    }
  }

  if(nudge.method == "replace"){
    if(0 %in% bounds.numeric){
      x <- ifelse(x == 0,
                  epsilon,
                  x)
    }
    if(1 %in% bounds.numeric){
      x <- ifelse(x == 1,
                  1-epsilon,
                  x)
    }
  } else {
    if(nudge.method == "add"){
      x <- x + epsilon
    } else {
      if(nudge.method == "subtract"){
        x <- x - epsilon
      } else {
        if(nudge.method == "smithson"){
          if(na.action == "ignore"){
            n <- length(x_na_rm)
          } else {
            n <- length(x)
          }
          x <- (x * (n-1) + epsilon) / n
        }
      }
    }
  }

  if(any(x[!is.na(x)] > 1)){
    warning("nudged 'x' contains values greater than 1")
  }
  if(any(x[!is.na(x)] < 0)){
    warning("nudged 'x' contains negative values")
  }

  x.trans <- switch(trans,
                    "identity" = x,
                    "logit" = qlogis(x),
                    "empirical_logit" = log((x + epsilon) / (1 - x + epsilon)),
                    "probit" = qnorm(x),
                    "log" = log(x),
                    "sqarcsin" = asin(sqrt(x)))
  x.trans <- structure(x.trans,
                       "adj.method" = list(
                         "epsilon" = epsilon,
                         "nudge.size" = nudge.size,
                         "nudge.method" = nudge.method,
                         "trans" = trans,
                         "bounds" = bounds,
                         "na.action" = na.action
                       )
  )
  return(x.trans)
}





