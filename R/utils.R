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
#' @description A wrapper function of \code{stats::optim()}, \code{stats::nlm()}, and \code{stats::nlminb()} with error handling. Intended as an interval function.
#' @param init initial values for the parameters to be optimized over
#' @param fn function for which to find the minimum
#' @param method method allowed by \code{stats::optim()}, as well as \code{nlm} and \code{nlminb}
#' @param lower,upper lower and upper bound
#' @param silent if \code{TRUE}, no message is returned. Default to \code{TRUE}.
#' @param ... extra arguments passed to \code{optim()}
#' @return a list
#' @export
optim2 <- function(init, fn, method, lower, upper, silent = TRUE, ...){
  mcall<-match.call()
  fun <- function(foo,...){
    if(method == "nlm"){
      out <- stats::nlm(p = init,
                 f = fn,
                 hessian = T,
                 ...)
      out$objective <- out$minimum
      out$counts <- out$iterations
      out$par <- out$estimate
      out$convergence <- ifelse(out$code > 2 , 1 , 0) # rename the convergence code for consistency
      out$message <- paste0("code = ", out$code)
    }
    if(method == "nlminb"){
      out <- stats::nlminb(start = init,
                    objective = fn,
                    lower = lower,
                    upper = upper,
                    ...)
      out$value <- out$objective
      out$counts <- out$evaluations
      out$hessian <- hessian(nll, out$par)
    } else {
      out <- stats::optim(
          par = init,
          fn = fn,
          method = method,
          upper = upper,
          lower = lower,
          ...)
    }
    return(out)
  }

  tryCatch(
    if(silent){
      suppressWarnings(fun(...))
    } else {
      fun(...)
    },
    error = function(e){
      if(!silent){
        message("Optimization failed: \n", deparse(mcall))
      }
      list("par" = rep(NA, length(init)),
           "value" = NA,
           "counts" = c("function" = NA, "gradient" = NA),
           "convergence" = 1,
           "message" = "Model did not converge",
           "hessian" = NA)
    })
}


.is_inst <- function(pkg, stop.if.false = FALSE, prompt = TRUE) {
  out <- nzchar(system.file(package = pkg))
  if(!out){
    message(pkg," is not installed. Install ",  pkg, " and dependencies?")
    if(prompt && menu(c("Yes","No")) == 1){
      install.packages(pkg)
      out <- nzchar(system.file(package = pkg))
    }
  }
  if(stop.if.false && !out){
    stop("Required package ", pkg, " not found")
  }
  return(out)
}




