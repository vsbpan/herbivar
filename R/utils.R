#' @export
combine_data_lists<-function(data_list,data_list2){
  for (i in seq_along(data_list2)){
    data_list2[[i]]<-data_list2[[i]][names(data_list[[1]])]
  }
  out<-c(data_list,data_list2)
  return(out)
}


.herb_data_check <- function(x, min.phi, allow.zero, allow.one = TRUE, coerce.to.zero = FALSE){
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
      if(coerce.to.zero){
        x[z] <- 0
        message(sum(z)," non-zero values less than min.phi detected and coerced to zero")
      } else {
        x <- x[!z]
        message(sum(z)," non-zero values less than min.phi detected and removed from data")
      }
    }
  } else {
    z <- x < min.phi
    if(any(z)){
        x <- x[!z]
        message(sum(z)," values less than min.phi detected and removed from data")
    }
  }



  if(!allow.one){
    o <- x == 1
    if(any(o)){
      message(sum(o), " ones detected and removed from data")
    }
    x <- x[!o]
  }
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


#' @title Maximum Likelihood Optimization With Error Handling
#' @description A wrapper function of \code{stats::optim()}, \code{stats::nlm()}, and \code{stats::nlminb()} with error handling. Intended as an interval function.
#' @param init initial values for the parameters to be optimized over
#' @param fn function for which to find the minimum
#' @param method method allowed by \code{stats::optim()}, as well as \code{nlm} and \code{nlminb}
#' @param lower,upper lower and upper bounds
#' @param silent if \code{TRUE}, no message is returned. Default to \code{TRUE}.
#' @param ... extra arguments passed to \code{optim()}, \code{nlm()}, or \code{nlminb()}
#' @return a list
#' @export
optim2 <- function(init, fn, method, lower, upper, silent = TRUE, ...){
  mcall<-match.call()
  stopifnot(!missing(method))
  stopifnot(!missing(lower))
  stopifnot(!missing(upper))

  fun <- function(foo,...){
    if(method == "nlm"){
      out <- stats::nlm(p = init,
                 f = fn,
                 hessian = TRUE,
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
          hessian = TRUE,
          ...)
    }

    out$se <- tryCatch(sqrt(diag(solve(out$hessian))),
                       error = function(e) {
                         rep(NA, length(init))
                       })

    if(any(is.na(out$se))){
      out$convergence <- 2
    }
    if(is.na(out$value) || !is.finite(out$value)){
      out$convergence <- 1
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
           "se" = rep(NA, length(init)),
           "counts" = c("function" = NA, "gradient" = NA),
           "convergence" = 1,
           "message" = "Model did not converge",
           "hessian" = NA)
    })
}


.is_inst <- function(pkg, stop.if.false = FALSE, prompt = FALSE) {
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

#' @title Bind a list of vector into data.frame
#' @description Bind a list of vector into a data.frame and keep track of list names and vector names.
#' @param x a list of vectors
#' @param margin the direction for which the vectors will be bound. 1 indicates rows (default) and 2 indicates columns.
#' @param keep_row_names if \code{TRUE} (default), the row names will be kept. Otherwise, they will be set to \code{NULL}.
#' @param row_names_as_col if a character string is provided, the row name will be added to the data.frame as the first column with that name. If \code{TRUE}, "rownames" will be used as the column name. Otherwise, the no column is added (default is \code{FALSE}).
#' @return a data.frame
bind_vec <- function(x,margin = 1L, keep_row_names = TRUE, row_names_as_col = FALSE){
  out <- as.data.frame(do.call("cbind",x))
  names(out) <- names(x)
  if(as.numeric(margin) == 1){
    out <- as.data.frame(t(out))
  }
  if(isTRUE(row_names_as_col)){
    row_names_as_col <- "rownames"
  }

  if(is.character(row_names_as_col)){
    out <- cbind("v" = rownames(out), out)
    names(out)[1] <- row_names_as_col
  }

  if(!keep_row_names){
    rownames(out) <- NULL
  }

  return(out)
}

#' @title First Order Second Moment Method
#' @description Compute the variance of a transformed random variable using first order Taylor's aproximation.
#' @param x The expectation of the random variable
#' @param var The variance of the random variable
#' @param trans A function by which the random variable is transformed
#' @details
#' The first order second moment method can be used to approximate the variance of a transformed random variable provided that we know the first derivitive of the transformation function, the mean, and variance of the random variable:
#' \deqn{\mathbb{Var}[f(X)] \approx f'(\mathbb{E}[X])^2 \mathbb{Var}[X]}
#'
#' @return a numeric value
FOSM <- function(x,var,trans){
  # First order taylor's approximation
  # Var[f(X)] \approx f'(E[X])^2  var[X]
  stopifnot(is.function(trans))
  fun.string <- deparse1(trans)
  variable_name <- gsub(".*function \\(|\\).*","",fun.string)
  fun.string <- gsub(variable_name,"x",fun.string)
  fun.string <-gsub(".*\\{| |\\}.*| ","",fun.string)
  fun.string <- gsub("function\\(x\\)","",fun.string)
  if(fun.string == "plogis(x)"){
    deriv.expression <- expression(1/x + 1/(1-x))
  } else {
    deriv.expression <- stats::D(parse(text = fun.string), "x")
  }

  var.trans <- eval(deriv.expression)^2 * var
  return(var.trans)
}
