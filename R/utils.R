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


.choose_theta_val <- function(param.vals,param.val.trans,optim.vars,theta.name){
  ifelse(
    is.na(param.vals[theta.name]),
    param.val.trans[[theta.name]](eval(parse(
      text =
        paste0("theta[",
               which(optim.vars==theta.name),"]")))),
    param.vals[theta.name])
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







