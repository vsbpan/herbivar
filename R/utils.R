# #' @export
# combine_data_lists<-function(data_list,data_list2){
#   for (i in seq_along(data_list2)){
#     data_list2[[i]]<-data_list2[[i]][names(data_list[[1]])]
#   }
#   out<-c(data_list,data_list2)
#   return(out)
# }


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
#' @description Compute the variance of a transformed random variable using first order Taylor's approximation (Delta method).
#' @param x The expectation of the random variable
#' @param var The variance of the random variable
#' @param trans A function by which the random variable is transformed
#' @details
#' The first order second moment method can be used to approximate the variance of a transformed random variable provided that we know the first derivative of the transformation function, the mean, and variance of the random variable:
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


#' @title Turn data frame to named vector
#' @description Turn data frame to named vector
#' @param x the data.frame
#' @param id_col the index or column name (as a character string) of the names
#' @param val_col the index or column name (as a character string) of the values
#' @return a named vector
as_named_vector <- function(x, id_col = 1, val_col = 2){
  z <- x[,val_col,drop = TRUE]
  names(z) <- x[,id_col,drop = TRUE]
  return(z)
}

#' @title Evenly spaced sequences between the range of a vector
#' @description Generate evenly spaced sequences between the range of a vector \code{x} using \code{seq()}
#' @param x the vector
#' @param length.out the length of output vector
#' @return a numeric vector
seq_interval <- function(x, length.out = 300){
  seq(min(x),max(x), length.out = length.out)
}


#' @title Ordered statistics
#' @description Find the nth ordered statistics
#' @param x the vector
#' @param n the nth order
#' @return an atomic numeric vector
ordered_stat <- function(x, n){
  x <- Rfast::Sort(x, na.last = NA)

  if(n %% 1 > 0){
    stop("n must be a non-zero integer")
  }

  if(n > 0){
    x[n]
  } else {
    x[length(x) + n]
  }
}

#' @title Reshape an array to long format
#' @description Reshape an array to long format where each row correspond to the index of a value in the nth dimension. The value is stored in the 'val' column.
#' @param x the array
#' @param drop if \code{TRUE}, default is \code{FALSE}, dimensions with length one is removed.
#' @return a data.frame with \code{prod(dim(x))} number of rows and \code{length(dim(x)) + 1} number of columns.
melt <- function (x, drop = FALSE){
  x.dim <- dim(x)
  l <- lapply(
    seq.int(length(x.dim)),
    function(i, name, dim_n){
      temp <- name[[i]]
      if(is.null(temp)){
        temp <- seq.int(dim_n[[i]])
      }
      return(temp)
    },
    name = dimnames(x),
    dim_n = x.dim
  )
  names(l) <- paste0("dim", seq_along(l))



  d <- do.call("expand.grid", l)
  if (drop) {
    d <- d[, x.dim > 1]
  }
  d$val <- c(x)
  return(d)
}


#' @title Reshape an array in long format to an array
#' @description Reshape an array in the long format where each row correspond to the index of a value in the nth dimension and the value is stored in the 'val' column into an nth dimensional array.
#' @param x a matrix or data.frame with named columns. The value of each cell is stored as a 'val' column.
#' @return an array
unmelt <- function(x){
  nm <- colnames(x)

  if(!("val" %in% nm)){
    stop("Expects 'val' as a column in 'x'.")
  }

  nm <- nm[nm != "val"]

  array(
    x[do.call("order", as.list(x[rev(nm)])), "val"],
    dim = Rfast::colMaxs(as.matrix(x[,nm]), value = TRUE)
  )
}


#' @title Length of unique elements
#' @description Find the length of unique elements
#' @param x a vector
#' @return a numeric value
unique_len <- function(x){
  length(unique(x))
}


#' @title Unscale variable
#' @description Generate a function from a vector that can back transform a Z-score transformed vector to its original scale.
#' @param x a vector
#' @return a function
unscale <- function(x){
  function(z){
    mean(x) + sd(x) * z
  }
}

#' @title Add quotes around a character string
#' @description Add quotes around a character string
#' @param x a vector
#' @return an atomic character vector
add_quote <- function(x){
  paste0("'",x,"'")
}

#' @title Plot colors
#' @description Display the color or hex code via \code{ggplot2::geom_tile()}
#' @param x a vector of hex code or acceptable ggplot color names
#' @param label if \code{NULL} (default), the supplied character string via \code{x} is overlayed on the colored tiles. If \code{NA} or \code{FALSE}, nothing is displayed
#' @param label_size The size of the label
#' @param boarder_color The color of the boarder of the tiles
#' @param background The color of the background
#' @param linear If \code{TRUE}, the tiles are displayed in one row. Default is \code{FALSe}.
#' @return A ggplot
plot_colors <- function(x, label = NULL, label_size = 3,
                        boarder_color = "black", background = "grey",
                        linear = FALSE){

  if(linear){
    n <- length(x)
    d <- expand.grid(
      "x" = seq_len(n),
      "y" = 1
    )
  } else {
    n <- ceiling(sqrt(length(x)))
    d <- expand.grid(
      "x" = seq_len(n),
      "y" = seq_len(n)
    )
  }

  if(is.null(label)){
    label <- x
  } else {
    if(all(isFALSE(label)) || all(is.na(label))){
      label <- ""
    }
  }

  if(is.null(boarder_color) || isFALSE(boarder_color)){
    boarder_color <- NA
  }

  n_no_col <- (nrow(d) - length(x))
  if(length(label) == 1){
    label <- rep(label, length(x))
  }

  d <- cbind(d, "hex" = c(x, rep(background,n_no_col)), "lab" = c(label, rep("",n_no_col)))

  g <- ggplot2::ggplot(d,
                  ggplot2::aes(x = x, y = -y)) +
    ggplot2::geom_tile(fill = d$hex, color = boarder_color) +
    ggplot2::geom_text(label = d$lab, size = label_size) +
    ggplot2::theme_void()

  return(g)
}


#' @title Bin values
#' @description Find the closest lower bin of a vector
#' @param x a vector of numeric values
#' @param by the spacing between the bins
#' @param value If \code{TRUE}, return the bin value, otherwise return the bin id
#' @param range A numeric vector specifying the lower and upper limit of the bins
#' @return a numeric vector
bin <- function(x, by = 1, value = TRUE, range = c(min(x), max(x))){
  z <- seq(range[1], range[2], by = by)
  out <- findInterval(x, z)
  if(value){
    return(z[out])
  } else {
    return(out)
  }
}



#' @title Generate uncertainty intervals for binary data using conjugate prior
#' @description Generate uncertainty intervals for binary data using conjugate prior (beta distribution)
#' @param alpha The number of successes
#' @param beta The number of failure
#' @param interval A numeric value indicating the intervals for upper and lower bounds (default to 0.95).
#' @param prior A numeric vector of the alpha and beta value of the conjugate prior. Default is the Bayes-Laplace prior (Tuyl et al. 2008) \eqn{\alpha = 1, \beta = 1}. Other popular choices include 'neutral' prior (Kerman 2011) \eqn{\alpha = 1/3, \beta = 1/3}, and the Jeffreys prior \eqn{\alpha = 0.5, \beta = 0.5}.
#' @references
#' Kerman, J. 2011. Neutral noninformative and informative conjugate beta and gamma prior distributions. Electronic Journal of Statistics 5:1450–1470.
#' Tuyl, F., R. Gerlach, and K. Mengersen. 2008. A Comparison of Bayes-Laplace, Jeffreys, and Other Priors: The Case of Zero Events. The American Statistician 62:40–44.
#' @return a numeric vector
beta_conjugate <- function(alpha, beta, interval = 0.95, prior = c(1,1)){
  a <- alpha + prior[1]
  b <- beta + prior[2]

  lower_prob <- (1 - interval) / 2
  upper_prob <- interval + lower_prob

  out <- data.frame("estimate" = a / (a+b),
                    "std" = sqrt(
                      (a * b) / ((a + b)^2 * (a + b + 1))
                    ),
                    "lower" = stats::qbeta(lower_prob, shape1 = a, shape2 = b),
                    "upper" = stats::qbeta(upper_prob, shape1 = a, shape2 = b)
                    )
  return(out)
}


#' @title Find nearest bin
#' @description Find the closest bin of a vector
#' @param x a vector of numeric values
#' @param grid the vector of the bins
#' @param value If \code{TRUE}, return the bin value, otherwise return the bin id
#' @return a numeric vector
nearest_bin <- function(x, grid, value = TRUE){
  grid_o <- order(grid)
  grid <- grid[grid_o]
  cuts <- c(-Inf, grid[-1]-diff(grid)/2, Inf)
  index <- findInterval(x, cuts)

  if(value){
    return(grid[index])
  } else {
    return(grid_o[index])
  }
}

#' @title Check if values overlaps zero
#' @description Check if values overlaps zero and returns a logical value
#' @param ... Numeric vectors of any length
#' @return an atomic logical value
overlaps_zero <- function(...){
  l <- list(...)


  ifelse(var(
    do.call("c", lapply(l, sign))
  ) == 0,
  FALSE,
  TRUE)
}


#' @title Pre-compile custom functions
#' @description Pre-compile all custom functions in the global environment
#' @param ... additional arguments
#' @return NULL
pre_cmp_fun <- function(...){
  invisible(sapply(
    ls(globalenv()),
    function(x){
      z <- get(x)
      if(is.function(z)){
        assign(x, compiler::cmpfun(z), envir = globalenv())
      }
    }
  ))
  message("Done!")
}






