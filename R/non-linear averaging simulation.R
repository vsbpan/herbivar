
derivative <- function(x, y) diff(y) / diff(x)


mid_pts <- function(x) x[-1] - diff(x)/2

mid_pts2 <- function(x) mid_pts(mid_pts(x))

derivative2 <- function(x,y) derivative(mid_pts(x), derivative(x, y))

count_inflection_pts<-function(v){
  num<-sum(diff(sign(v[abs(v) > 0.5])) != 0 )
  return(num)
}

center_mass<-function(x,y){
  y<-y-min(y)
  sum(x*y/sum(y))
}

#' @title Non-linear Averaging Simulation
#' @description Estimate Jensen's effect via direct simulation or Jensen's effect potential via second order Taylor approximation.
#' @param model A function or a fitted model that can be handled by the package 'insight'.
#' @param x A vector of numeric value used for the simulation.
#' @param RV.name An optional character value naming the predictor if \code{model} is a fitted model, telling the model which predictor to use for generating predictions.
#' @param nboot An integer value indicating the number of bootstraps for estimating the uncertainty in the Jensen's effect due to the uncertainty of the function. Only relevant for when \code{model} is fitted. Provided \code{model} is treated as exact when it is a function.
#' @param summarise If \code{TRUE}, the output will be summarized. Otherwise, the raw simulated draws will be returned.
#' @param plot If \code{TRUE}, a plot will be generated.
#' @param spaghetti If \code{TRUE} (default), a spaghetti plot will be generated.
#' @details
#'
#' The Jensen's effect \eqn{J} at \eqn{\mathbb{E}[X])} is calculated exactly as:
#'
#' \deqn{J = \mathbb{E}[f(X)] - f(\mathbb{E}[X]).}
#'
#' Jensen's effect is related to the local acceleration as we can see from the second order Taylor approximation:
#'
#' \deqn{J(x) = \frac{1}{2}\frac{d^2}{dx^2} f(x) \mathbb{Var}[X] + \mathcal{O}(x^3).}
#'
#' Therefore, Jensen's effect potential can be expressed in terms of the local second derivative of \eqn{f(x)}, where multiplication by \eqn{\frac{1}{2}\mathbb{Var}[X]} gives us an approximation of \eqn{J(x)}.
#'
#' @return
#'
#' \code{JE} generates a vector of simulated Jensen's effect.
#' \code{JE2} generates a data.frame of simulated Jensen's effect potential (second derivative or acceleration) of each unique value of \code{x}.
#'
#' @rdname JE
#' @examples
#' f <- function(x){
#'  x^3 - 0.5*x^2 + x
#' }
#'
#' x <- rnorm(100, 0, 1)
#' y <- f(x)
#' plot(y~x)
#'
#' m <- lm(y ~ x + I(x^3))
#' y.mpred <- predict(m)
#' plot(y.mpred ~ x)
#'
#' JE(f, x, plot = TRUE)
#' JE(m, x, "x", nboot = 10, plot = TRUE)
#'
#' JE2(f,x, "x", nboot = 10, plot = TRUE)
#' JE2(m,x, "x", nboot = 10, plot = TRUE)
#'
#' @export
JE <- function(model, x, RV.name = NULL, nboot = 1L, summarise = TRUE, plot = FALSE){
  x <- sort(x)
  x.bar <- mean(x)

  if(is.function(model)){
    y <- model(x)
    y.xbar <- model(x.bar)
    y.bar <- mean(y)
  } else {
    .is_inst("insight",stop.if.false = TRUE)
    if(is.null(getS3method("predict", class(model)[1],optional = TRUE))){
      warning(paste0("S3 method for predict() not found for class '", class(m)[1]),"'. Double check model prediction.")
    }
    new.data <- data.frame(c(x.bar,x))
    if(is.null(RV.name)){
      stop("Missing RV.name.")
    }
    if(!RV.name %in% names(insight::get_predictors(model)) || length(RV.name) != 1){
      stop("Must provide only one variable in the fitted model.")
    }
    names(new.data) <- RV.name

    if(!is.numeric(nboot) || nboot < 1){
      stop("'nboot' must be a positive integer.")
    }
    if(nboot == 1){
      y <- as.numeric(insight::get_predicted(model, data = new.data, verbose = FALSE))
      y.xbar <- y[1]
      y <- y[-1]
      y.bar <- mean(y)
    } else {
      y <- as.data.frame(
        insight::get_predicted(model,
                               data = new.data,
                               verbose = FALSE,
                               iterations = nboot))
      y <- y[,grepl("iter",names(y))]
      y.xbar <- as.numeric(y[1,])
      y <- y[-1,]
      y.bar <- colMeans(y)
    }
  }

  J <- y.bar - y.xbar

  if(plot){
    if(is.data.frame(y)){
      y <- rowMeans(y)
      y.bar <- mean(y.bar)
      y.xbar <- mean(y.xbar)
    }
    plot(y~x,type="l",xlab = "x", ylab = "y")
    graphics::points(y.bar~x.bar,cex=1.4)
    graphics::points(y.xbar~x.bar,pch=19,cex=1.4)
    graphics::segments(x0 = x.bar, x1 = x.bar, y0 = y.bar, y1 = y.xbar, col = "blue",lwd = 3)
    graphics::legend("bottomleft",legend = c("E[f(X)]", "f(E[X])"),
           col = c("black","black"), pch = c(1,19))
  }

  if(summarise){
    return(c("mean"= mean(J),"se" = sd(J),quantile(J,probs = c(0.025,0.975))))
  } else {
    return(unname(J))
  }
}

#' @rdname JE
JE2 <- function(model, x, RV.name = NULL, nboot = 1L, summarise = TRUE,
                plot = FALSE, spaghetti = TRUE){
  x <- sort(unique(x))

  if(is.function(model)){
    y <- model(x)
  } else {
    .is_inst("insight",stop.if.false = TRUE)
    if(is.null(getS3method("predict", class(model)[1],optional = TRUE))){
      warning(paste0("S3 method for predict() not found for class '", class(m)[1]),"'. Double check model prediction.")
    }

    new.data <- data.frame(x)
    if(is.null(RV.name)){
      stop("Missing RV.name.")
    }
    if(!RV.name %in% names(insight::get_predictors(model)) || length(RV.name) != 1){
      stop("Must provide only one variable in the fitted model.")
    }
    names(new.data) <- RV.name

    if(!is.numeric(nboot) || nboot < 1){
      stop("'nboot' must be a positive integer.")
    }
    if(nboot == 1){
      y <- as.numeric(insight::get_predicted(model, data = new.data, verbose = FALSE))
    } else {
      y <- as.data.frame(
        insight::get_predicted(model,
                               data = new.data,
                               verbose = FALSE,
                               iterations = nboot))
      y <- y[,grepl("iter",names(y))]
    }
  }

  if(is.null(dim(y))){
    y <- as.data.frame(y)
  }

  mat <- vapply(y,
                FUN = function(y,x){
                  derivative2(x,y)
                },
                x = x,
                FUN.VALUE = numeric(length(x)-2))
  x.md2 <- mid_pts(mid_pts(x))

  if(plot){
    alpha.val <- 1 / (ncol(mat)^0.3)
    for (i in seq_len(ncol(mat))){
      if(i == 1){
        plot(mat[,i]~x.md2, type = "l",
             xlab = "x",
             ylab = "f''(x)",
             ylim = range(mat),
             col = scales::alpha(colour = "grey", 0))
        if(!spaghetti){
          break
        }
      }
      graphics::lines(mat[,i]~x.md2, col = scales::alpha(colour = "blue", alpha.val))
    }
    graphics::abline(h = 0, col = "black", lwd = 2)
    graphics::lines(apply(mat,1,function(x){quantile(x,0.025)}) ~ x.md2, col = "red",
                    lwd = 2, l = "dashed")
    graphics::lines(apply(mat,1,function(x){quantile(x,0.975)}) ~ x.md2, col = "red",
                    lwd = 2, lty = "dashed")
    graphics::lines(rowMeans(mat)~x.md2, col = "red", lwd = 2.4)
  }

  if(summarise){
    a <- t(apply(mat, 1, function(x){
      c("mean" = mean(x), "se" = sd(x), quantile(x, prob = c(0.025, 0.975)))
    }))
  } else {
    a <- unname(mat)
  }
  return(data.frame("x" = x.md2, "acceleration" = a))
}





