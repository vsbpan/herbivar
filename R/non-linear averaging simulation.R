
derivative <- function(x, y) diff(y) / diff(x)


mid_pts <- function(x) x[-1] - diff(x)/2

derivative2 <- function(x,y) derivative(mid_pts(x), derivative(x, y))

count_inflection_pts<-function(v){
  num<-sum(diff(sign(v[abs(v) > 0.5])) != 0 )
  return(num)
}

center_mass<-function(x,y){
  y<-y-min(y)
  sum(x*y/sum(y))
}

#' @title
#' @description
#' @param model
#' @param x
#' @param RV.name
#' @param nboot
#' @param summarise
#' @param plot
#' @param spaghetti
#' @details
#' @return
#' @rdname JE
#' @examples
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
      graphics::lines(mat[,i]~x.md2, col = scales::alpha(colour = "black", alpha.val))
    }
    graphics::abline(h = 0, col = "blue", lwd = 2)
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




