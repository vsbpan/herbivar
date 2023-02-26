#'
slide_show<-function(data_list,expression,sleep.time=1,length=1){
  for (i in seq_len(length)){
    eval(parse(text=expression))
    Sys.sleep(sleep.time)
    cat(i,"/",length,      "\r")
  }
}


#' @title Adjust Proportion Data
#' @description Adjusts proportion data to get rid of 0's and 1's, and optionally transforms the resulting adjusted proportion.
#' @param x a vector of numeric data bounded between 0 and 1 inclusive.
#' @param trans a character value indicating the transformation to be applied to the nudged data. Valid options are "identity" (default-- no transformation), "logit", "empirical_logit", "probit", "log", and "sqarcsin". See details.
#' @param nudge.method a character value indicating the method of getting rid of 0's and 1's. Valid options are "replace" (default), "none", "smithson", "add", "subtract", and "drop". see details.
#' @param nudge.size a numeric value of the size of a small nudge (\eqn{\epsilon}). Valid character values can also be supplied to choose the method for estimating \eqn{\epsilon}. Valid methods are "macmillan" (default), "warton_min", and "warton_max". see details.
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
#' \code{none}: No method applied. Can be used along with the empirical logit transformation.
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
                        nudge.method = c("replace","none","smithson", "add", "subtract", "drop"),
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



#' @title Summarise vector
#' @description Handy function to summaries a vector of numeric values
#' @param x A vector or matrix of numeric values
#' @param interval The density interval used to find the upper and lower intervals. Default is 0.95.
#' @param na.rm If \code{TRUE} (default is \code{FALSE}), missing values will be removed.
#' @return a vector of numeric values summarising the vector \code{x}.
summarise_vec <- function(x, interval = 0.95, na.rm = FALSE){
  x <- c(x)
  if(na.rm){
    x <- x[!is.na(x)]
  }
  lower.prob <- (1-interval)/2
  upper.prob <- 1 - lower.prob
  if(any(is.na(x))){
    out <- c("mean" = NA_real_,
             "median" = NA_real_,
             "sd" = NA_real_,
             "lower" = NA_real_,
             "upper" = NA_real_)
  } else {
    out <- c("mean" = mean(x),
             "median" = median(x),
             "sd" = sd(x),
             "lower" = unname(quantile(x, prob = lower.prob)),
             "upper" = unname(quantile(x, prob = upper.prob)))
  }
  return(out)
}

#' @title Get Biplot or Triplot
#' @description Generate biplot or triplot using ggplot2.
#' @param x an object that is supported by \code{vegan::scores()}
#' @param choices a vector of length 2 defining the axes to plot
#' @param scaling scaling argument passed to \code{vegan::scores()}
#' @param display a vector of characters defining what to plot
#' @param ellipse a character string passed to the \code{type} argument of \code{ggplot2::stat_ellipse()}. If \code{NA}, no ellipse is plotted.
#' @param group a vector in the same length and order as the sites data used to fit the model that is used to color the site points
#' @param sites_alpha the transparency of the sites plotted as points
#' @param species_color the color of the species labels
#' @return a ggplot object
#' @examples
#'
#' # Simulate some data
#' probed.data <- do.call("rbind", lapply(1:3,
#' function(z){
#' bind_vec(lapply(seq_len(100),
#' function(i){
#' x <- rpois(n = 1000, lambda = z)
#' c(probe_dist(x), "z" = z)
#' }))
#' }))
#' std_probes <- scale(probed.data[,-which(names(probed.data) %in% c("z","N","min"))])
#' probed.data$z <- as.factor(probed.data$z)
#'
#' # Fit a simple RDA model
#' m <- vegan::rda(std_probes ~ z, data = probed.data)
#'
#' get_biplot(m, group = probed.data$z) # Triplot
#'
#' # The function returns a ggplot object, so you can modify the plot as you would normally with ggplot()
#' get_biplot(m, group = probed.data$z) + ggplot2::labs(color = "z") # Change color label
#'
#' get_biplot(m, group = probed.data$z, ellipse = "t") + ggplot2::scale_color_brewer(type = "qual")
#'
get_biplot <- function(x, choices = c(1,2), scaling = 2,
                       display = c("sites", "species", "biplot", "centroids"),
                       ellipse = NA,
                       group = NULL,
                       sites_alpha = 1,
                       species_color = "violetred"){
  display <- match.arg(display,several.ok = TRUE)
  if(.is_inst("vegan",stop.if.false = TRUE)){
    s <- vegan::scores(x,choices = choices,scaling = scaling)
  }
  name <- names(s)
  if(is.null(name)){
    name <- "sites"
    s <- list(s)
  }
  s <- lapply(seq_along(name), function(i,s,name){
    x <- as.data.frame(s[[i]])
    x <- cbind(x, "a" = rownames(x))
    names(x)[length(x)] <- name[i]
    x
  }, s = s, name = name)
  names(s) <- name
  dim1 <- names(s$sites)[1]
  dim2 <- names(s$sites)[2]
  s$dummy <- data.frame(1,2)
  names(s$dummy) <- c(dim1, dim2)

  egv <- vegan::eigenvals(x)
  prec_var <- signif((egv / sum(egv) * 100)[choices], digits = 2)

  g <- s$dummy %>%
    ggplot2::ggplot(ggplot2::aes_string(paste0("x = ", dim1),paste0("y = ", dim2))) +
    ggplot2::geom_vline(ggplot2::aes(xintercept=0),linetype="dashed",color="grey",size=1) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=0),linetype="dashed",color="grey",size=1) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = paste0(dim1, " (",prec_var[1],"%)"), y= paste0(dim2, " (",prec_var[2],"%)"))

  if("sites" %in% display){
    if(!is.null(group)){
      s$sites <- cbind(s$sites, "group" = group)
      g <- g + ggplot2::geom_point(data = s$sites, ggplot2::aes(color = group), alpha = sites_alpha)
    } else {
      g <- g + ggplot2::geom_point(data = s$sites,color = "deepskyblue", alpha = sites_alpha)
    }
  }

  if(!is.na(ellipse)){
    if(!is.null(group)){
      s$sites <- cbind(s$sites, "z" = group)
      g <- g + ggplot2::stat_ellipse(data = s$sites, ggplot2::aes(group = z), color = "navy",
                                     linetype = 4, size = 1)
    } else {
      g <- g + ggplot2::stat_ellipse(data = s$sites, color = "navy",
                                     linetype = 4, size = 1)
    }
  }

  if("species" %in% display && !is.null(s$species)){
    g <- g + ggplot2::geom_text(data = s$species,
                                ggplot2::aes(label = species), color = species_color)
    # maybe include ggrepel::geom_text_repel()
  }

  if("biplot" %in% display && !is.null(s$biplot)){
    s$biplot$x <- 0
    s$biplot$y <- 0
    s$biplot$xend <- s$biplot[,dim1]
    s$biplot$yend <-s$biplot[,dim2]
    g <- g + ggplot2::geom_text(data = s$biplot,
                                ggplot2::aes(label = biplot), color = "black") +
      ggplot2::geom_segment(data = s$biplot, ggplot2::aes(x = x,
                                                          y = y,
                                                          xend = xend,
                                                          yend = yend),
                            arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
                            color = "black",
                            size = 1)
  }

  if("centroids" %in% display && !is.null(s$centroids)){
    g <- g + ggplot2::geom_text(data = s$centroids,
                                ggplot2::aes(label = centroids), color = "darkolivegreen")
  }
  return(g)
}


#' @title Hellinger Transformation
#' @description  The Hellinger transformation for a community matrix where the relative species abundance at each site (row) is square root transformed. If no species is observed at a site, the function returns a zero for the site.
#' @param x a community matrix where each column is a species and each row is a site
#' @return a transformed community matrix
hellinger_trans <- function(x){
  rs <- rowSums(x)
  rs[rs == 0] <- 1
  sqrt(x/rs)
}

#' @title Partial R2 of Fitted Model Predictor
#' @description Generic method of finding the partial R2 of a specified predictor variable by finding the change in r2 when the predictor variable is set to zero.
#' @param object a supported model fit object
#' @param var the name of the variable for which to assess the partial R2
#' @param ... additional parameters
#' @return the partial r2
#' @export
r2_partial <- function(object, var, ...){
  UseMethod("r2_partial")
}




#' @title Partial R2 of Fitted 'brmsfit' Model Predictor
#' @description Calculate the Bayesian partial R2 of a specified predictor variable by finding the change in r2 when the predictor variable is set to zero.
#' @param object A 'brmsfit' object
#' @param var the name of the variable for which to assess the partial R2
#' @param summary if \code{TRUE} (default), the posterior draws are summarized
#' @param robust if \code{TRUE} (default is FALSE), the median instead of the mean is returned as the estimate
#' @param probs The lower and upper interval of posterior summary
#' @param ... additional arguments
#' @return a matrix of partial r2 draws or a vector of partial r2 summary.
#' @rdname r2_partial.brmsfit
#' @export
r2_partial.brmsfit<-function(object, var, summary = TRUE,
                             robust = FALSE, probs = c(0.025,0.975),
                             ...){
  .is_inst("brms", stop.if.false = TRUE)
  .is_inst("matrixStats", stop.if.false = TRUE)
  stopifnot(inherits(object, "brmsfit"))
  pred.data <- object$data
  var.names <- names(object$data)
  if(!is.null(stats::formula(object)$response)){
    warning("Multivariate model is untested. Proceed with caution.")
  }

  var <- match.arg(var, var.names[-1])

  pred.data[,match(var,var.names)] <- 0

  y <- brms::get_y(object)
  ypred <- brms::posterior_epred(object, newdata = pred.data,
                                 allow_new_levels = TRUE)
  ypred.og <- brms::posterior_epred(object)
  var_ypred <- matrixStats::rowVars(ypred)
  var_ypred.og <- matrixStats::rowVars(ypred.og)
  var_e <- matrixStats::rowVars(-1 * sweep(ypred, 2, y))
  var_e.og <- matrixStats::rowVars(-1 * sweep(ypred.og, 2, y))
  partial_r2 <- as.matrix(var_ypred.og / (var_ypred.og + var_e.og) - var_ypred /(var_ypred + var_e) )

  if(summary){
    partial_r2 <- brms::posterior_summary(partial_r2, probs = probs, robust = robust)
  }
  rownames(partial_r2) <- "partial_R2"
  return(partial_r2)
}

