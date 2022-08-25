#' @title  Calculate Coefficient of Variation (CV)
#' @description calculate the CV for a vector of data
#' @param  x
#'  A vector of numeric values
#' @param  method
#'  A value indicating the type of CV estimator
#' @param  na.rm
#'  A logical value indicating whether to ignore \code{NA} values in the \code{x} vector
#' @return
#'  A numeric value
#' @details
#'  The method option selects one of four CV estimators investigated by Yang et al. (2020).
#'
#'  **CV1**: This is the conventional CV estimator:
#'  \deqn{CV_1 = \frac{\sigma}{\mu}}
#'
#'  **CV2**: This assumes the sample is normally distributed and corrects for the bias using the sample size as in Sokal and Rohlf (1995):
#'  \deqn{CV_2 = CV_1 + \frac{CV_1}{4N}}
#'
#'  **CV3**: This implements Breunig's (2001) method:
#'  \deqn{CV_3 = CV_1 \sqrt{1- \frac{CV_1}{N}(3CV_1 -2\gamma_1)}}
#'  \eqn{\gamma_1} is Pearson's measure of skewness
#'
#'  **CV4**: This implements Bao's (2009) method:
#'  \deqn{CV_4 = CV_1 + \frac{CV_1^3}{N}+ \frac{CV_1}{4N}+ \frac{CV_1^2\gamma_1}{2N} + \frac{CV_1\gamma_2}{8N}}
#'  \eqn{\gamma_2} is Pearson's measure of kurtosis
#'
#' @references
#'
#' Bao, Y. 2009. Finite-Sample Moments of the Coefficient of Variation. Econometric Theory 25:291–297.
#'
#' Sokal, R. R., and J. F. Rohlf. 1995. Biometry: The Principles and Practices of Statistics in Biological Research. Third edition. Freeman, New York.
#'
#' Yang, J., J. Lu, Y. Chen, E. Yan, J. Hu, X. Wang, and G. Shen. 2020. Large Underestimation of Intraspecific Trait Variation and Its Improvements. Frontiers in Plant Science 11.
#'
#' @export
cv<-function(x, method = c(1,2,3,4),na.rm = FALSE){
  if(isTRUE(any(x < 0))) {
    warning("Data contain non-positive values")
  }
  if(na.rm){
    x <- x[!is.na(x)]
  }
  cv1<-sd(x)/mean(x)
  if(method[1] == 1) {
    return(cv1)
  } else
    if(method[1] == 2) {
      cv2<- cv1 + cv1 / (4 * length(x))
      return(cv2)
    } else
      if(method[1] == 3) {
        cv3 <- cv1*sqrt(1-cv1/length(x)*(3*cv1 - 2 * Skew(x)))
        return(cv3)
      } else
        if(method[1] == 4){
          N <- length(x)
          cv4 <- cv1 + cv1^3/N + cv1/(4*N) + cv1^2*Skew(x)/(2*N) + cv1*Kurt(x)/(8*N)
          return(cv4)
        }
}

#' @title  Calculate Aggregation Coefficient
#' @description calculate the aggregation \eqn{J} for a given mean and variance or vector of data
#' @param x a vector of numeric data
#' @param mean mean of data. Ignored if \code{x} is supplied
#' @param var variance of data. Ignored if \code{x} is supplied
#' @param na.rm a logical value indicating whether to drop NA values
#' @return
#'  A numeric value
#' @details
#' The aggregation coefficient \eqn{J} is calculated as
#'  \deqn{J = \frac{\sigma^2 - \mu}{\mu^2}}
#'  \eqn{J = 0} corresponds to a possion distribution (random)
#'  \eqn{J < 0} corresponds to a more even distribution (underdispersed)
#'  \eqn{J > 0} corresponds to a more aggregated distribution (overdispersed)
#'  @export

J.index <- function(x=NULL,mean=NULL,var=NULL, na.rm = FALSE){
  if(any(x < 0)) {
    warning("Data contain non-positive values")
  }
  if(!is.null(x)){
   var <- var(x, na.rm = na.rm)
   mean <- mean(x,na.rm = na.rm)
  }
  (var - mean) / (mean)^2
}


#' @title  Calculate Coefficient of Dispersion (CD)
#' @description calculate CD for a given mean and variance or vector of data
#' @param x a vector of numeric data
#' @param robust if \code{TRUE}, a more robust for of CD is calculated. Default is \code{FALSE}
#' @param na.rm a logical value indicating whether to drop \code{NA} values
#' @return
#'  A numeric value
#' @details
#' The coefficient of dispersion \eqn{CD}, a.k.a. variance mean ratio, is calculated as
#'  \deqn{CD = \frac{\sigma^2}{\mu}}
#' The robust version is calculated as
#' \deqn{CD = \frac{1}{n} \frac{\sum_i^n |m - x_i|}{m}}
#' where \eqn{m} is the median of the data
#' @export
cd <- function(x, robust = FALSE, na.rm = FALSE){
  if(any(x < 0)) {
    warning("Data contain non-positive values")
  }
  if(robust){
    m <- median(x, na.rm = na.rm)
    mean(abs(x - m),na.rm = na.rm) / m
  } else {
    var(x,na.rm = na.rm) / mean(x, na.rm = na.rm)
  }
}


#' @title Calculate Standard Error of The Mean (SE)
#' @description calculate standard error for a vector of data
#' @param x a vector of numeric data
#' @param na.rm a logical value indicating whether to drop \code{NA} values
#' @return
#'  A numeric value
#' @details
#' The standard error of the mean \eqn{SE}, is calculated as
#'  \deqn{SE(X) = \sqrt{\frac{\matbb{Var}[X]}{n}}}
#' @export
se <- function(x, na.rm = FALSE){
  if(na.rm){
    x <- x[!is.na(x)]
  }
  sd(x)/sqrt(length(x))
}

#' @title Get Sample Size
#' @description Basically an alias for \code{length()}
#' @param x a vector of numeric data
#' @param na.rm a logical value indicating whether to drop \code{NA} values
#' @return
#'  A numeric value
#' @export
N <- function(x, na.rm = FALSE){
  if(na.rm){
    x <- x[!is.na(x)]
  }
  return(length(x))
}

#' @title Get Median
#' @description Find median of a vector of numeric data
#' @param x a vector of numeric data
#' @param na.rm a logical value indicating whether to drop \code{NA} values
#' @return
#'  A numeric value
#' @export
q50 <- function(x, na.rm = FALSE){
  median(x, na.rm = na.rm)
}

#' @title Get 25th Percentile
#' @description Find the lower quantile of a vector of numeric data
#' @param x a vector of numeric data
#' @param na.rm a logical value indicating whether to drop \code{NA} values
#' @return
#'  A numeric value
#' @export
q25 <- function(x, na.rm = FALSE){
  if(na.rm){
    x <- x[!is.na(x)]
  }
  return(as.numeric(quantile(x,probs = 0.25)))
}

#' @title Get 75th Percentile
#' @description Find the upper quantile of a vector of numeric data
#' @param x a vector of numeric data
#' @param na.rm a logical value indicating whether to drop \code{NA} values
#' @return
#'  A numeric value
#' @export
q75 <- function(x, na.rm = FALSE){
  if(na.rm){
    x <- x[!is.na(x)]
  }
  return(as.numeric(quantile(x,probs = 0.75)))
}

#' @title Lorenz Asymmetry Coefficient
#' @description Calculate Lorenz Asymmetry Coefficient (also known as \eqn{S}) for a vector of data. The statistic identifies which quantile of \code{x} contribute most to the total inequality.
#' @param x a vector of non-negative numeric values
#' @param n a vector of frequencies the same length as \code{x}. If \code{NULL} (default), the frequencies are assumed to be 1 (as in a vector of untabulated raw data).
#' @param na.rm a logical value indicating whether to drop \code{NA} values
#' @param interval a logical value indicating whether to return an interval when there are x values exactly equal to the mean. If \code{FALSE} (default), the midpoint is returned.
#' @note Part of the code for this function is adapted from Achim Zeileis's \code{ineq::Lasym()} version 0.2-13.
#' @return a single numeric value if \eqn{x_i \neq \hat\mu} for all values of \code{x} or if \code{interval = FALSE}. Otherwise, a vector of length two is returned.
#' @details
#' The Lorenz Asymmetry Coefficient \eqn{S} is a complementary statistic to the Gini coefficient that describes the Lorenz curve. Whereas the Gini coefficient is a measure of inequality, the asymmetry coefficient is is a measure of which half of the \code{x} quantile contribute most to the inequality. If \eqn{S = 1}, the Lorenz curve is symmetrical. If \eqn{S < 1}, the point where the Lorenz curve is parallel to the line of equality is below the axis of symmetry. This means that the lower half of the population contribute most to the inequality. If \eqn{S > 1}, the point where the Lorenz curve is parallelt to the line of equality is above the axis of symmetry. This means that the upper half of the population contribute most to the inequality. See example plot.
#'
#'
#' For a vector of ordered non-negative value \eqn{(x_1, x_2, ..., x_m, x_{m+1}, ..., x_n)}, the sample Lorenz Asymmetry is defined as:
#' \deqn{S = F(\hat\mu) + L(\hat\mu)}
#'
#' where
#' \deqn{\delta = \frac{\hat\mu - x_m}{x_{m+1} - x_m},}
#'
#' \deqn{F(\hat\mu) = \frac{m + \delta}{n},}
#'
#' \deqn{L(\hat\mu) = \frac{L_m + \delta x_{m+1}}{L_m}.}
#'
#' \eqn{\hat\mu} is the mean of \eqn{x}, \eqn{n} is the sample size, \eqn{m} is the number of elements where \eqn{x_m < \hat\mu}, and \eqn{L_i = \sum_{j=1}^i x_j}.
#'
#'
#' If one or more values of \eqn{x_i = \hat\mu}, the closed interval is calculated instead:
#'
#' \deqn{[ \frac{m}{n} + \frac{L_m}{L_n}, \frac{m + a}{n} + \frac{L_{m+a}}{L_n} ],}
#' where \eqn{a} is the number of elements of the vector \eqn{x} that satisfies \eqn{x_i = \hat\mu}.
#'
#' @references
#' Damgaard, C., and J. Weiner. 2000. Describing Inequality in Plant Size or Fecundity. Ecology 81:1139–1142.
#'
#' @examples
#' x <- rbeta(1000,1,0.2)
#' hist(x)
#' plot(lorenz_curve(x))
#' abline(a=1,b=-1) # Line of symmetry
#' lac(x) # S = 0.48
#'
#' x2 <- rallo(1000,a=3)
#' hist(x2)
#' plot(lorenz_curve(x2))
#' abline(a=1,b=-1) # Line of symmetry
#' lac(x2) # S = 1.24

#'
#' @export
lac <- function(x, n = NULL, interval = FALSE, na.rm = FALSE) {
  if(!is.numeric(x)){
    stop("x must be a numeric vector")
  }
  if(na.rm){
    x <- x[!is.na(x)]
  }
  if(any(is.na(x)) || any(x < 0)){
    return(NA_real_)
  }
  if(!is.null(n)){
    x <- rep(x,n)
  }
  x <- sort(x)
  mu <- mean(x)
  xlow <- x < mu
  xeq <- x == mu
  m <- sum(xlow)
  n <- length(x)
  Lm <- sum(x[xlow])
  Ln <- sum(x)


  if (any(xeq)) {
    a <- sum(xeq)
    Lma <- sum(x[xlow | xeq])
    s <- c(m/n + Lm/Ln, (m + a)/n + Lma/Ln)
    if (!interval)
      s <- mean(s) #hmmm
  } else {
    xm <- max(x[xlow])
    xm1 <- min(x[!xlow])
    delta <- (mu - xm)/(xm1 - xm)
    s <- (m + delta)/n + (Lm + delta * xm1)/Ln
  }
  return(s)
}


#' @title Generate Summary Statistics For A Vector of Data
#' @description Get values for a set of statistical probes.
#' @param x a vector of numeric data
#' @param probes a vector of function names, character string, or functions that can be found with \code{match.fun()}. Default probes are \code{mean()}, \code{var()}, \code{cv()}, \code{Gini()}, \code{max()}, \code{min()}, \code{median()}, \code{Skew()}, \code{Kurt()}, \code{q25()}, \code{q75()}, \code{N()}, and \code{Lasym()}.
#' @param na.rm a logical value indicating whether to drop \code{NA} values
#' @return
#'  A named vector of numeric values
#' @examples
#' x <- runif(100,0.1,1)
#' probe_distribution(x)
#' probe_distribution(x, probes = "cv")
#' probe_distribution(x, probes = c("sd", "length","se"))
#' probe_distribution(x, probes = c(function(v){sd(v)/length(v)^0.5}))
#' @export
probe_distribution<-function(x, probes = c("mean","var","cv","Gini",
                                           "max","min","median", "Skew",
                                           "Kurt","q25","q75","N","lac"),
                             na.rm = FALSE){
  if(na.rm){
    x <- x[!is.na(x)]
  }

  probe.val<-simplify2array(lapply(probes, FUN = function(fun,x){
    match.fun(fun)(x)
  }, x = x))
  names(probe.val) <- probes

  return(probe.val)
}


#' @export
apply_probes<-function(data.list, probes = c("mean","var","cv","Gini",
                                             "max","min","median", "Skew",
                                             "Kurt","q25","q75","N","lac"),
                       trivial.rm = TRUE, na.rm = FALSE, add.id = TRUE){
  if(trivial.rm){
    data.list <- data.list[
      simplify2array(lapply(data.list, FUN = function(x){
        sum(x) > 0
      }))
    ] # All zeros gets thrown out
  }

  out <- simplify2array(lapply(data.list,function(x, na.rm, probes) {
    probe_distribution(x, na.rm = na.rm, probes = probes)
  }, na.rm = na.rm, probes = probes))

  if(length(probes) > 1){
    out <- t(out)
  }

  out <- as.data.frame(out,row.names = NULL)
  out$id<-names(data.list)

  colnames(out)[seq_along(probes)] <- do.call("c",lapply(probes,FUN = function(x){
    if(!is.character(x)){
      gsub(" ","",paste0(deparse(match.fun(probes)),collapse = ""))
    } else {
      x
    }
  }))

  if(
    add.id &&
    all(
      grepl("--",unique(out$id)),
      na.rm = TRUE)
  ){
    out$surveyID<-gsub("--.*","",out$id) # Add survey id if data.list is detected at the plant level
  }
  rownames(out) <- NULL

  return(out)
}


#' @export
plot_distributions<-function(data.list,type=c("ecdf","hist"),
                             by=0.1,ecdf.xlim=c(0,1),...){
  if(any(type%in%"ecdf")){
    for (i in seq_along(data.list)){
      data.list[[i]] %>%
        round(digits = floor(log10(by)*-1)) %>%
        ecdf() %>%
        plot(col=i, add = c(i>1),xlim = ecdf.xlim,...)
    }
  }
  if(any(type%in%"hist")){
    data.tally.list<-vector(mode="list",length=length(data.list))
    for (i in seq_along(data.list)){
      data.tally.list[[i]]<-graphics::hist(data.list[[i]],
                                 breaks = seq(0,1,by=by),
                                 plot = F)
      data.tally.list[[i]]<-data.frame(
        x=data.tally.list[[i]]$mids,
        y=data.tally.list[[i]]$counts,
        type=names(data.list)[i])
    }
    tally.data<-do.call("rbind",data.tally.list)
    g<-tally.data %>%
      dplyr::group_by(type) %>%
      dplyr::mutate(y=y/sum(y)) %>%
      ggplot2::ggplot(ggplot2::aes(x=x,y=(y)^0.1,group=type))+
      ggplot2::geom_col(ggplot2::aes(fill=type,color=NULL),alpha=0.5,position = "dodge")+
      ggplot2::geom_line(stat="smooth",ggplot2::aes(color=type),size=1)+
      ggplot2::labs(y="Pr(Prop. herbivory)^0.1",
           x="Prop. herbivory",
           color="Data",
           fill="Data")+
      ggplot2::theme_bw(base_size=15)
    g2<-tally.data %>%
      dplyr::group_by(type) %>%
      dplyr::mutate(y=y/sum(y)) %>%
      ggplot2::ggplot(ggplot2::aes(x=x,y=(y),group=type))+
      ggplot2::geom_col(ggplot2::aes(fill=type,color=NULL),alpha=0.5,position = "dodge")+
      ggplot2::labs(y="Pr(Prop. herbivory)",
           x="Prop. herbivory",
           color="Data",
           fill="Data")+
      ggplot2::theme_bw(base_size=15)
    return(ggpubr::ggarrange(g,g2,common.legend = T,legend = "top"))
  }
}

#' @export
compare_dist<-function(data.list = NULL, obs.index = 1, pred.index = 2,
                             test = c("ks","kl","ad"), digits = 2, kl.by = 0.1){
  test <- match.arg(test, several.ok = TRUE)
  obs.data<-data.list[[obs.index]]
  pred.data<-data.list[[pred.index]]

  if(any(test%in%"ks")){
    ks.out<-suppressWarnings(ks.test(
       round(obs.data, digits = digits),
       round(pred.data, digits = digits)
    ))
    ks.out<-c("ks.D" = unname(ks.out$statistic),"ks.P" = ks.out$p.value)
  } else {
    ks.out <- NULL
  }
  if(any(test%in%"kl")){
    kl.div<-suppressMessages(KL(
      x = rbind(
        density(obs.data,from = 0,to = 1,bw = kl.by)$y,
        density(pred.data,from = 0,to = 1,bw = kl.by)$y
      ),
      unit = "log",
      est.prob = "empirical"))
    kl.out<-c("kl.div"=unname(kl.div))
  } else {
    kl.out<-NULL
  }
  if(any(test%in%"ad")){
    ad.out<-ad.test(
      round(obs.data, digits = digits),
      round(pred.data, digits = digits))
    ad.out<-c("ad"= (ad.out$ad)[1,1],"ad.P" = (ad.out$ad)[1,3])
  } else {
    ad.out<-NULL
  }

  out<-list(
    "ks"=ks.out,
    "kl"=kl.out,
    "ad"=ad.out
  )
  return(out)
}


#' @export
get_dist_test_sim<-function(allo.fit.list,test = c("ks","kl","ad"),
                            nboot=1,n.sim=NULL,digits=2){
  test <- match.arg(test)
  out.list<-vector(mode="list",length=nboot)
  obs.out<-lapply(allo.fit.list,function(x){
    x$data
  })
  names(obs.out)<-do.call("c",
                          lapply(allo.fit.list,function(x){
                            x$id
                          }))

  n.rows<-length(allo.fit.list)
  test.list<-vector(mode = "list",length=n.rows)

  for (j in seq_len(nboot)){
    herb.data.bind<-c(
      list("obs" = obs.out),
      get_data_sim(
               allo.fit.list = allo.fit.list,
               n.sim = n.sim,
               digits = digits,
               return.obs = F)
      )
    for (i in seq_len(n.rows)){
      test.list[[i]]<-tryCatch(compare_dist(purrr::map(herb.data.bind,i),
                                                  obs.index = 1,
                                                  pred.index = 2,
                                                  test=test,
                                                  digits = digits)[[
                                                    match(test, c("ks","kl","ad"))
                                                  ]],
                               error = function(e) NA_real_)
      cat(j,"-",i,"\r","\r")
    }
    out.list[[j]]<-do.call("rbind",test.list)
  }
  out<-simplify2array(out.list)

  if(!any(is.null(names(obs.out)))){
    dimnames(out)[[1]] <- names(obs.out)
  }
  return(out)
}


#' @title Generate Neutral Predictions From Fitted Model
#' @description This is a convenience function that generates predicted distribution from a list of fitted models
#' @param allo.fit.list A list of "allo_herb_fit" objects.
#' @param n.sim A numeric value indicating number of points to draw from the neutral distribution for each "allo_herb_fit" object. If \code{NULL}, the original sample size of the data will be used.
#' @param digits A numeric value indicating the number of digits to round data values to
#' @param return.obs A logical value indicating whether to return observed data as a list in the out put along with the predicted draws
#' @param obs.raw A logical value indicating whether to round the observed data; ignored if \code{return.obs} is set to \code{FALSE}.
#' @param  new.param A vector of named values used to make new predictions. If a parameter is set to \code{NA}, the parameter value is extracted from the "allo_herb_fit" object.
#' @returns A list of list of vectors of numeric values.
#' @export
get_data_sim<-function(allo.fit.list,n.sim = NULL,digits=2,
                       return.obs = TRUE, obs.raw = FALSE,
                       new.param = c("mean.phi.T" = NA,
                                     "min.phi" = NA,
                                     "max.phi" = NA,
                                     "a" = NA)
                       ){
  if(any(
    unlist(
      lapply(allo.fit.list,function(x){
        !inherits(x,"allo_herb_fit")
        })
    )
  )){
    stop("allo.fit.list contains object of non-'allo_herb_fit' class.")
  }

    sim.fit.out<-lapply(allo.fit.list,function(x){
      predicted.draws<-round(
        ralloT(n = ifelse(is.null(n.sim),
                          length(x$data),
                          n.sim),
               mean.phi.T = .choose_new_theta_val(x,new.param,"mean.phi.T"),
               min.phi = .choose_new_theta_val(x,new.param,"min.phi"),
               max.phi = .choose_new_theta_val(x,new.param,"max.phi"),
               a = .choose_new_theta_val(x,new.param,"a")
        ),
        digits = digits)
      return(predicted.draws)
    })
    names(sim.fit.out)<-do.call("c",
                                lapply(allo.fit.list,function(x){
                                  x$id
                                }))
    if(return.obs){
      if(obs.raw){
        obs.out<-lapply(allo.fit.list,function(x){
          x$data
        })
      } else {
        obs.out<-lapply(allo.fit.list,function(x){
          round(x$data, digits = digits)
        })
      }
      names(obs.out) <- names(sim.fit.out)
      herb.data.bind<-c(
        list("obs"=obs.out),
        list("fitted.pred"=sim.fit.out))
    } else {
      herb.data.bind<-list("fitted.pred"=sim.fit.out)
    }
    return(herb.data.bind)
}




#' @title Calculate Lorenz Curve
#' @description Calculate Lorenz Curve for a vector of data with known frequencies
#' @param x a vector of non-negative numeric data
#' @param n a vector of frequencies the same length as \code{x}. If \code{NULL} (default), the frequencies are assumed to be 1 (as in a vector of untabulated raw data).
#' @param na.rm a logical value indicating whether to drop \code{NA} values
#' @param nboot an integer value indicating the number of bootstraps to perform to estimate the uncertainty (default to 100). If set below 0 or \code{FALSE}, no bootstraps will be performed. If set to \code{TRUE}, 100 bootstraps will be performed.
#' @param return.array a logical value indicating whether to save the bootstrapped Lorenz Curves. Default is \code{FALSE}. Must set to \code{TRUE} for spaghetti plot in \code{plot.lc()}.
#' @return
#' an object of class "lc" and "list" with three slots:
#'
#' \code{lc}: a data.frame of the exact Lorenz Curve calculated from the supplied data
#'
#' \code{lc_boot_summary}: a data.frame of the mean and standard error of the Lorenz Curve. If no bootstrap is performed, \code{NULL} is returned.
#'
#' \code{boot_array}: an array of the individual Lorenz Curve of each bootstrapped dataset. If \code{return.array = FALSE}, \code{NULL} is returned.
#' @export
lorenz_curve <- function(x, n = NULL, na.rm = FALSE, nboot = 100, return.array = FALSE){
  if(!is.numeric(x)){
    stop("x must be a numeric vector")
  }
  if(na.rm){
    x <- x[!is.na(x)]
  }
  if(any(is.na(x)) || any(x < 0)){
    return(NA_real_)
  }
  if(!is.null(n)){
    x <- rep(x,n)
  }
  x <- sort(x)
  n <- rep(1, length(x))

  p <- c(0,cumsum(n)/sum(n))
  mat.exact <- cbind("p" = p,"L" = c(0,cumsum(x)/sum(x)))

  if(nboot > 0){
    if(is.logical(nboot)){
      nboot <- 100
    }
    boot.array <- simplify2array(lapply(seq_len(nboot), function(i,x){
      x.boot <- sample(x,size = length(x), replace = TRUE)
      x.boot <- sort(x.boot)
      return(cbind("p" = p,"L" = c(0,cumsum(x.boot)/sum(x.boot))))
    }, x = x))

    Lc.boot.summary <- cbind(p, t(apply(boot.array[,2,],1,function(x){
      c(mean(x),sd(x))
    })))
    colnames(Lc.boot.summary) <- c("p","L.mean","L.se")
    Lc.boot.summary <- as.data.frame(Lc.boot.summary)
  } else {
    Lc.boot.summary <- NULL
    boot.array <- NULL
  }

  if(!return.array){
    boot.array <- NULL
  }


  out <- list("lc" = as.data.frame(mat.exact),
              "lc_boot_summary" = Lc.boot.summary,
              "boot_array" = boot.array)

  out <- structure(out,
                   class = c("lc","list"))
  return(out)
}


#' @title Plot Lorenz Curve Using Base Graphics
#' @description Plot Lorenz Curve for an object of class 'lc'. If available, uncertainty from the bootstrap simulation is also plotted either as 95% confidence intervals or as a spaghetti plot.
#' @param x an object of class 'lc' generated by \code{lorenz_curve()}
#' @param spaghetti if \code{TRUE}, a spaghetti plot of all simulations in \code{boot_array} is plotted. If an integer value is supplied, the number of curves up to the total number of simulations in \code{boot_array} will be plotted. If \code{FALSE} (default), the 95% confidence interval will be plotted instead.
#' @param main the title of the plot
#' @param xlab,ylab the x or y axes label
#' @param lc.lwd,identity.lwd the thickness of the Lorenz Curve or the identity line (representing perfect equality).
#' @param lc.col,identity.col the color of the Lorenz Curve or the identity line (representing perfect equality).
#' @param ... extra arguments passed to \code{plot.default()}
#' @return
#' \code{NULL}
#' @export
plot.lc <- function(x, spaghetti = FALSE, main = "Lorenz curve",
                    xlab = "Proportion of x", ylab = "Cumulative proprotion x",
                    lc.lwd = 2, lc.col = "red", identity.lwd = 3, identity.col = "black",
                    ...){
  if(!inherits(x,"lc")){
    stop("x must be of class 'lc'")
  }

  plot(x = x$lc$p, y = x$lc$p,
       type = "l", col = identity.col, lwd = identity.lwd,
       xlab = xlab, ylab = ylab, main = main)

  if(spaghetti > 0){
    if(is.null(x$boot_array)){
      stop("'boot_array' is empty. Set return.array = TRUE in lorenz_curve()")
    }
    nboot <- dim(x$boot_array)[3]
    if(is.logical(spaghetti)){
      spaghetti <- nboot
    }
    ndraws <- ifelse(spaghetti < nboot, spaghetti, nboot)

    for (i in seq_len(ndraws)){
      graphics::lines(x$boot_array[,"p",i],x$boot_array[,"L",i], col = "grey")
    }
    graphics::lines(x = x$lc$p, y = x$lc$L, col = lc.col,lwd = lc.lwd)
  } else {
    if(!is.null(x$lc_boot_summary)){
      graphics::lines(x = x$lc_boot_summary$p,
                      y = x$lc_boot_summary$L.mean + x$lc_boot_summary$L.se * 1.96,
            col = "grey", lty = "dashed", lwd = lc.lwd)
      graphics::lines(x = x$lc_boot_summary$p,
                      y = x$lc_boot_summary$L.mean - x$lc_boot_summary$L.se * 1.96,
            col = "grey", lty = "dashed", lwd = lc.lwd)
      graphics::lines(x = x$lc_boot_summary$p,
                      y = x$lc_boot_summary$L.mean,
            col = lc.col,lwd = lc.lwd)
    } else {
      graphics::lines(x = x$lc$p, y = x$lc$L, col = lc.col,lwd = lc.lwd)
    }
  }
}

