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

J.index <- function(x=NULL,mean=NULL,var=NULL, na.rm = F){
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
cd <- function(x, robust = F, na.rm = F){
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


n <- function(x, na.rm = FALSE){
  if(na.rm){
    x <- x[!is.na(x)]
  }
  return(length(x))
}

q50 <- function(x, na.rm = FALSE){
  median(x, na.rm = na.rm)
}


q25 <- function(x, na.rm = FALSE){
  if(na.rm){
    x <- x[!is.na(x)]
  }
  return(as.numeric(quantile(x,probs = 0.25)))
}


q75 <- function(x, na.rm = FALSE){
  if(na.rm){
    x <- x[!is.na(x)]
  }
  return(as.numeric(quantile(x,probs = 0.75)))
}



#' @export
probe_distribution<-function(x, probes = c("mean","var","cv","Gini",
                                           "max","min","median", "Skew",
                                           "Kurt","q25","q75","n","Lasym"),
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
                                             "Kurt","q25","q75","n","Lasym"),
                       trivial.rm = TRUE, na.rm = FALSE, add.id = TRUE){
  out <- t(simplify2array(lapply(data.list,function(x, na.rm, probes) {
    probe_distribution(x, na.rm = na.rm, probes = probes)
  }, na.rm = na.rm, probes = probes)))
  out$id<-names(data.list)
  if(
    add.id &&
    all(
      grepl("--",unique(out$id)),
      na.rm = T)
  ){
    out$surveyID<-gsub("--.*","",out$id) # Add survey id if data.list is detected at the plant level
  }
  if(trivial.rm){
    out<-out[!out$mean==0,] # All zeros gets thrown out
  }
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
      data.tally.list[[i]]<-hist(data.list[[i]],
                                 breaks = seq(0,1,by=by),
                                 plot = F)
      data.tally.list[[i]]<-data.frame(
        x=data.tally.list[[i]]$mids,
        y=data.tally.list[[i]]$counts,
        type=names(data.list)[i])
    }
    tally.data<-do.call("rbind",data.tally.list)
    g<-tally.data %>%
      group_by(type) %>%
      mutate(y=y/sum(y)) %>%
      ggplot(aes(x=x,y=(y)^0.1,group=type))+
      geom_col(aes(fill=type,color=NULL),alpha=0.5,position = "dodge")+
      geom_line(stat="smooth",aes(color=type),size=1)+
      labs(y="Pr(Prop. herbivory)^0.1",
           x="Prop. herbivory",
           color="Data",
           fill="Data")+
      theme_bw(base_size=15)
    g2<-tally.data %>%
      group_by(type) %>%
      mutate(y=y/sum(y)) %>%
      ggplot(aes(x=x,y=(y),group=type))+
      geom_col(aes(fill=type,color=NULL),alpha=0.5,position = "dodge")+
      labs(y="Pr(Prop. herbivory)",
           x="Prop. herbivory",
           color="Data",
           fill="Data")+
      theme_bw(base_size=15)
    return(ggarrange(g,g2,common.legend = T,legend = "top"))
  }
}

#' @export
compare_dist<-function(data.list=NULL,obs.index=1,pred.index=2,
                             test=c("ks","kl","ad"),digits=1,kl.by=0.1){
  if(!is.null(data.list)){
    obs.data<-data.list[[obs.index]]
    pred.data<-data.list[[pred.index]]
  }

  if(any(test%in%"ks")){
    ks.out<-suppressWarnings(ks.test(
      obs.data %>% round(digits = digits),
      pred.data %>% round(digits = digits)
    ))
    ks.out<-c("ks.D" = unname(ks.out$statistic),"ks.P" = ks.out$p.value)
  } else {
    ks.out<-NULL
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
      obs.data %>% round(digits = digits),
      pred.data  %>% round(digits = digits))
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
                            n.boot=1,n.sim=NULL,digits=2){
  test <- test[1]
  out.list<-vector(mode="list",length=n.boot)
  obs.out<-lapply(allo.fit.list,function(x){
    x$data
  })
  names(obs.out)<-do.call("c",
                          lapply(allo.fit.list,function(x){
                            x$id
                          }))

  n.rows<-length(allo.fit.list)
  test.list<-vector(mode = "list",length=n.rows)

  for (j in seq_len(n.boot)){
    herb.data.bind<-c(
      list("obs"=obs.out),
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
                               error = function(e) NA)
      cat(j,"-",i,"\r")
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
#' @description generate neutral predictions from fitten model
#' This is a convenience function generates predicted distribution from a list of fitted models
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
      lapply(allo.fit.list,function(x){ !inherits(x,"allo_herb_fit")})
    )
  )){
    stop("allo.fit.list contains object of none 'allo_herb_fit' class.")
  }

    sim.fit.out<-lapply(allo.fit.list,function(x){
      predicted.draws<-round(ralloT(n = ifelse(is.null(n.sim),
                        length(x$data),
                        n.sim),
             mean.phi.T = ifelse(
               is.na(new.param["mean.phi.T"]),
               ifelse(
                 is.na(x$param.vals["mean.phi.T"]),
                 match.fun(x$param.val.trans[["mean.phi.T"]])(x$par[(x$theta.names == "mean.phi.T")]),
                 x$param.vals["mean.phi.T"]
               ),
               new.param["mean.phi.T"]
             ),
             min.phi = ifelse(
               is.na(new.param["min.phi"]),
               ifelse(
                 is.na(x$param.vals["min.phi"]),
                 match.fun(x$param.val.trans[["min.phi"]])(x$par[(x$theta.names == "min.phi")]),
                 x$param.vals["min.phi"]
               ),
               new.param["min.phi"]
             ),
             max.phi = ifelse(
               is.na(new.param["max.phi"]),
               ifelse(
                 is.na(x$param.vals["max.phi"]),
                 match.fun(x$param.val.trans[["max.phi"]])(x$par[(x$theta.names == "max.phi")]),
                 x$param.vals["max.phi"]
               ),
               new.param["max.phi"]
             ),
             a = ifelse(
               is.na(new.param["a"]),
               ifelse(
                 is.na(x$param.vals["a"]),
                 match.fun(x$param.val.trans[["a"]])(x$par[(x$theta.names == "a")]),
                 x$param.vals["a"]
               ),
               new.param["a"]
             )
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

