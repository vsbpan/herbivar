### Export functions from other packages

#' @export
ad.test <- function(...) {
  .is_inst("kSamples",stop.if.false = TRUE)
  kSamples::ad.test(...)
}

#' @export
hessian <- function(...) {
  .is_inst("numDeriv",stop.if.false = TRUE)
  numDeriv::hessian(...)
}

#' @export
KL<-function(x, test.na = TRUE, unit = "log", est.prob = NULL) {
  .is_inst("philentropy",stop.if.false = TRUE)
  philentropy::KL(x, test.na = test.na, unit = unit, est.prob = est.prob)
}

#' @export
makeCluster <- function(...) {parallel::makeCluster(...)}

#' @export
stopCluster<-function(...){parallel::stopCluster(...)}


#' @export
Gini<-function(x, ...) {
  if(any(is.na(x))){
    if(na.rm){
      x <- x[!is.na(x)]
    } else {
      return(NA_real_)
    }
  }
  if(any(x < 0)){
    warning("Vector contains negative one or more values. Interpert with caution.")
  }
  DescTools::Gini(x, ...)
}


#' @export
Skew<-function(x, weights = NULL, na.rm = FALSE, method = 1, ...) {
  DescTools::Skew(x, weights = weights, na.rm = na.rm, method = method, ...)
}

#' @export
Kurt<-function(x, weights = NULL, na.rm = FALSE, method = 1, ...) {
  DescTools::Kurt(x, weights = weights, na.rm = na.rm, method = method, ...) + 3
}

#' @export
clarkevans <- function(...){
  spatstat.explore::clarkevans(...)
}

#' @export
clarkevans.test <- function(...){
  spatstat.explore::clarkevans.test(...)
}

#' @export
hopskel <- function(...){
  spatstat.explore::hopskel(...)
}

#' @export
allstats <- function(...){
  spatstat.explore::allstats(...)
}

#' @export
Lest <- function(...){
  spatstat.explore::Lest(...)
}

#' @export
Kest <- function(...){
  spatstat.explore::Kest(...)
}

#' @export
Gest <- function(...){
  spatstat.explore::Gest(...)
}

#' @export
Fest <- function(...){
  spatstat.explore::Fest(...)
}

#' @export
Jest <- function(...){
  spatstat.explore::Jest(...)
}

#' @export
envelope.ppp <- function(...){
  spatstat.explore::envelope.ppp(...)
}

