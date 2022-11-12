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
  if(any(x < 0)){
    warning("Vector contains negative one or more values. Interpert with caution.")
  }
  DescTools::Gini(x, ...)
}


#' @export
Skew<-function(x, weights = NULL, na.rm = FALSE, method = 1, ...) {
  DescTools::Skew(x, weights = NULL, na.rm = FALSE, method = 1, ...)
}

#' @export
Kurt<-function(x, weights = NULL, na.rm = FALSE, method = 1, ...) {
  DescTools::Kurt(x, weights = NULL, na.rm = FALSE, method = 1, ...) + 3
}

#' @export
clarkevans <- function(...){
  spatstat.core::clarkevans(...)
}

#' @export
clarkevans.test <- function(...){
  spatstat.core::clarkevans.test(...)
}

#' @export
hopskel <- function(...){
  spatstat.core::hopskel(...)
}

#' @export
allstats <- function(...){
  spatstat.core::allstats(...)
}

#' @export
Lest <- function(...){
  spatstat.core::Lest(...)
}

#' @export
Kest <- function(...){
  spatstat.core::Kest(...)
}

#' @export
Gest <- function(...){
  spatstat.core::Gest(...)
}

#' @export
Fest <- function(...){
  spatstat.core::Fest(...)
}

#' @export
Jest <- function(...){
  spatstat.core::Jest(...)
}

#' @export
envelope.ppp <- function(...){
  spatstat.core::envelope.ppp(...)
}

