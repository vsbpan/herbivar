### Export functions from other packages

#' @export
ad.test <- function(...) {kSamples::ad.test(...)}

#' @export
hessian <- function(...) {numDeriv::hessian(...)}

#' @export
KL<-function(x,test.na=TRUE,unit="log",est.prob=NULL) {philentropy::KL(x,test.na=test.na,unit = unit,est.prob = est.prob)}

#' @export
makeCluster <- function(...) {parallel::makeCluster(...)}

#' @export
stopCluster<-function(...){parallel::stopCluster(...)}

#' @export
ess_bulk<-function(...) {rstan::ess_bulk(...)}

#' @export
ess_tail<-function(...) {rstan::ess_tail(...)}

#' @export
Rhat<-function(...) {rstan::Rhat(...)}

#' @export
Gini<-function(x, ...) {
  if(any(x < 0)){
    warning("Vector contains negative one or more values. Interpert with caution.")
  }
  ineq::Gini(x, ...)
}

#' @export
Lasym<-function(x, ...) {
  if(any(x < 0)){
    warning("Vector contains negative one or more values. Interpert with caution.")
  }
  ineq::Lasym(x, ...)
}

#' @export
Atkinson<-function(...) {ineq::Atkinson(...)}

#' @export
Theil<-function(...) {ineq::Theil(...)}


#' @export
ineq<-function(...) {ineq::ineq(...)}


#' @export
Kolm<-function(...) {ineq::Kolm(...)}

#' @export
RS<-function(...) {ineq::RS(...)}

#' @export
entropy<-function(...) {ineq::entropy(...)}


#' @export
Skew<-function(...) {DescTools::Skew(...)}

#' @export
Kurt<-function(...) {DescTools::Kurt(...)}

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
