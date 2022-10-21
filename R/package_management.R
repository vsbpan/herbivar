herbivar.version <- function(x){
  print(utils::packageDate("herbivar"))
  message("Still in development; use with caution!")
}

.onAttach <- function(...){
  heart <- ifelse(rbinom(1,1,0.1) == 1, " <3", "")
  packageStartupMessage("herbivar version date: ",utils::packageDate("herbivar"),"\n",
                        "Package still in development; use with caution!",heart)
}

detach.herbivar <- function(x){
  message("detatching package. . . ")
  if("herbivar" %in% (.packages())){
    detach("package:herbivar",unload=TRUE)
    message("goodbye!")
  } else {
    message("package is not loaded and therefore cannot be unloaded.")
  }
}

reinstall.herbivar <- function(...){
  load.herbivar <- "herbivar" %in% (.packages())
  detach.herbivar()
  install.packages("C:/R Projects/Package Building/herbivar_0.1.0.tar.gz", repos = NULL, type="source")
  .rs.restartR()
  if(load.herbivar){
    library(herbivar)
  }
}

### To do:
## Patch a = 1, a = 2, add documentation
## Vignette for neutral model & leaf scan pipeline
## Documentation & more commentary
## Add leaf number support
## Add in co-limitation theory
## Add levene test, correlation least squares, Brown-Forsthe test, Bayesian test for heteroskedasticity (Dumitrscu 2019), dglm (?)
## Kaiser–Meyer–Olkin Test




# f <- function(x){
#   a<-ralloT(1000,x)
#   c(Skew2(a)^2,-Kurt2(a))
# }
#
#
#
# out <-vapply(seq(0.01,1.2,by=0.05),FUN = f,FUN.VALUE = numeric(2))
#
#
# plot(t(out))
# abline(a = -1, b = -1)
#
# Skew2 <- function(x){
#   mu <- mean(x)
#   mean((x - mu)^3 / var(x)^1.5)
# }
# Kurt2 <- function(x){
#   mu <- mean(x)
#   mean((x - mu)^4 / var(x)^2)
# }
#
# Skew2(runif(1000))
