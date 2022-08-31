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
  detach("package:herbivar",unload=T)
  message("goodbye!")
}

reinstall.herbivar <- function(...){
  detach.herbivar()
  install.packages("C:/R Projects/Package Building/herbivar_0.1.0.tar.gz", repos = NULL, type="source")
  .rs.restartR()
  library(herbivar)
}

### To do:
## Vignette for neutral model & leaf scan pipeline
## Documentation & more commentary
## Add leaf number support
## Add lambda argument in alloT functions
## Add in co-limitation theory
## Add in jiterate and aiterate
## Ad levene test, correlation least squares, Brown-Forsthe test, Bayesian test for heteroskedasticity (Dumitrscu 2019), dglm (?)
## Padé approximant for JE stuff


