#install.packages("C:/R Projects/Package Building/herbivar_0.1.0.tar.gz", repos = NULL, type="source")
#detach("package:herbivar",unload=T)
# Version date  "2022-08-07 11:12:31 EDT"
#Sys.time()
herbivar.version <- function(x){
  print("2022-08-07 11:12:31 EDT")
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

### To do:
## Add leaf number support
## Add lambda argument in alloT functions
## Add in co-limitation theory
## Add in Jensen's sim and aiterate


#### End ####




