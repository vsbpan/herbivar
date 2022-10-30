herbivar.version <- function(x){
  print(utils::packageDate("herbivar"))
  message("Still in development; use with caution!")
}

.onAttach <- function(...){
  conflicts <- .herbivar_conflicts(message = FALSE)
  conflict_msg <- ifelse(
    conflicts == "",
    "",
    paste0("\n \n","Conflicts: \n", conflicts)
  )
  heart <- ifelse(rbinom(1,1,0.1) == 1, " <3", "")
  packageStartupMessage("herbivar version date: ",utils::packageDate("herbivar"),"\n",
                        "Package still in development; use with caution!",heart, conflict_msg)
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

reinstall.herbivar <- function(package_path = "C:/R Projects/Package Building/herbivar_0.1.0.tar.gz"){
  load.herbivar <- "herbivar" %in% (.packages())
  detach.herbivar()
  install.packages(package_path, repos = NULL, type="source")
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




.herbivar_dependencies <- function (include_self = TRUE){
  raw <- paste(utils::packageDescription("herbivar")$Imports,
                utils::packageDescription("herbivar")$Depends,sep = ",")
  imports <- strsplit(raw, ",")[[1]]
  parsed <- gsub("^\\s+|\\s+$", "", imports)
  names <- vapply(strsplit(parsed, "\\s+"), "[[", 1, FUN.VALUE = character(1))
  if (include_self) {
    names <- c(names, "herbivar")
  }
  names
}

.herbivar_conflicts <- function (message = TRUE) {
  envs <- grep("^package:", search(), value = TRUE)
  envs <- purrr::set_names(envs)

  objs <- lapply(envs, function(env){
    x <- ls(pos = env)
    if (identical(env, "package:dplyr")) {
      x <- setdiff(x, c("intersect", "setdiff", "setequal",
                        "union"))
    }
    x
  })

  if (length(objs) == 0){
    objs <- objs
  } else {
    stacked <- utils::stack(objs)
    objs <- tapply(as.character(stacked$ind), stacked$values, list)
  }

  conflicts <- purrr::keep(objs, ~length(.x) > 1)
  tidy_names <- paste0("package:", .herbivar_dependencies())
  conflicts <- purrr::keep(conflicts, ~any(.x %in% tidy_names))
  conflict_funs <- purrr::imap(conflicts, .confirm_conflict)
  x <- purrr::compact(conflict_funs)
  if (length(x) == 0)
    return("")
  pkgs <- x %>% purrr::map(~gsub("^package:", "", .))
  others <- pkgs %>% purrr::map(`[`, -1)
  other_calls <- purrr::map2_chr(others, names(others), ~paste0(.x,
                                                                "::", .y, "()", collapse = ", "))
  winner <- pkgs %>% purrr::map_chr(1)
  funs <- format(paste0(winner, "::", paste0(names(x),
                                             "()")))
  bullets <- paste0(funs,
                    " masks ", other_calls, collapse = "\n")
  if(message){
    message(bullets)
  }
  invisible(bullets)
}

.confirm_conflict <- function (packages, name){
  objs <- packages %>% purrr::map(~get(name, pos = .)) %>%
    purrr::keep(is.function)
  if (length(objs) <= 1)
    return()
  objs <- objs[!duplicated(objs)]
  packages <- packages[!duplicated(packages)]
  if (length(objs) == 1)
    return()
  packages
}
