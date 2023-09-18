library(herbivar)
f <- function(x){
  a<-ralloT(1000,x)
  b<-1-a
  c(Skew2(b),cv(b),mean(a))
}



out <-vapply(seq(0.01,1.2,by=0.05),FUN = f,FUN.VALUE = numeric(3))

t(out)

plot(t(out))
abline(a = -1, b = -1)

Skew2 <- function(x){
  mu <- mean(x)
  mean((x - mu)^3 / var(x)^1.5)
}
Kurt2 <- function(x){
  mu <- mean(x)
  mean((x - mu)^4 / var(x)^2)
}


file_path <- system.file("extdata/mock_leaf8.png",package="herbivar")
img <- load.image(file_path) %>%
  crop(empty.rm = "auto", cut_edge = 10) %>%
  add_px_size("dpi:300") %>%
  thin(1)

array2Image(img) %>% pliman::image_binary()

array2Image(img) %>%
  pliman::analyze_objects(index = "R",invert=TRUE,watershed = TRUE)


img2 <-array2Image(img)# %>% EBImage::otsu()

px1 <-threshold2(as.cimg(img[,,,1]),thr = 0.3)
px2 <-threshold2(as.cimg(img[,,,2]),thr = 0.3)

plot(px1)
plot(px2)

(px1 | !px2) %>% plot()


ind <- read.csv(file = system.file("indexes.csv", package = "pliman",
                                   mustWork = TRUE), header = T, sep = ";")



library(herbivar)

img <- image_example() %>% thin(3)

color_index <- function(img, index = "all",plot = TRUE){
  .is_inst("pliman",stop.if.false = TRUE,prompt = TRUE)
  ind <- read.csv(file = system.file("indexes.csv", package = "pliman",
                                     mustWork = TRUE), header = T, sep = ";")
  ind$Equation <- gsub("B","B(.)",gsub("G","G(.)",gsub("R","R(.)",ind$Equation)))
  max2 <- function(...){
    slice_eval(...,FUN = "max")
  }
  min2 <- function(...){
    slice_eval(...,FUN = "min")
  }
  ind$Equation <- gsub("min","min2",gsub("max","max2",ind$Equation))
  if(index == "all"){
    index <- ind$Index
  }
  iml<-lapply(seq_along(index), function(i,index){
    form <- ind[ind$Index == index[i],"Equation"]
    imeval(img,~eval(parse(text=form)))
  }, index = index)
  names(iml) <- index
  iml <- as.imlist(iml)
  if(plot){
    herbivar::plot.imlist(iml)
  }
  return(iml)
}

color_index(img)
img <- thin(img,100)

all.equal(slice_eval(img,FUN = "max"),slice_eval(R(img),B(img),G(img),FUN = "max"))

img2 <- array(c(img,img),dim = c(dim(img)[1],dim(img)[2],2,3))


a <- EBImage::abind(img2,img2)
slice_eval(img2,img2,FUN = "max")


px <-threshold2((img),thr="otsu")


debug(immask)
threshold2(grayscale(img),thr=0.3) %>% plot()

immask((img),threshold2(grayscale(img),thr=0.4),background = "hotpink") %>% plot()

immask2(grayscale(img),threshold2(grayscale(img),thr=0.3),background = 1) %>% plot()

plot(img)
img[]







get_lambda_p<-function(p){
  -log(1-p)
}

plot(get_lambda_p(seq(0,1,0.001)))


library(herbivar)
img <- load.image("C:/R Projects/misc analysis/sylvie_images/IMG_20220809_12r.png")
plot(img)


img <- thin(img,3) %>% crop(empty.rm = "auto",cut_edge = 10)
plot(img)
color_index(img)



immask(img,
       threshold2(
         color_index(img,index = "GR",plot = FALSE)[[1]],
         thr = "otsu"
         ) |
         !threshold2(
           color_index(img,index = "SCI",plot = FALSE)[[1]],
           thr = "99.5"
         )
       ) %>% plot()


library(herbivar)

img <- load.image("misc/test1.png")
plot(img2)
img2 <- load.image("misc/test2.png")

crop_leaf(img2) %>%
  imager::boundary() %>%
  crop(y=2500:2636) %>%
  plot()

crop_leaf(img2) %>%
  imager::boundary() %>%
  crop(y=0:636) %>%
  plot()


d<-imager::imeval(crop_leaf(img2), ~ !is.na(.)) %>%
  imager::boundary() %>%
  .[,,1,1] %>%
  which(arr.ind = TRUE) %>%
  dist(upper = FALSE,diag = FALSE)

head(d)



hist(d)


analyze_holes

img3 <-crop_leaf(img2)

img4 <- rbind(
  rep(NA,ncol(img3)+2),
  cbind(rep(NA,nrow(img3)),img3[,,1,1],rep(NA,nrow(img3))),
  rep(NA,ncol(img3)+2)
) %>% as.cimg()

px <-imager::imeval(img4, ~ !is.na(.)) %>%
  as.cimg()



plot(px)

iml <- imager::split_connected(px,high_connectivity = TRUE)
iml[[lapply(iml, function(x){
  sum(x)
}) %>% simplify2array() %>%
  which.max()]] %>%
  imager::boundary() %>%
  .[,,1,1] %>%
  which(arr.ind = TRUE) %>%
  dist() %>%
  c() %>%
  max()

crop_leaf(img2) %>% leaf_length()

img<-herbivar::image_example(type = "p")
library(herbivar)

object <- img

grow_edge(img,1) %>% plot()








library(herbivar)

x <- seq(0.00001,1.5,by=0.00001)

rallo(10000,0.00001,2,a = 0.5) %>% hist()


rallo(10000, 0.0000001,2)
herbivar:::.allometry.herb.quasi.sim



for (i in c(1, 3, 10, 30, 100, 300,1000, 3000)){

}

hist(ralloT(10000,0.3))

library(herbivar)
library(tidyverse)

a<-ralloT(100000,min.phi = 0.2,max.phi = 0.23,lambda = 2) %>% hist(nclass = 200)
data.frame("x" = a$mids, "y" = a$counts/sum(a$counts)) %>%
  ggplot(aes(x=x,y=y)) +
  geom_col(color = "navy", fill = "steelblue") +
  theme_bw(base_size = 15) +
  labs(x = expression("Proportion herbivory"~phi[T]), y = expression(P(Binned~phi[T])),
       subtitle = expression(lambda~"="~2~","~phi[m]~"="~0.2~","~phi[M]~"="~0.23))

f(50)


x<-rallo(10000,max.phi = 10)

x[x>1] <- 1
hist(x)

f(10,max.phi = 40,min.phi = 0.001,a = 1.55)


f <- function(i, min.phi = 0.00001, max.phi = 100,a=14/9){
  n.sim <- 10000
  lambda <- i
  k <- rpois(n.sim, lambda)
  phi_T <- vapply(X = k,
                  FUN = function(x) sum(
                    rallo(n = x,min.phi = min.phi, max.phi = max.phi,a=a)),
                  FUN.VALUE = numeric(1))
  phi_T[phi_T > 1] <- 1
  phi_T[phi_T < 0.005] <- 0
  phi_T %>% hist(main = i,nclass= 100)
  mean(phi_T)
  invisible(phi_T)
}

?hist

library(herbivar)

debug(herbivar:::dalloT)
debug(herbivar:::.dalloT.cond.k.fft.conv)

dalloT(c(0,0.5,1),lambda = 1)

#1-sum(dpois(0:50,1)%*%t(cond.prob.mat[-nrow(cond.prob.mat),]))




nll <- function(theta){
  -sum(dcb(
    data.vec,
    lambda = plogis(theta[1]),
    log = TRUE
  ))
}

nll(c(1,4))

debug(optim2)

optim(c(0.1), fn = nll, lower = c(-Inf), upper = c(Inf), method = "BFGS")

data.vec <- rbeta(100, 1,1)
data.vec <- rcb(100, 0.1)

fit_bite_size(rallo(1000),family = "all",min.phi = 0.005, method = "BFGS")


optim2(init = c(1,1),
       fn = function(theta){
         -sum(dbeta(x = data.vec,
                    shape1 = exp(theta[1]),
                    shape2 = exp(theta[2]),
                    log = T))
       }, method = method,
       hessian = T)





herbivar::adjust_prop










x<-ralloT(100,mean.phi.T = 0.1)
a$optim.param$method
a<-fit_allo_herb(x)
fit_generic_null(x)
a

debug(fit_bite_size)




fit_bite_size(rallo(1000,max.phi = 0.5),min.phi = 0.005)


x <- rallo(100,max.phi = 0.90)

exp(seq(0.001,5,by = 0.1)) %>%
  vapply(nll,numeric(1)) %>% plot()



optim2(init = exp(11),
       fn = nll,
       method = "L-BFGS-B",
       lower = exp(0.005),
       upper = exp(10))



nll<-function(theta){
  -sum(dallo(x,max.phi = exp(theta), log = TRUE))
}

.herb_data_check <- herbivar:::.herb_data_check

debug(fit_bite_size)



?quadprog::solve.QP()


x

log(x)

x <- rallo(1000,min.phi = 1,max.phi = 10,a = 14/9)
o<-order(x)
x2 <- x[o]
vals <- unique(x2)
rval <- cumsum(tabulate(match(x2,vals)))/length(x2)
plot(log(1-rval)~log(vals))
lines(log(pallo(vals, min.phi = 1, max.phi = 10, a = 14/9,lower.tail = FALSE)) ~ log(vals))

pallo(4,min.phi = 1, max.phi = 10, a = 14/9,lower.tail = FALSE)
palloT(10,mean.phi.T = 0.4,max.phi = 6)

pallo





D(quote(x/100),"x")

X <- c(1,2,3)

vcv <- matrix(c(rbeta(9,1,1)),ncol = 3)
vcv

trans <-  c(
  "theta" = function(x) {plogis(x)},
  "meanlog" = function(x) {x},
  "sdlog" = function(x) {exp(z)})

x <- c(1,2,3)
var <- c()

trans$sdlog

parse(text = match.fun(trans$sdlog))







x <- ralloT(100,lambda = 3)

fit <- fit_allo_herb(x)
fit2 <- fit_generic_null(x)
fit
fit2
debug(coef.generic_null_fit)

coef(fit2,se = TRUE)

coef(fit,se = TRUE)



survival_plot(rgamma(1000,1,1))
survival_plot(rbeta(1000,1,1))
survival_plot(rlnorm(1000,1,1))
survival_plot(rallo(1000,max.phi = 5))
survival_plot()




allo_probed <- lapply(seq_len(100), function(i){
  x <- ralloT(1000, lambda = 3, max.phi = 1)
  probe_distribution(x)
}) %>%
  do.call("rbind",.)

allo_probed2 <- lapply(seq_len(100), function(i){
  x <- ralloT(1000, lambda = 3, max.phi = 0.5)
  probe_distribution(x)
}) %>%
  do.call("rbind",.)


probed.data <- rbind(
  data.frame("distribution" = "max.phi = 1", allo_probed),
  data.frame("distribution" = "max.phi = 0.5", allo_probed2)
)

vegan::rda(probed.data[,-1]) %>%
  get_biplot(group = probed.data$distribution)


vegan::rda(probed.data[,-1] ~ probed.data$distribution)


get_data_sim(allo_fit)


get_data_sim



lapply(seq_len(n.rows), function(i) {
  out <- tryCatch(compare_dist(purrr::map(herb.data.bind,
                                   i), obs.index = 1, pred.index = 2, test = test,
                        digits = digits)[[match(test, c("ks", "kl",
                                                        "ad"))]],
           error = function(e) NA_real_)
  cat(j, "-", i, "\r", "\r")
  return(out)
})

a<-chisq.test(x = rpois(10,1), y = rpois(10,1))
a$p.value

suppressMessages(KL(
  x = rbind(
  p,
  p2
  ),
  unit = "log",
  est.prob = "empirical"))

suppressMessages(KL(
  x = rbind(
    obs.data,
    pred.data
  ),
  unit = "log",
  est.prob = "empirical"))


?density
chisq.test(
  rbind(density(obs.data,
                from = 0, to = 1, bw = kl.by)$y,
        table(pred.data, from = 0, to = 1, bw = kl.by)$y)
)
?chisq.test


p

fit_allo_herb(x,
              optim.vars = c("mean.phi.T","a"),
              init = c(10, 1),
              by = 0.0001,
              method = "Nelder-Mead",
              cores = 4)




plot_distributions(list(x,x2), type = "ecdf")

P <- 1:10/sum(1:10)
Q <- 20:29/sum(20:29)
x <- rbind(P,Q)

plot(P)
plot(Q)

sum(P)
sum(Q)

KL(x,)

compare_dist(list("obs" = rbeta(10000,1,1),
                       "pred" = rbeta(10000,1,1)),
             bin_size = 0.01)

debug(compare_dist)

library(herbivar)
x <- ralloT(10000, lambda = 3)
x2 <- ralloT(10000, lambda = 1)


par(mfrow=c(1,2))


par(cex = 0.9)
par(fig=c(0.25,1,0,1))
survival_plot(x[x>0])
par(fig=c(0,0.25,0,1), new = TRUE)
plot(c(mean(x == 0)) ~ 1,
     ylim = c(0,1), xaxt='n',
     ylab = "P(X=0)", xlab = "")



par(cex = 0.9)
par(fig=c(0.25,1,0,1), new = TRUE)
survival_plot(x2[x2>0], add=TRUE)
par(fig=c(0,0.25,0,1), new = TRUE)
points(c(mean(x2 == 0)) ~ 1,
     ylim = c(0,1), xaxt='n',
     ylab = "P(X=0)", xlab = "")







library(R1magic)


N <- 100
# Sparse components
K <- 4
#  Up to Measurements  > K LOG (N/K)
M <- 40
# Measurement Matrix (Random Sampling Sampling)
phi <- GaussianMatrix(N,M)
# R1magic generate random signal
xorg <- sparseSignal(N, K, nlev=1e-3)
y <- phi %*% xorg# generate measurement
T.mat <- diag(N)# Do identity transform
p <- matrix(0, N, 1)# initial guess
# R1magic Convex Minimization ! (unoptimized-parameter)
ll <- solveL1(phi, y, T.mat, p)
ll
x1 <- ll$estimate

plot(1:100, seq(0.011,1.1,0.011), type ="n",xlab="",ylab="")

title(main="Random Sparse Signal Recovery",
      xlab="Signal Component",ylab="Spike Value")
lines(1:100, xorg , col = "red")
lines(1:100, x1, col = "blue", cex = 1.5)



plot(y)


library(herbivar)
library(tidyverse)

saveRDS(z, "real_herbivory_data.rds")



d$x
z <- list("field" = d$x, "lab" = lab_herb)


d <- readRDS("misc/z.rds")
lab_herb <- readRDS("misc/x.rds")



white_herb <-d$x[d$morph == "w"]
purple_herb <-d$x[d$morph == "p"]


probes_d %>%
  group_by(source,type) %>%
  summarise_all(mean)



# Cohen, J. E., and M. Xu. 2015. Random sampling of skewed distributions implies Taylor’s power law of fluctuation scaling. Proceedings of the National Academy of Sciences 112:7749–7754.




TL_approx <- function(x, na.rm = FALSE){
  if(na.rm){
    x <- x[!is.na(x)]
  }
  gamma1 <- Skew(x, method = 1)
  CV <- cv(x, method = "standard")
  kappa <- Kurt(x, method = 1)
  b.hat <- gamma1 / CV
  a.log <- log(var(x)) - gamma1 / CV * log(mean(x))
  b.se <- sqrt((kappa -1 - gamma1^2) / ((length(x) - 2) * CV^2))
  return(c("b.hat" = b.hat, "b.se" = b.se, "b.lower" = b.hat - b.se * 1.96,
           "b.upper" = b.hat + b.se * 1.96, "a.log" = a.log))
}


test_TL_approx <- function(x){
  if(is.character(x)){
    x <- parse(text = x)
  } else {
    x <- substitute(x)
  }
  analytical_out <- TL_approx(eval(x))
  MC_out <- replicate(1000, eval(x)) %>%
    apply(2, function(z){
      c("v" = log(var(z)), "m" = log(mean(z)))
    }) %>%
    t() %>%
    as.data.frame() %>%
    lm() %>%
    coef()
  list("analytical" = analytical_out, "MC" = MC_out)
}

undebug(test_TL_approx)



test_TL_approx(ralloT(1000,0.05))




library(herbivar)


library(Rcpp)



plot_distributions(
 list(
    rnorm(100),
  rnorm(100),
   "z" = rnorm(100)
 ), xlim = c(-10,10), return_transformed = TRUE, by = 1
)







