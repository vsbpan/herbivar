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






