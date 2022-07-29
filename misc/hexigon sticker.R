library(hexSticker)
library(tidyverse)
library(magick)




# imgurl <- "C:/R Projects/Package Building/VarEcol/misc//white oak.png"
# a <- image_read(imgurl)
# b <- image_fill(a, "black", point = "+500 + 400", fuzz = 100)
# plot(b)
#
# magick::image_write(b,"C:/R Projects/Package Building/VarEcol/misc/white oak cleaned.png")
#img <- image_read("C:/R Projects/Package Building/VarEcol/misc/white oak cleaned.png")


# library(GauPro)
# x <- c(0,9.80377537198365,19.1958010196686,36.5704460535198,43.9040231285617+5,62.5593277160078+5,79.3195570120588+5,91.8250818038359+5,95.464633195661+5,100+5)
# y <- cos(x/5) * 1.5^(-(x-70)^2/100)
# y[6] <- 0.69
# y[7] <- -0.92
# plot(x,y)
# gp <- GauPro(x, y, parallel=FALSE)
# self <- gp
# minx <- min(self$X)
# maxx <- max(self$X)
# x1 <- minx - 0.1 * (maxx - minx)
# x2 <- maxx + 0.1 * (maxx - minx)
# x <- seq(x1, x2, length.out = nn)
# px <- self$pred(x, covmat = T)
# Sigma.try <- try(newy <- MASS::mvrnorm(n = n2, mu = px$mean,
#                                          Sigma = px$cov))
#
# miny <- min(newy)
# maxy <- max(newy)
#
# d <- as.data.frame(t(newy) %>% cbind(x))
# d <- gather(d,key="sim",value = "y", -x)
#
# p <- d %>% ggplot(aes(x=x,y=y,group=sim))+geom_line(size=0.3, color = "darkgrey")+theme_void()+ theme_transparent()
# p
# ggsave("misc/gp.png",p,dpi = 800)

img <- image_read("C:/R Projects/Package Building/herbivar/misc/herbivar_merged_image.png")

sticker(img, package="herbivar",
        p_size=20,
        s_x=1,
        s_y=0.58,
        s_width=2,
        s_height=2.1,
        h_fill = "cadetblue",
        h_color = "darkseagreen4",
        spotlight = T,
        l_y = 1.35,
        l_x = 1,
        l_width = 2.7,
        l_height = 20,
        l_alpha = 0.35,
        filename="misc/herbivar_sticker.png")


