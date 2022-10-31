#' @title Convert Array to Class 'Image'
#' @description Convert array formatted as cimg or pixset to an object of class 'Image' from the package \code{EBImage}.
#' @param x A four dimensional array (as in the format in the package \code{imager}), a cimg, or pixset object.
#' @return an 'Image' object
array2Image <- function(x){
  if(dim(x)[3] > 1){
    warning("More than one depth layer detected. Choosing only the first one. ")
  }
  if(dim(x)[4] > 1){
    if(dim(x)[4] != 3){
      warning("Number of color channels should be 1 or 3 only.")
    }
    color <- "Color"
  } else {
    color <- "Grayscale"
  }
  out <- EBImage::Image(data = x[,,1,],
                        dim = dim(x)[-3],
                        colormode = color)
  attr(out,"px.size") <- attr(x,"px.size")
  return(out)
}

#' @title Convert Array to Class 'Image'
#' @description Convert an object of class 'Image' from the package \code{EBImage} to an array formatted as cimg or pixset.
#' @param x An 'Image' object
#' @param imager_class If \code{TRUE} (default), convert the array into a cimg or pixset depending on whether the array elements are numeric or logical.
#' @return A four dimensional array (as in the format in the package \code{imager}), a cimg, or pixset object.
Image2array <- function(x, imager_class = TRUE){
  spec <- dim(x)[3]
  spec <- ifelse(is.na(spec),1,spec)

  if(spec > 1){
    if(spec != 3){
      warning("Number of color channels should be 1 or 3 only.")
    }
  }
  out <- array(x[], dim = c(dim(x)[1:2],1,spec))
  if(imager_class){
    if(typeof(out) == "double"){
      out <- imager::as.cimg(out)
    } else if(typeof(out) == "logical"){
      out <- imager::as.pixset(out)
    }
  }
  attr(out,"px.size") <- attr(x,"px.size")
  return(out)
}


#' @title Mask Unselected Pixels
#' @description Set pixels in an image not selected in the pixset as an user defined value
#' @param object a cimg or array for which the values of unselected pixel are set
#' @param pixset a pixset object defining the set of pixels to be retained
#' @param background a numeric vector of equal length as the number of color spectra, a character string to be parsed by \code{col2rgb()}, or \code{NA} (default).
#' @return a masked object
#' @examples
#' img <- image_example()
#' px <- threshold2(img,thr = 0.3)
#' plot(px)
#'
#' immask(img,px,"hotpink") %>% plot()
#' immask(img, px) %>% plot()
#'
#' immask(grayscale(img),B(px), 0) %>% plot()
#'
immask <- function(object, pixset, background = NA){
  img.spec <- imager::spectrum(object)
  mask.spec <- imager::spectrum(pixset)
  if(is.character(background)){
    if(img.spec==3){
      background <- c(grDevices::col2rgb(background))/255
    } else {
      message("Background color taken as luminance")
      background <- crossprod((grDevices::col2rgb(background))/255, c(0.3,0.59,0.11))
    }

  }

  if(img.spec==3){
    if(mask.spec == 3){
      col.index <- c(1,2,3)
      if(length(background) < 3){
        background <- rep(background,3)
      }
    } else {
      col.index <- c(1,1,1)
    }
    object[,,,1][pixset[,,,col.index[1]]] <- background[1]
    object[,,,2][pixset[,,,col.index[2]]] <- background[2]
    object[,,,3][pixset[,,,col.index[3]]] <- background[3]
  } else {
    if(mask.spec > 1){
      warning("Only fist color channel of pixset used to subset from image.")
    }
    if(length(background)>1){
      warning("Only first elemenet of 'background' used. ")
    }
    object[,,,1][pixset[,,,1]] <- background[1]
  }

  return(object)
}


#' @title Evaluate Function For Each Height, Width, And Depth Element
#' @description Evaluate a function for each height, width, and depth element share among a set of provided objects or color spectrum arrays.
#' @param ... a set of array, cimg, or pixset
#' @param FUN a function to be evaluated, supported by \code{match.fun()}.
#' @return a cimg object
#' @examples
#' img <- image_example() %>% thin(3)
#'
#' # Find the maximum RGB value
#' out <- slice_eval(G(img),R(img),B(img),FUN = "max")
#' out
#' # same thing
#' all.equal(
#'   out,
#'   slice_eval(img,FUN = "max"))
#'
slice_eval <- function(..., FUN){
  a <- EBImage::abind(...)
  out <- array(
    apply(
      a,
      MARGIN = c(1,2,3),
      FUN = FUN
    ),
    dim =c(dim(a)[1:3],1)
  )
  as.cimg(out)
}

#' @title Invert colors
#' @description Invert colors of an image. Works with any object type that is supported by \code{imager::imeval()}.
#' @param object an image, pixset, or imlist
#' @return the same object with inverted colors
#' @export
invert <- function(object){
  imeval(object, ~max_scale(object)-.)
}




#' @title Convert split_herb To Point Pattern
#' @description converts a 'split_herb' object to a 'ppp' object. When \code{pt.type = "centroid"}, centroids that are rejected by \code{ppp()} because they are out of bound are nudged to the nearest within bound location. Usually, these points are on or close to the edge already.
#' @param split_herb an object of class 'split_herb'.
#' @param pt.type An atomic character indicating what to turn into a point. Acceptable values are "centroid" for centroids, and "px" for pixels.
#' @return an object of class 'ppp' with pixels as the unit
#' @export
split2ppp<-function(split_herb,pt.type = c("centroid","px")){
  if(!inherits(split_herb,"split_herb")){
    stop("Object not of class 'split_herb'")
  }
  pt.type <- match.arg(pt.type)

  if(pt.type == "centroid"){
    ppp.out<-spatstat.geom::ppp(x = split_herb$hole_centroid[,2],
                                y= split_herb$hole_centroid[,1],
                                mask = !is.na((split_herb$px)[,,1,1]),
                                unitname = "pixels")
    rejects<- attr(ppp.out,"rejects")
    if(!is.null(rejects)){
      index<-which(!is.na(split2mat(split_herb)),arr.ind = T)
      split_herb$hole_centroid[
        (split_herb$hole_centroid) %>%
          {(.[,2] %in% rejects$x &
              .[,1] %in% rejects$y)},] <- index[cbind(rejects$x,rejects$y) %>%
                                                  apply(
                                                    MARGIN = 1,
                                                    FUN = function(x){
                                                      which.min((index[,1] - x[1])^2 +
                                                                  (index[,2] - x[2])^2)}),]

      ppp.out<-spatstat.geom::ppp(x = split_herb$hole_centroid[,2],
                                  y= split_herb$hole_centroid[,1],
                                  mask = !is.na((split_herb$px)[,,1,1]),
                                  unitname = "pixels")
      rejects<- attr(ppp.out,"rejects")
      status <- is.null(rejects)
      warning(paste0("Nudging rejects to nearest valid pixel . . . ",
                     ifelse(status,"success!","failed!")),call. = F)
      # if(!status && ask.nudge){
      #   message("Select new coordinates for rejected points.")
      #   plot(split_herb$px)
      #   points(split_herb$hole_centroid[,2]~split_herb$hole_centroid[,1],
      #          col="blue",pch=19)
      #   points(rejects$y~rejects$x,col="red",pch=19)
      #   selected.pts <- locator(length(rejects$x))
      #   print(selected.pts)
      #   split_herb$hole_centroid[
      #     (split_herb$hole_centroid) %>%
      #       {(.[,2] %in% rejects$x &
      #           .[,1] %in% rejects$y)},] <- cbind(selected.pts$y,selected.pts$x)
      #   ppp.out<-spatstat.geom::ppp(x = split_herb$hole_centroid[,2],
      #                               y= split_herb$hole_centroid[,1],
      #                               mask = !is.na((split_herb$px)[,,1,1]), unitname = "pixels")
      #   status<- is.null(attr(ppp.out,"rejects"))
      #   warning(paste0("Selected points is valid . . . ",ifelse(status,"success!","failed!")),call. = F)
      # }
    }
  } else if(pt.type == "px"){
    ppp.out <- mat2ppp((split_herb$px)[,,1,1])
  }
  attr(ppp.out,"px.size") <- split_herb$px.size
  return(ppp.out)
}

#' @title Convert split_herb Into Matrix
#' @description converts split_herb into matrix
#' @param split_herb an object of class 'split_herb'
#' @return a matrix
#' @export
split2mat <- function(split_herb){
  if(!inherits(split_herb,"split_herb")){
    stop("Object not of class 'split_herb'")
  }
  mat.out <- split_herb$px[,,1,1]
  attr(mat.out,"px.size") <- split_herb$px.size
  return(mat.out)
}


#' @author Simon Barthelme
#' @note The code for the function obtained from Simon Barthelme's vignette \link{https://cran.r-project.org/web/packages/imager/vignettes/gettingstarted.html}. accessed date 2022-07-25.
get_centers <- function(im, thr="99%", approx.res = 10000, km = c(1L, 2L)){
  dt <- imhessian(im) %$%
    { xx*yy - xy^2 } %>%
    threshold2(thr, approx.res = approx.res, km = km) %>%
    label()
  as.data.frame(dt) %>% subset(value>0) %>% dplyr::group_by(value) %>% dplyr::summarise(mx=mean(x),my=mean(y))
}


#' @author Simon Barthelme
#' @note The code for the function obtained from Simon Barthelme's vignette \link{https://cran.r-project.org/web/packages/imager/vignettes/gettingstarted.html}. accessed date 2022-07-25.
hessdet <- function(im,scale=1){
  isoblur(im,scale) %>% imhessian %$% { scale^2*(xx*yy - xy^2) }
}


#' @title Attach Pixel Size Attribute
#' @description Attach "px.size" attribute to an object.
#' @param object any object
#' @param px.size a list, a character string, or NULL that is passed to \code{px_size_calc()}.
#' @return The original object with a "px.size" attribute
#' @export
add_px_size <- function(object,px.size){
  attr(object, "px.size") <- px_size_calc(object,px.size)
  return(object)
}

#' @title Calculate Pixel Size
#' @description Calculate pixel size from known parameters or image
#' @param object a matrix, array, cimg, or pixset that is optionally passed to measure image size. Default is NULL.
#' @param px.size A list of the correct format as would be returned by this function, \code{NA}, \code{NULL}, or a character string. see details.
#' @details
#' The character string must be in the general format "reference unit:value". Other special characters "=" and "-" are also acceptable to delineate "reference unit" from "value".
#'
#' The "reference unit" specifies the unit of "value". This can be "dpi", "ppi", "mm_x" or "cm_x". Here, "x" (must be separated by the character "_") is a placeholder for a number.
#'
#' The "value" specifies the image resolution (for dpi or ppi) or how many pixels ("mm_x" or "cm_x") are represented by the reference unit. "width" or "height" may be supplied in place of a number to indicate that the width or height of the image is to be measured from the supplied object instead.
#'
#' @return
#' If \code{px.size} is \code{NA} or \code{NULL}, \code{NULL} is returned. Otherwise, a list is returned.
#'
#' \code{size}: the length of the side of a single pixel
#'
#' \code{unit}: unit of size. Currently, only mm is supported.
#' @examples
#' # Image resolution is 150 dpi or 150 ppi
#' px_size_calc(NULL, px.size = "dpi:150")
#' px_size_calc(NULL, px.size = "ppi-150")
#'
#' # 1234 pixels is known to be 33 mm long.
#' px_size_calc(NULL, px.size = "mm_33:1234")
#'
#'
#' # image width (100 pixels) is 33 cm
#' px_size_calc(matrix(rnorm(10000),nrow = 100), px.size = "cm_33:width")
#' px_size_calc(NULL, px.size = "cm_33=100") # same thing
#'
#' @export
px_size_calc <- function(object = NULL,px.size){
  if(is.list(px.size)){
    return(px.size)
  } else {
    if(is.null(px.size) || is.na(px.size)){
      size <- NA
    } else {
      px.size.vec <- unlist(strsplit(px.size, ":|=|-"))
      if(length(px.size.vec) == 2){
        unit.vec <- unlist(strsplit(px.size.vec[1], "_| "))
        unit <- unit.vec[1]
        val <- px.size.vec[2]
        unit.val <- as.numeric(ifelse(length(unit.vec) == 2, unit.vec[2], NA))
      } else {
        if(length(px.size.vec) > 2){
          warning(paste0("Only first two splitted strings used: ",
                         "'",px.size.vec[1],"' ","'",px.size.vec[2]),"'")
        } else if(length(px.size.vec) == 1){
          stop("Missing character between unit and number. e.g. 'dpi:150', 'mm_2.7-593', 'cm_5=width'")
        }
      }

      if(val == "width"){
        if(is.null(object)){
          stop("No objected provided to measure width")
        }
        val <- nrow(object)
      } else if(val == "height"){
        if(is.null(object)){
          stop("No objected provided to measure height")
        }
        val <- ncol(object)
      } else {
        val <- as.numeric(val)
        if(!is.numeric(val) || is.na(val)){
          stop(paste("Invalid pixel number:",px.size.vec[2]))
        }
      }

      if(unit == "cm"){
        val <- val / 10
        unit <- "mm"
      }

      if(unit == "dpi" || unit == "ppi"){
        size <- 25.39999 / val
        if(length(unit.vec) > 1){
          warning(paste("Ignored unrecognized string"))
        }
      } else if(unit == "mm"){
        if(is.null(unit.val) || !is.numeric(unit.val)){
          stop("Invalid unit value")
        }
        size <- unit.val / val
      } else {
        stop(paste("Invalid unit:",unit))
      }
    }

    if(is.na(size)){
      return(NULL)
    } else{
      return(list(
        "size"= size,
        "unit" = "mm"))
    }}
}



#' @title Maximum Pixel Value of Image
#' @description Get maximum pixel value of image. Useful for doing arithmetics with image matrix.
#' @param object a cimg, pixset, matrix, or array
#' @return an atomic integer value
#' @export
max_scale <- function(object){
  val <- max(c(object),na.rm = TRUE)
  return(ifelse(val > 1, 3L, 1L))
}



#' @title Convert Matrix Into Point Pattern
#' @description convert matrix into point pattern. Non-leaf parts, remaining leaf parts, and eaten leaf parts of the image matrix must be indicated by \code{NA}s, 0, and any number greater than 0 respectively.
#' @param mat a matrix
#' @return a 'ppp' object with pixels as the unit.
mat2ppp <- function(mat){
  if(!inherits(mat, "matrix")){
    stop("Object must be a matrix")
  }
  indices<-as.data.frame(which(mat > 0,arr.ind = T))
  X <- spatstat.geom::ppp(x = indices$col,
                          y=indices$row,
                          mask = !is.na(mat),
                          unitname = "pixels")
  attr(X,"px.size") <- attr(mat,"px.size")
  return(X)
}


#' @title Reduce Pixel Density
#' @description reduce pixel density
#' @param object an object of class 'cimg', 'pixset', 'matrix', or 'array'
#' @param thin a numeric value indicating the square root of the factor by which to reduce pixel density. Default is 3, which reduces pixel density by 9 times.
#' @return an object of the same class
#' @export
thin <- function(object, thin = 3){
  px.size <- attr(object,"px.size")
  indices<-c(T,rep(F,(thin-1)))
  x<-rep(indices,floor(dim(object)[1]/thin))
  y<-rep(indices,floor(dim(object)[2]/thin))

  if(is.matrix(object)){
    out<-object[x,y]
  } else if(is.cimg(object) || is.pixset(object)){
    out<-object[x,y,,,drop=FALSE]
  }
  px.size$size <- px.size$size * thin

  if(length(px.size$size)==0){
    px.size <- NULL
  }
  attr(out,"px.size") <- px.size

  return(out)
}


#' @title Crop Image
#' @description crop image
#' @param object a 'cimg', 'pixset', or 'array' object to be cropped.
#' @param x,y a vector of integers indicating which row or column to keep. If set to \code{NULL}, all rows or clumns will be retained.
#' @param empty.rm white regions of the image as indicated by the supplied value will be automatically cropped out up to a width (set by \code{pad}) on the image boundaries. Only supported for grayscale images unless set to "auto". Ignored if set to \code{NULL} (default). Images cropped via \code{crop_leaf()} has white space represented by \code{NA}, otherwise white space is usually represented by 1. If set to "auto", the color of the white space will be automatically selected via \code{kmeans}. The function will throw a warning if the object contains \code{NA} but \code{empty.rm} is not set to \code{NULL} or \code{NA}.
#' @param pad the number of pixels on the image boundaries when using \code{empty.rm}. Ignored if \code{empty.rm} is set to \code{FALSE}.
#' @param cut_edge the number of pixels to cut off from all image boundaries. No pixels are removed if set to \code{NA} (default). Overwrites \code{pad}.
#' @param km an integer vector of length two. The first value selects the the boundary rank in increasing means identified by kmeans clustering (defaults to 1, selecting the boundary that separates the two clusters with the lowest mean). The second value selects the k number of clusters to form (defaults to 2).
#' @return a cropped object of the same supplied object class.
#' @export
crop <- function(object, x = NULL, y = NULL, empty.rm = NULL,
                 pad = 10, cut_edge = NA, km = c(2,3)){
  if(!(is.cimg(object) || is.pixset(object) || is.array(object))){
    stop("Unsupported object type. Object must be of class 'cimg', 'array', or 'pixset'")
  }

  n.x<- nrow(object)
  n.y<- ncol(object)

  if(!is.null(empty.rm)){
    if(!is.na(empty.rm) && tolower(empty.rm) != "auto"){
      if(!is.numeric(empty.rm)){
        stop("'empty.rm' must be a numeric value, 'auto', NA, or NULL. ")
      }
      if(imager::spectrum(object) > 1 ){
        stop("empty.rm not supported for multi-spectrum objects")
      }
    }
    if(!is.na(empty.rm) && any(is.na(object))){
      if(tolower(empty.rm) != "auto"){
        warning("NA values detected; empty.rm defaults to NA.")
      }
      empty.rm <- NA
    }

    if(!is.null(x) || !is.null(y)){
      warning("Ignored x and y settings when removing empty spaces")
    }

    if(!is.na(empty.rm) && tolower(empty.rm) == "auto"){
      gs_obj <- grayscale(object)
      white_space_thr <- cut_kmeans(c(gs_obj), km = km)
      x.range <- range(which(rowSums((gs_obj > white_space_thr)) != n.y))
      y.range <- range(which(colSums((gs_obj > white_space_thr)) != n.x))
    } else {
      if(is.na(empty.rm)){
        x.range <- range(which(rowSums(is.na(object)) != n.y))
        y.range <- range(which(colSums(is.na(object)) != n.x))
      } else {
        x.range <- range(which(rowSums((object == empty.rm)) != n.y))
        y.range <- range(which(colSums((object == empty.rm)) != n.x))
      }
    }

    x.range[1] <- ifelse((x.range[1] - pad) < 1, 1, x.range[1] - pad)
    x.range[2] <- ifelse((x.range[2] + pad) > n.x, n.x, x.range[2] + pad)
    x <- seq(x.range[1],x.range[2],by=1)

    y.range[1] <- ifelse((y.range[1] - pad) < 1, 1, y.range[1] - pad)
    y.range[2] <- ifelse((y.range[2] + pad) > n.y, n.y, y.range[2] + pad)
    y <- seq(y.range[1],y.range[2],by=1)
  } else {
    if(is.null(x)){
      x <- seq_len(n.x)
    }
    if(is.null(y)){
      y <- seq_len(n.y)
    }
  }
  if(!is.na(cut_edge) && is.numeric(cut_edge)){
    x <- x[-c(seq_len(cut_edge),seq(from = n.x-cut_edge, to = n.x, by = 1))]
    y <- y[-c(seq_len(cut_edge),seq(from = n.y-cut_edge, to = n.y, by = 1))]
  }
  out<-object[x,y,,,drop=FALSE]
  attr(out,"px.size") <- attr(object,"px.size")
  return(out)
}

#' @title Find Threshold For k Clusters Using k-Means
#' @description find threshold for k clusters using k-means
#' @param x a numeric vector. \code{NA}s are automatically removed.
#' @param km an integer vector of length two. The first value selects the the boundary rank in increasing means identified by kmeans clustering (defaults to 1, selecting the boundary that separates the two clusters with the lowest mean). The second value selects the k number of clusters to form (defaults to 2).
#' @return an atomic numeric value of boundary threshold
#' @note The code for the function is adapted from Simon Barthelme's \code{imager:::cut.kmeans()}. version 0.42.13.
cut_kmeans <- function (x, km = c(1L,2L)) {
  if(km[1] >= km[2]){
    stop("Chosen boundary must be less than the number of clusters")
  }
  if(km[2] < 2){
    stop("Number of clusters must be greater than two")
  }
  x <- x[!is.na(x)]

  km_out <- kmeans(x, km[2])
  index <- which(km_out$centers == sort(km_out$centers,decreasing = TRUE)[km[2]])
  thr <- max(x[km_out$cluster == index])

  return(thr)
}

#' @title Load Example Leaf Image
#' @description A convenience function for loading example image built in the package.
#' @param type a character string indicating the example image to be loaded
#' @return an object of class 'cimg'
image_example <- function(type = c("raw","processed")){
  type <- match.arg(type)
  file <- switch(type,
                 "raw" = "extdata/raw_scan.jpg",
                 "processed" = "extdata/mock_leaf8.png")
  file_path <- system.file(file,package = "herbivar")
  img <- imager::load.image(file_path)
  return(img)
}



#' @title Image Color Index From pliman
#' @description Return and display different transformed colored image, useful for thresholding. Indices are imported from the pliman package. Custom color index can be calculated using \code{imager::imeval()}.
#' @details This is an experimental function for internal use at the moment. Credit is due to the package Pliman.
#' @param img a cimg object with three color spectra
#' @param index a character string indicating the indices to use. If "all", all indices will be chosen
#' @param plot if \code{TRUE} (default), plot the transformed images
#' @return an imlist object
color_index <- function(img, index = "all",plot = TRUE){
  .is_inst("pliman",stop.if.false = TRUE,prompt = TRUE)
  stopifnot(imager::spectrum(img)==3)
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
  } else {
    if(!any(index %in% ind$Index)){
      stop("No valid index selected. Try: ", paste0(ind$Index, collapse = ", "))
    }
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

