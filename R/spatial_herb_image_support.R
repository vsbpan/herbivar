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

#' @title Compute Watershed Transform With Grayscale And Color Image Support
#' @description Internal function of \code{clean_leaf()}. A wrapper of \code{watershed()} and \code{threshold()} that handles both grayscale and color images. If available, the blue spectrum of the transformed image is used for thresholding, as blue is a rare color among plats.
#' @param seed an image
#' @param p priority map
#' @param spectrum number of spectrum.
#' @param thr a threshold, either numeric, or "auto", or string for quantities.
#' @param approx.res the number of pixels used in threshold approximation. Default is 10000.
#' @param km an integer vector of length two. The first value selects the the boundary rank in increasing means identified by kmeans clustering (defaults to 1, selecting the boundary that separates the two clusters with the lowest mean). The second value selects the k number of clusters to form (defaults to 2). Ignored when not auto-thresholding.
#' @return a pixset
#' @export
watershed2 <- function(seed, p, spectrum, thr = "auto", approx.res = 10000, km = c(1L, 2L)){
  ws <- watershed(seed,p)
  if(spectrum == 3){
    ws <- as.cimg(ws[,,,3]) # Blue gray scale works better for leaves because blue is rare
  } else if(spectrum > 1){
    ws <- as.cimg(ws[,,,1])
  }
  ws <- threshold2(ws, thr = thr, approx.res = approx.res, km = km,
                   thr.exact = FALSE, return.thr.only = FALSE)
  return(ws)
}


#' @title Threshold Foreground With Automatic Threshold Selection
#' @description Internal function of \code{clean_leaf()}. Basically acts as a wrapper for \code{threshold()} with added threshold reporting and error handling. When foreground threshold is selected, either manually or automatically, above 100%, the threshold is reduced by 10% recursively up to five times until foreground threshold is less than 100%. Automatically selected threshold by \code{threshold()}that is below 40% maybe too low and is increased by 30% automatically.
#' @param img a cimg object
#' @param fg.thresh a threshold for turning the cimg pixels into a selected pixset. Can be either numeric, a string, or "auto". When set to "auto", the threshold is automatically selected by the internal function \code{threshold2()}.
#' @param fg.adjust used to adjust the automatic threshold. If the auto-threshold is at k, the effective threshold will be at \code{fg.adjust} * k.
#' @param approx.res the number of pixels used in threshold approximation. Default is 10000.
#' @param km an integer vector of length two. The first value selects the the boundary rank in increasing means identified by kmeans clustering (defaults to 1, selecting the boundary that separates the two clusters with the lowest mean). The second value selects the k number of clusters to form (defaults to 2). Ignored when not auto-thresholding.
#' @param fg.thr an internal argument to pass calculated threshold to reduce redundant computation. Must be a numeric value.
#' @return a pixset with the selected pixels
#' @export
fg_threshold <- function(img, fg.thresh = "auto", fg.adjust = 1, approx.res = 10000,
                         km = c(1L, 2L), fg.thr = NULL){
  if(is.null(fg.thr)){
    fg.thr <- threshold2(img,thr = as.character(fg.thresh),adjust = fg.adjust,
                         approx.res = approx.res, km = km,
                         return.thr.only = TRUE, thr.exact = FALSE)
    if(fg.thr >= 1){
      for(i in seq_len(5)){
        fg.thr <- fg.thr * 0.9
        if(fg.thr < 1){
          break
        }
      }
    }
    message(paste0("Foregound threshold selected: ",signif(fg.thr*100,4),"%"))

    if(fg.thresh == "auto" && fg.thr < 0.4){
      message("Foreground threshold might be too low. Automatically increasing fg.adjust by 30%.")
      fg.thr <- fg.thr * 1.3
      message(paste0("New foregound threshold selected: ",signif(fg.thr*100,4),"%"))
    }
  }

  fg <- threshold2(img, thr = fg.thr, return.thr.only = FALSE, thr.exact = TRUE)

  return(fg)
}

#' @title Maximum Pixel Value of Image
#' @description Get maximum pixel value of image. Useful for doing arithmetics with image matrix.
#' @param object a cimg, pixset, matrix, or array
#' @return an atomic integer value
#' @export
max_scale <- function(object){
  val <- max(c(object),na.rm = T)
  return(ifelse(val > 1, 3L, 1L))
}


#' @title Select Pixels of Uneaten Leaf From Scanned Image
#' @description get a binary-ized image of a leaf scan, removing background artifacts and eaten sections of the leaf.
#' @param object a cimg, pixset, or matrix, or array.
#' @param fg.thresh a threshold for selecting the foreground in the watershed transform. Can be either numeric, a string, or "auto". When set to "auto", the threshold is automatically selected by the internal function \code{threshold2()}. Default is "auto".
#' @param bg.thresh a threshold for selecting the background in the watershed transform. Can be either numeric, a string, or "auto". When set to "auto", the threshold is automatically selected by the internal function \code{threshold2()}. Default is "10%".
#' @param fg.adjust used to adjust the automatic threshold. If the auto-threshold is at k, the effective threshold will be at \code{fg.adjust} * k.
#' @param blur.size a positive numeric value indicating the size of the median filter used to remove small speckles. Default is 2.
#' @param plot if \code{TRUE}, a plot of the original image and binary-ized image side-by-side will be plotted. Maybe slow.
#' @param save.plot.path the file path as a character string for saving the quality check plot. If set to \code{NA} (default), no plot will be saved.
#' @param save.plot.size an integer vector of length two defining the size of the quality check plot in the number of pixels. If set to "original", the size of the original image will be retained (same height and double the width).
#' @param return.cimg if \code{TRUE} (default), the function returns a 'cimg' object. Otherwise, the function returns a 'pixset'.
#' @param approx.res the number of pixels used in threshold approximation. Default is 10000.
#' @param km an integer vector of length two. The first value selects the the boundary rank in increasing means identified by kmeans clustering (defaults to 1, selecting the boundary that separates the two clusters with the lowest mean). The second value selects the k number of clusters to form (defaults to 2). Ignored when not auto-thresholding.
#' @param ws.thresh a threshold for selecting leaf pixels from a watershed transformed image. Can be either numeric, a string, or "auto". When set to "auto", the threshold is automatically selected by the internal function \code{threshold2()}. Default is "auto".
#' @param px.size value passed to \code{px_size_calc()}. When set to \code{NA} (default), the 'px.size' attribute is extracted from the supplied \code{object}.
#' @details
#'
#' @return A 'cimg' object if \code{return.cimg} is set to \code{TRUE}, otherwise, a 'pixset'.
#' @note The code for the watershed method is adapted from Simon Barthelme's vignette \link{https://cran.r-project.org/web/packages/imager/vignettes/pixsets.html}. accessed date 2022-07-25.
clean_leaf <- function(object, fg.thresh = "auto", bg.thresh = "10%",
                       fg.adjust = 1, blur.size = 2, plot = TRUE,
                       save.plot.path = NA, save.plot.size = "original",
                       return.cimg = TRUE, approx.res = 10000, km = c(1L, 2L),
                       ws.thresh = "auto", px.size = NA){
  if(is.na(px.size)){
    px.size <- attr(object,"px.size")
  }

  if(!(is.cimg(object) || is.pixset(object) || is.array(object) || is.matrix(object))){
    stop("Unsupported object type.'object' must be a cimg, pixset, array, or matrix.")
  } else {
    object <- as.cimg(object)
  }

  bg <- (!threshold2(object,thr = bg.thresh, approx.res = approx.res, km = km))
  fg <- fg_threshold(img,fg.thresh = fg.thresh,fg.adjust = fg.adjust,
                     approx.res = approx.res, km = km)
  fg.thr <- attr(fg, "threshold")

  seed <- bg+2*fg
  p <- 1/(1+enorm(imgradient(object,"xy")))
  spec <- imager::spectrum(object)

  ws <- watershed2(seed = seed, p = p, spectrum = spec, thr = ws.thresh,
                   approx.res = approx.res, km = km)

  for(i in seq_len(9)){
    if( sum(ws) / (nPix(ws) * max_scale(ws)) > 0.3 ){
      break
    } else {
      if(i < 6){
        #Try different seeds
        ws <- watershed2(seed = seed, p = p, spectrum = spec, thr = ws.thresh,
                         approx.res = approx.res, km = km)
      } else {
        fg.thr <- fg.thr * (1 + 0.1)

        if(fg.thr >= 1) {
          break
        }

        fg <- fg_threshold(object, fg.thr = fg.thr)
        seed <- bg+2*fg
        p <- 1/(1+enorm(imgradient(object,"xy")))
        ws <- watershed2(seed = seed, p = p, spectrum = spec, thr = ws.thresh,
                         approx.res = approx.res, km = km)
      }
    }
  }
  ws <- medianblur(ws,blur.size)

  if(!is.na(save.plot.path)){
    par.default <- par("mfrow","mar")
    if(save.plot.size == "original"){
      save.plot.size <- c(nrow(object)*2, ncol(object))
    }
    png(save.plot.path,width = save.plot.size[1],
        height = save.plot.size[2])
    par(mfrow=c(1,2),mar=c(4,2,0.5,0.5))
    plot(object)
    plot(ws)
    dev.off()
    par(par.default)
  }
  if(plot){
    par.default <- par("mfrow","mar")
    par(mfrow=c(1,2),mar=c(4,2,0.5,0.5))
    plot(object)
    plot(ws)
    par(par.default)
  }
  if(return.cimg){
    ws <- as.cimg(ws)
  }
  attr(ws,"px.size") <- px_size_calc(ws, px.size)
  invisible(ws)
}



#' @title Crop Out Non-Leaf Pixels
#' @description Crop out non-leaf pixels marked by the color blue as \code{NA}s. Non-leaf pixels need to be cropped out for further analyses.
#' @param object a cimg, pixset, or matrix, or array to be cropped.
#' @param as.mat if \code{TRUE}, the function returns a matrix of the pixel values of the \code{rgb.index} spectrum. If \code{FALSE} (default), the function returns a grayscale 'cimg' object.
#' @param invalid the value to replace non-leaf pixels as. Default is \code{NA}.
#' @param rgb.index an integer value \code{1L}, \code{2L}, or \code{3L}, indicating red, green, or blue colors respectively to be identified as non-leaf pixels.
#' @param px.size value passed to \code{px_size_calc()}. When set to \code{NA} (default), the 'px.size' attribute is extracted from the supplied \code{object}.
#'
#' @return a cropped 'cimg' object
#' @export
crop_leaf <- function(object, as.mat = FALSE, invalid = NA, rgb.index = 3L, px.size = NA){
  if(is.na(px.size)){
    px.size <- attr(object,"px.size")
  }

  if(!(is.cimg(object) || is.pixset(object) || is.array(object) || is.matrix(object))){
    stop("Unsupported object type.'object' must be a cimg, pixset, array, or matrix.")
  } else {
    object <- as.cimg(object)
  }

  if(imager::spectrum(object) == 1){ #Handle grayscale images
    rgb.index <- 1
    other.col.check <- TRUE
  } else {
    other.colors <- seq_len(3)[-rgb.index]
    other.col.check <- object[,,,other.colors[1]] < 0.2  & object[,,,other.colors[2]] < 0.2
  }
  object[object[,,,rgb.index] > 0.6 & other.col.check] <- invalid
  object[object[,,,rgb.index] < 0.95] <- 0
  object[object[,,,rgb.index] >= 0.95] <- 1
  if(as.mat){
    out <- object[,,,rgb.index]
  } else {
    if(spectrum(object) > 1){
      out<-grayscale(object)
    } else {
      out <- object
    }
  }
  attr(out,"px.size") <- px_size_calc(out, px.size)
  return(out)
}


#' @title Calculate Leaf Herbivory
#' @description calculate leaf herbivory in proportion, pixels, and mm2. Non-leaf parts, remaining leaf parts, and eaten leaf parts of the image matrix must be indicated by \code{NA}s, 0, and 1 (or 3) respectively.
#' @param object an object of type 'split_herb', 'matrix', 'array ', 'cimg', or 'pixset'.
#' @param type a vector of character string indicating what values to return. Partial string matching supported. "proportion" and "percent" returns relative proportion herbivory. "pixels" and "px" returns the number of pixels eaten. "mm2" and "real" returns the leaf area eaten in mm2.
#' @param px.size value passed to \code{px_size_calc()}. When set to \code{NA} (default), the 'px.size' attribute is extracted from the supplied \code{object}.
#' @return a vector of proportion, number of pixels, and milometer squared area of herbivory. The object must have a "px.size" attribute or the \code{px.size} argument must be defined to return values for mm2 leaf area removed.
#' @export
leaf_herb <- function(object, type = c("proportion","percent",
                                    "pixels","px","mm2","real"),
                      px.size = NA){
  if(inherits(object, "split_herb")){
    mat <- object$px
  } else {
    mat <- object
  }

  if(!(is.cimg(mat) || is.pixset(mat) || is.matrix(mat) || is.array(mat))){
    stop("Unsupported object type")
  }

  valid.type = c("proportion","percent",
                 "pixels","px","mm2","real")
  if(is.na(px.size)){
    px.size <- attr(mat,"px.size")
  }
  if("all"%in%type){
    type <- c("prop","px","real")
  }
  type <- valid.type[pmatch(type,valid.type)]

  if(is.array(mat) && any(dim(mat)[c(-1,-2)] > 1)){
    stop("Unexpected dimension present. Matrix should be n X m")
  }

  na.v <- is.na(c(mat))
  if(!any(na.v)){
    warning("No NA detected; make sure the boundaries are cropped")
  }

  crp.mat <- c(mat)[!na.v]
  max.scale <- max_scale(mat)

  if(any(c("proportion","percent") %in% type)){
    prop.herb <- mean(crp.mat) / max.scale
  } else {
    prop.herb <- NA
  }
  if(any(c("px","pixels") %in% type)){
    px.herb <- sum(crp.mat) / max.scale
  } else {
    px.herb <- NA
  }
  if(any(c("mm2","real") %in% type)){
    px.size_mm <- px_size_calc(mat, px.size)
    real.herb <- sum(crp.mat) / max.scale * (px.size_mm$size)^2
  } else{
    real.herb <- NA
  }
  out <- c("prop"=prop.herb,"px"=px.herb,"mm2"=real.herb)
  return(out[!is.na(out)])
}


#' @title Calculate Leaf Area
#' @description calculate leaf area pixels and mm2. Non-leaf parts of the image matrix must be indicated by \code{NA}s.
#' @param object an object of type 'split_herb', 'matrix', 'array ', 'cimg', or 'pixset'.
#' @param px.size value passed to \code{px_size_calc()}. When set to \code{NA} (default), the 'px.size' attribute is extracted from the supplied \code{object}.
#' @return a vector of the number of pixels and milometer squared area of leaf area. The object must have a "px.size" attribute or the \code{px.size} argument must be defined to return values for mm2 leaf area.
#' @export
leaf_area <- function(object, px.size = NA){
  if(is.na(px.size)){
    px.size <- attr(object,"px.size")
  }

  if(inherits(object, "split_herb")){
    mat <- object$px
  } else {
    mat <- object
  }

  if(is.array(mat) && any(dim(mat)[c(-1,-2)] > 1)){
    stop("Unexpected dimension present. Matrix should be n X m")
  }

  if(!(is.cimg(mat) || is.pixset(mat) || is.matrix(mat) || is.array(mat))){
    stop("Unsupported object type")
  }

  na.v <- is.na(c(mat))
  if(!any(na.v)){
    warning("No NA detected; make sure the boundaries are cropped")
  }

  px.tot <- sum(!na.v)

  px.size_mm <- px_size_calc(mat, px.size)
  real.tot <- px.tot * (px.size_mm$size)^2
  return(c("px"=px.tot,"mm2"=real.tot))
}


#' @title Convert Matrix Into Point Pattern
#' @description convert matrix into point pattern. Non-leaf parts, remaining leaf parts, and eaten leaf parts of the image matrix must be indicated by \code{NA}s, 0, and any number greater than 0 respectively.
#' @param mat a matrix
#' @return a 'ppp' object with pixels as the unit.
mat2ppp <- function(mat){
  if(inherits(mat, "matrix")){
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
  px.size$size <- px.size$size * thin^2

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
#' @param empty.rm white regions of the image as indicated by the supplied value will be automatically cropped out up to a width (set by \code{pad}) on the image boundaries. Only supported for grayscale images unless set to "auto". Ignored if set to \code{NULL} (default). Images cropped via \code{crop_leaf()} has white space represented by \code{NA}, otherwise white space is usually represented by 1. If set to "auto", the color of the white space will be automatically selected via \code{kmeans}. The function will throw a warning if the object contains \code{NA} but \code{empty.rm} is not set to \code{NULL} or \code{NA} .
#' @param pad the number of pixels on the image boundaries when using \code{empty.rm}. Ignored if \code{empty.rm} is set to \code{FALSE}.
#' @param cut_edge the number of pixels to cut off from all image boundaries. No pixels are removed if set to \code{NA} (default). Overwrites \code{pad}.
#' @return a cropped object of the same supplied object class.
#' @export
crop <- function(object, x = NULL, y = NULL, empty.rm = NULL, pad = 10, cut_edge = NA){
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

    if(!is.na(empty.rm) && tolower(empty.rm) == "auto"){
      gs_obj <- grayscale(object)
      white_space_thr <- cut_kmeans(c(gs_obj), km = c(2,3))
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

#' @title Threshold Image With Missing Values And Error Handling
#' @description Basically the same thing as \code{imager::threshold()}, but handles missing values in the cimg object, preserves "px.size" attribute, increases the approximation resolution, and handles error generated by k-means approximation.
#' @details see \code{?imager::threshold()}
#' @param im the image
#' @param thr a threshold, either numeric, or "auto", or string for quantities. If \code{thr.exact} is \code{TRUE}, the value must be numeric.
#' @param approx skip pixels when computing quantiles in large images (defaults to \code{TRUE})
#' @param adjust use to adjust the automatic threshold: if the auto-threshold is at k, effective threshold will be at adjust*k (default 1)
#' @param approx.res the number of pixels used in threshold approximation. Default is 10000.
#' @param km an integer vector of length two. The first value selects the the boundary rank in increasing means identified by kmeans clustering (defaults to 1, selecting the boundary that separates the two clusters with the lowest mean). The second value selects the k number of clusters to form (defaults to 2). Ignored when not auto-thresholding.
#' @param return.thr.only if \code{TRUE} (default to \code{FALSE}), return the calculated numeric threshold only, otherwise a selected pixset of the image will be returned.
#' @param thr.exact if \code{TRUE}, the exact threshold supplied by \code{thr} is used. Used internally to reduce redundant computation.
#' @return a pixet with select pixels if \code{return.thr.only = FLASE}, otherwise, a numeric value is returned
#' @note The code for much of the function is obtained from Simon Barthelme's \code{imager::threshold()}. version 0.42.13.
#' @export
threshold2<-function (im, thr = "auto", approx = TRUE, adjust = 1,
                      approx.res = 10000, km = c(1L, 2L),
                      return.thr.only = FALSE, thr.exact = FALSE){
  if(!thr.exact){
    # Calculate 'thr' if not using exact supplied value
    if (is.character(thr)) {
      if (nPix(im) > approx.res && approx) {
        v <- im[round(seq(1, nPix(im), l = approx.res))]
      }
      else {
        v <- im
      }
      if (tolower(thr) == "auto") {
        thr <- tryCatch(
          cut_kmeans(c(v), km = km) * adjust,
          error = function(e){
            cut_kmeans(c(im), km = km) * adjust
          })
      }
      else {
        if(adjust != 1){
          warning("'adjust' ignored when 'thr' not set to 'auto'")
        }
        regexp.num <- "\\d+(\\.\\d*)?|\\.\\d+"
        qt <- stringr::str_match(thr, regexp.num)[, 1] %>%
          as.numeric
        thr <- quantile(v, qt/100,na.rm = TRUE)
      }
    }
  }
  if(return.thr.only){
    return(thr)
  } else {
    a <- im > thr
    attr(a, "threshold") <- thr
    attr(a,"px.size") <- attr(im,"px.size")
    return(a)
  }
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



#' @title Parse Holes And Extract Within Leaf Herbivory Data
#' @description Split pixels into different herbivory regions ("holes") and measure hole and leaf level data. A pre-processing step for use with subsequent functions such as \code{split2ppp()}, \code{nn_dist()}, \code{pair_dist()}, and
#' @param object A cimg, pixset, matrix, or array that has been cropped via \code{crop_leaf()}
#' @param min_prop a numeric value between 0 and 1 indicating the minimum proportion herbivory threshold for inclusion. This gets rid of artifacts in the image that are too small to actually be real herbivory or biologically non-trivial herbivory. Default is 0.0001.
#' @param plot if \code{TRUE} (default is \code{FALSE}), each hole splitted from the image will be plotted side-by-side. If more than 20 holes are detected, the user will be prompted to decide whether to proceed with plotting. Can be slow if the number of holes and image resolution is high.
#' @param tolerance a numeric value between 0 and 1 determining if two pixels belong to the same region. Defaults to 0.1.
#' @param silent if \code{TRUE} (default is \code{FALSE}), suppress messages.
#' @param image_id a character string supplied to name the returned object for book keeping.
#' @return
#' an object of class 'split_herb' and 'list'.
#'
#' \code{n_holes}: The number of holes. Holes with proportion area less than \code{min_prop} are excluded.
#'
#' \code{hole_prop}: The proportion area of each hole. Returns NA if \code{n_hole} is 0.
#'
#' \code{leaf_area}: The leaf area in pixels and mm2.
#'
#' \code{is.margin}: A vector of logical values indicating whether each parsed hole is marginal or internal. Returns NA if \code{n_hole} is 0.
#'
#' \code{hole_perimeter}: A \code{n_hole} X 2 matrix of values of the perimeter of each parsed hole returned by the function \code{hole_perimeter()}. The \code{non_leaf_border} argument is set to \code{TRUE}. Returns NA if \code{n_hole} is 0.
#'
#' \code{hole_centroid}: A \code{n_hole} X 2 matrix of the the centroid coordinates of each parsed hole. Returns NA if \code{n_hole} is 0.
#'
#' \code{px}: A binary-ized pixset of the supplied \code{object} whose value is above 0.8.
#'
#' \code{min_prop}: The supplied \code{min_prop} argument
#'
#' \code{px.size}: The length in millimeter of a pixel. If unsuccessfully extracted from the \code{object} \code{px.size} attribute, only results in pixels will be returned for other function outputs.
#'
#' \code{imlist}: A list of images of class \code{imlist} of each parsed hole.
#'
#' \code{image_id}: The supplied \code{image_id} argument
#'
#' @export
analyze_holes <- function(object, min_prop = 10^-4, plot = FALSE,
                          tolerance = 0.1, silent = FALSE, image_id = NA){
  if(!any(is.na(c(object)))){
      warning("No NA detected; make sure the boundaries are cropped")
  }
  px.size <- attr(object,"px.size")

  if(is.array(object) || is.matrix(object)){
    object <- as.cimg(object)
  }

  if(is.cimg(object)){
    px <- object > 0.8
  } else {
    if(is.pixset(object)){
      px <- object
    } else {
      stop("Unsupported object type.")
    }
  }

  if(imager::spectrum(px) > 1){
    warning("Image not in grayscale. Converting to grayscale.")
    px <- grayscale(px)
  }
  attr(px,"px.size") <- px.size

  lab <- label(px, tolerance = tolerance) * as.cimg(px)
  split.set <- unique(lab) %>% {
    .[. > 0 & !is.na(.)]
  }
  imlist <- split.set %>% map_il(function(v) lab == v) %>%
    lapply(function(x){
      attr(x,"px.size") <- px.size
      return(as.cimg(x))
    }) %>% as.imlist(check=FALSE)

  hole_prop <- simplify2array(lapply(imlist,
                                     function(x) {
                                       leaf_herb(x,type="proportion")
                                     }))
  imlist <- imlist[hole_prop > min_prop]
  n_holes <- length(imlist)

  if(plot && n_holes > 20){
    plot <- ifelse(
      menu(c("Yes","No"),
           title = paste0(n_holes,
                          " holes detected; plot will be slow! Proceed with plot?")) == 1,
      TRUE,
      FALSE)
  }
  if(plot && n_holes > 0){
    imlist %>% map_il(function(x){
      x %>% colorise(.,
                     px = ~. > 0.99,
                     col = "violet")
    }) %>%
      spatstat.geom::plot.imlist(plotcommand="plot")
  }

  if(n_holes > 0){
    hole_prop <- hole_prop[hole_prop > min_prop]

    is_margin <- simplify2array(lapply(
      imlist,
      FUN = function(x){
        is.margin(x)
      }
    ))

    hole_peri <- do.call("rbind",
                         lapply(imlist,FUN = function(x) {
                           hole_perimeter(x,non_leaf_border = TRUE)
                         }))

    hole_cen <- do.call("rbind",
                        lapply(
                          imlist,
                          FUN = function(x) {
                            (colMeans(which(x[, , 1, 1] > 0.99, arr.ind = T)))
                          }
                        ))

  } else {
    hole_prop <- NA
    is_margin <- NA
    hole_peri <- NA
    hole_cen <- NA
  }

  out<-list(
    "n_holes" = n_holes,
    "hole_prop" =  hole_prop,
    "leaf_area" = leaf_area(px),
    "is.margin" = is_margin,
    "hole_perimeter" = hole_peri,
    "hole_centroid" = hole_cen,
    "px" = px,
    "min_prop" = min_prop,
    "px.size" = px.size,
    "imlist" = imlist,
    "image_id" = image_id
  )
  if(!silent){
    cat(paste0("Number of holes detected: ",out$n_holes,"\n"))
  }

  out <- structure(out,
                   class = c("split_herb", "list"))

  invisible(out)
}





#' @title Display A List of Images Using Base Graphics
#' @description Plot different regions of herbivory form an object of class 'split_herb'
#' @param x a 'split_herb' object
#' @param prompt if \code{TRUE} (default), display a prompt for user to select whether to proceed with plotting when the number of holes is above 20.
#' @param ... extra arguments passed to \code{spatstat.geom::plot.imlist()}
#' @return NULL
plot.split_herb <- function(x, prompt = TRUE, ...){
  if(!inherits(x, "split_herb")){
    stop("Object must of of class 'split_herb'.")
  }

  if(x$n_holes > 20 && prompt){
    plot <- ifelse(
      menu(c("Yes","No"),
           title = paste0(x$n_holes,
                          " holes detected; plot will be slow! Proceed with plot?")) == 1,
      TRUE,
      FALSE)
  }
  if(x$n_holes > 0){
    x$imlist %>% map_il(function(img){
      img %>% colorise(.,
                     px = ~. > 0.99,
                     col = "violet")
    }) %>%
      spatstat.geom::plot.imlist(plotcommand="plot", ...)
  }
}


nn_dist <- function(object, # expect a matrix with two columns or split_herb
                    alternative = c("two.sided","clustered", "regular",
                                    "greater","lesser"),
                    silent = FALSE){
  if(inherits(object, "split_herb")){
    r<-dist(object$hole_centroid) %>%
      as.matrix() %>% apply(.,1,function(x){
        min(x[x>0])
      })
    rho <- object$n_holes / object$leaf_area["px"]
    names(rho) <- NULL
    r_E_bar <- 0.5 / sqrt(rho)
    r_A_bar<- mean(r)
    R <- r_A_bar / r_E_bar
    SE_r_E_bar <- 0.26136 / sqrt(object$n_holes * rho)
    c <- (r_A_bar - r_E_bar) / SE_r_E_bar

    alternative <- alternative[1]
    if(alternative == "two.sided"){
      p <- pnorm(abs(c),lower.tail = F)*2
    } else if(alternative == "clustered" || alternative == "lesser"){
      p <- pnorm(c,lower.tail = T)
    } else if(alternative == "regular" || alternative == "greater"){
      p <- pnorm(c,lower.tail = F)
    } else {
      stop("Acceptable 'alternative' values are: ",
           paste0(c("two.sided","clustered", "regular","greater","lesser"),
                  collapse = ", "))
    }

    test_stat<-c("rho" = rho, "rA" = r_A_bar, "rE" = r_E_bar,
                 "SErE" = SE_r_E_bar, "R" = R, "c" = c, "P" = p)
    if(!silent){
      cat("Clark-Evans Nearest Neighbor Z-Test","\n","\n")
      print(signif(test_stat[seq_len(4)],4))
      cat("\n")
      print(signif(test_stat[5:7],4))
      message(ifelse(p < 0.05,
                     ifelse( R > 1,
                             "Significant regularity",
                             "Significant aggregation"
                     ),
                     "No deviation from randomness detected"
      ))
      message("No correction supported at present")
    }
    invisible(
      list("test_stat" = test_stat,
           "r" = r)
    )
  } else {
    if(!(is.matrix(object) || is.array(object) || is.data.frame(object))){
      stop("Object type not supported")
    } else {
      if((dim(object)[2] != 2))
        stop("Object needs to be a n x 2 matrix")
    }
    r<-dist(object) %>%
      as.matrix() %>% apply(.,1,function(x){
        min(x[x>0])
      })
    return(r)
  }
}

#' @title Find All Pairwise Distance
#' @description Get all pairwise distance between centroids
#' @return an n choose 2 by 3 matrix of pairwise distance, where n is the number of centroids.
#' @param object An object of class 'split_herb', or a n by 2 matrix of centroid position as a matrix, array or data.frame.
#' @export
pair_dist<-function(object){
  if(inherits(object,"split_herb")){
    object <- object$hole_centroid
  } else {
    if(!(is.matrix(object) || is.array(object) || is.data.frame(object))){
      stop("Object type not supported")
    } else {
      if((dim(object)[2] != 2))
        stop("Object needs to be a n x 2 matrix")
    }
  }
  cen.dist<- dist(object)
  out <- data.frame(t(combn(attr(cen.dist,"Size"),2)),as.vector(cen.dist))
  names(out) <- c("cent1","cent2","dist")
  return(out)
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
  pt.type <- pmatch(pt.type[1],c("centroid","px"))
  if(is.na(pt.type)){
    stop("No matched 'pt.type' option.")
  }
  if(pt.type == 1){
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
  } else if(pt.type == 2){
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


#' @title Check Margin
#' @description Check if the selected pixels (must have value 1) is on the leaf margin (i.e. bordering \code{NA}s).
#' @param mat a cimg, pixset, matrix, or array
#' @return an atomic logical value
#' @export
is.margin <- function(mat){
  if(!is.matrix(mat)){
    if(is.cimg(mat) || is.array(mat) || is.pixset(mat)){
      if(imager::spectrum(mat) != 1){
        mat <- grayscale(mat)
      }
      mat <- mat[,,1,1]
    } else {
      stop("Unsupported object type")
    }
  }

  max.scale <- max_scale(mat)

  if(max.scale > 1){
    warning("Maximum pixel value is not 1. Try thresholding the image.")
  }

  if(!any(is.na(c(mat)))){
    warning("No NA detected; make sure the boundaries are cropped")
  }

  mat1 <- mat[-1,]
  mat2 <- mat[-nrow(mat),]
  mat3 <- mat[,-1]
  mat4 <- mat[,-ncol(mat)]

  out <- any(as.vector(is.na(mat1) & (mat2 > 0.99)),na.rm = T) ||
    any(as.vector(is.na(mat2) & (mat1 > 0.99)),na.rm = T) ||
    any(as.vector(is.na(mat3) & (mat4 > 0.99)),na.rm = T) ||
    any(as.vector(is.na(mat4) & (mat3 > 0.99)),na.rm = T)
  return(out)
}


#' @title Calculate Perimeter
#' @description Find the perimeter of selected pixels (must have value 1) including or excluding non-leaf borders.
#' @param mat a cimg, pixset, matrix, or array
#' @param non_leaf_border if \code{TRUE} (default is \code{FALSE}), return perimeter that includes non-leaf borders (edges that separate eaten leaf area and non-leaf area).
#' @param px.size value passed to \code{px_size_calc()}. When set to \code{NA} (default), the 'px.size' attribute is extracted from the supplied \code{mat} object.
#' @return a vector of the perimeter in the number of pixels and millimeters. The object must have a "px.size" attribute or the \code{px.size} argument must be defined to return values for mm perimeter
#' @export
hole_perimeter <- function(mat, non_leaf_border = FALSE, px.size = NA){
  if(is.na(px.size)){
    px.size <- attr(mat,"px.size")
  }
  px.size_mm <- px_size_calc(mat, px.size)

  if(!is.matrix(mat)){
    if(is.cimg(mat) || is.array(mat) || is.pixset(mat)){
      if(imager::spectrum(mat) != 1){
        mat <- grayscale(mat)
      }
      mat <- (mat[,,1,1])
    } else {
      stop("Unsupported object type")
    }
  }
  max.scale <- max_scale(mat)
  mat_cimg <- as.cimg(mat)

  if(max.scale > 1){
    warning("Maximum pixel value is not 1. Try thresholding the image.")
  }

  if(!any(is.na(c(mat)))){
    warning("No NA detected; make sure the boundaries are cropped")
  }

  perim_total <- imager::imeval(mat_cimg,~.> 0.99 & !is.na(.)) %>%
    imager::boundary() %>%
    as.vector() %>%
    sum()

  if(!non_leaf_border){
    mat1 <- mat[-1,]
    mat2 <- mat[-nrow(mat),]
    mat3 <- mat[,-1]
    mat4 <- mat[,-ncol(mat)]

    border_NA_perim <- sum((((is.na(mat1) & (mat2 > 0.99)) | (is.na(mat2) & (mat1 > 0.99)))[,-1]) | (((is.na(mat3) & (mat4 > 0.99)) | (is.na(mat4) & (mat3 > 0.99)))[-1,]),na.rm = T)

    perim_out <- perim_total - border_NA_perim
  } else {
    perim_out <- perim_total
  }
  return(c("px" = perim_out, "mm2" = perim_out * px.size_mm$size))
}


#' @title Print Value
#' @description print values from a 'split_herb' object
#' @param x a 'split_herb' object
#' @param ... additional arguments
#' @param digits number of digits to display. Default is four.
#' @return a vector of numeric values
#' @export
print.split_herb <- function(x, ..., digits = 4){
    cat("Image ID = ", x$image_id ,"\n", sep = "")
    cat("Minimum hole detection threshold: ", x$min_prop, "\n", sep = "")
    cat("\n")
    print(signif(c("total_holes" = x$n_holes, "toal_margin" = sum(x$is.margin), "prop_margin" = sum(x$is.margin)/x$n_holes),digits))
    cat("\n")
    if(!is.null(x$px.size)){
      print(
        signif(
          c("prop_herb" = sum(x$hole_prop),
            "leaf_area_mm2"= unname(x$leaf_area["mm2"]),
            "tot_herb_perim_mm2" = tryCatch(sum(x$hole_perimeter[,"mm2"]),
                                           error = function(e){
                                             NA
                                           })
          ),
          digits))
    } else {
      print(
        signif(
          c("prop_herb" = sum(x$hole_prop),
            "leaf_area_px"= unname(x$leaf_area["px"]),
            "tot_herb_perim_px" = tryCatch(sum(x$hole_perimeter[,"px"]),
                                           error = function(e){
                                             NA
                                           })
            ),
          digits))
    }
}

#' @title Invert colors
#' @description Invert colors of an image. Works with any object type that is supported by \code{imager::imeval()}.
#' @param object an image, pixset, or imlist
#' @return the same object with inverted colors
#' @export
invert <- function(object){
  imeval(object, ~max_scale(object)-.)
}


#' @title Plot one dimensional image
#' @description plot one dimensional image with pixel value plotted as a line. Internal function of \code{plot.cimg()}.
#' @param x the image
#' @param ... additional arguments passed to \code{plot.default()}
#' @return NULL
#' @note The code for the function is obtained from Simon Barthelme's \code{imager:::plot.singleton()}. version 0.42.13.
#' @export
plot.singleton <- function (x, ...){
  varying <- if (width(x) == 1)
    "y"
  else "x"
  l <- max(dim(x)[1:2])
  if (imager::spectrum(x) == 1) {
    plot(1:l, as.vector(x), xlab = varying, ylab = "Pixel value",
         type = "l", ...)
  }
  else if (imager::spectrum(x) == 3) {
    ylim <- range(x)
    plot(1:l, 1:l, type = "n", xlab = varying, ylim = ylim,
         ylab = "Pixel value", ...)
    cols <- c("red", "green", "blue")
    for (i in 1:3) {
      graphics::lines(1:l, as.vector(channel(x, i)), type = "l",
                      col = cols[i])
    }
  }
  else {
    stop("Unsupported image format")
  }
}


#' @title Plot Image Using Base Graphics Using Missing Value Handling
#' @description plot image the same way as \code{imager::plot.cimg()}. The only difference is that the \code{rescale} is automatically set to \code{FALSE} when the image contains missing pixel values. see \code{?imager::plot.cim()} for more details.
#' @param x the image
#' @param frame which frame to display, if the image has depth > 1
#' @param xlim x plot limits (default: 1 to width)
#' @param ylim y plot limits (default: 1 to height)
#' @param ylab x axis label
#' @param ylab y axis label
#' @param rescale rescale pixel values so that their range is [0,1]
#' @param colourscale,colorscale an optional colour scale (default is gray or rgb)
#' @param interpolate should the image be plotted with antialiasing (default \code{TRUE})
#' @param axes whether to draw axes (default \code{TRUE})
#' @param main main title
#' @param xaxs,yaxs The style of axis interval calculation to be used for the axes. See \code{?par}
#' @param col.na which colour to use for \code{NA} values, as R rgb code. The default is "rgb(0,0,0,0)", which corresponds to a fully transparent colour.
#' @param ... additional arguments passed to \code{plot.default()}
#' @return NULL
#' @note The code for the function is obtained from Simon Barthelme's \code{imager:::plot.singleton()}. version 0.42.13.
#' @export
plot.cimg <- function(x, frame, xlim = c(1, width(x)),
                      ylim = c(height(x), 1), xlab = "x", ylab = "y",
                      rescale = TRUE, colourscale = NULL,
                      colorscale = NULL, interpolate = TRUE, axes = TRUE, main = "",
                      xaxs = "i", yaxs = "i", asp = 1, col.na = rgb(0, 0, 0, 0),
                      ...) {
  v <- unique(c(x))
  if(length(v[!is.na(v)]) < 2){
    rescale <- FALSE
  }

  if (nPix(x) == 0)
    stop("Empty image")
  im <- x
  if (depth(im) > 1) {
    if (missing(frame)) {
      warning("Showing first frame")
      frame <- 1
    }
    im <- frame(x, frame)
  }
  if (1 %in% dim(im)[1:2]) {
    plot.singleton(im, ...)
  }
  else {
    if (is.character(asp) && asp == "varying") {
      plot(1, 1, xlim = xlim, ylim = ylim, xlab = xlab,
           ylab = ylab, type = "n", xaxs = xaxs, yaxs = yaxs,
           axes = axes, ...)
      as.raster(im, rescale = rescale, colorscale = colorscale,
                colourscale = colourscale, col.na = col.na) %>%
        rasterImage(1, height(im), width(im), 1, interpolate = interpolate)
      title(main)
    }
    else if (is.numeric(asp)) {
      plot.new()
      plot.window(xlim = xlim, ylim = ylim, asp = asp,
                  xaxs = xaxs, yaxs = yaxs, ...)
      rst <- as.raster(im, rescale = rescale, colorscale = colorscale,
                       colourscale = colourscale, col.na = col.na)
      rasterImage(rst, 1, nrow(rst), ncol(rst), 1, interpolate = interpolate)
      title(main)
      if (axes) {
        axis(1)
        axis(2)
      }
    }
    else {
      stop("Invalid value for parameter asp")
    }
  }
  invisible(x)
}


# ss <- function(k, mat.dist, method = "silhouette", nboot = 10){
#   ss <- seq_len(nboot)
#
#   for (i in seq_len(nboot)){
#     km <- kmeans(mat.dist, k)
#     if(method == "silhouette"){
#       ss[i] <- mean(silhouette(km$cluster,mat.dist)[,3])
#     } else if(method == "elbo"){
#       ss[i] <- km$tot.withinss
#     }
#   }
#   return(mean(ss))
# }
#
# ss_plot <- function(mat, k, method = "silhouette", plot = TRUE, nboot = 10){
#   mat.dist <- dist(mat)
#   ss.v<- vapply(k, FUN = function(k,mat.dist, method, nboot){
#     ss(k,mat.dist, method = method, nboot = nboot)
#   }, FUN.VALUE = numeric(1),
#   mat.dist = mat.dist,
#   method = method,
#   nboot = nboot)
#   if(method == "elbo"){
#     ss.v <- cumsum(ss.v)/k
#   }
#   d <- data.frame("k" = k, "ss" = ss.v)
#   if(plot){
#     plot(d$ss~d$k)
#   }
#   invisible(d)
# }

