
#' @title Remove punctuation and change all characters to lower case
#' @description
#' Remove punctuation and change all characters to lower case in supplied character vector
#' @param x a character vector
#' @return a character vector
clean_name <- function(x){
  tolower(gsub("[[:punct:]]| ","",x))
}

#' @title Validate file extension
#' @description
#' Check if the supplied file names have the correct file extensions
#' @param x a vector of file names
#' @param check_if a vector of acceptable file types. Defaults are \code{c("jpg", "jpeg")}
#' @param rm If \code{TRUE} (default is \code{FALSE}), remove file names that do not match the acceptable file extensions.
#' @param stop If \code{TRUE} (default), throw an error when unexpected file types are encountered. Otherwise, a warning is returned.
#' @return a vector of characters
validate_file_type <- function(x, check_if = c("jpg","jpeg"), rm = FALSE, stop = TRUE){
  out <- gsub(".*\\.","",gsub(".*/","",x))
  if(!isFALSE(check_if) && !all(out %in% check_if)){

    unexpected_file_types <- unique(out[!out %in% check_if])
    msg <- paste0("Unexpected file type detected: ", paste0(unexpected_file_types, collapse = ", "))

    if(stop){
      stop(msg)
    } else {
      warning(msg)
    }

    if(rm){
      out <- out[out %in% check_if]
    }
  }

  return(out)
}

#' @title Download images from URL
#' @description
#' Download images from a named vector of URL character string.
#' @param x a named vector with URLs as the values and file name (without the extension) as the name of the values.
#' @param save_dir a path to where the files will be saved
#' @param enforce_file_extension File extension to be enforced (default is \code{"jpg"}). Encountering unexpected file types will result in an error.
#' @param skip a vector of character string of file names in \code{x} that should be skipped.
#' @param confirm if \code{TRUE} (default), a prompt will appear to check if the user wants to continue with the download.

get_images <- function(x, save_dir,
                       enforce_file_extension = "jpg",
                       confirm = TRUE,
                       skip = NULL){

  if(length(unique(names(x))) != length(x)){
    stop("Missing unique names for each element of 'x'.")
  }

  if(confirm){
    choice <- utils::menu(
      choices = c("yes","no"),
      title = paste0("You are attempting to download ",
                     length(x),
                     " images.", " Continue downloading?"))
    if(!choice == 1){
      message("\nExit download")
      return(invisible(NULL))
    }
  }

  if(!is.null(skip)){
    skip <- names(x) %in% skip
    message(paste0("\n",sum(skip), " whitelisted images skipped."))
    x <- x[!skip]
  }

  n <- length(x)

  if(!is.character(enforce_file_extension)){
    if(enforce_file_extension == "jpg"){
      check_if <- c("jpg", "jpeg")
    } else {
      check_if <- enforce_file_extension
    }
    validate_file_type(x, check_if = check_if)
  }
  extn <- rep(enforce_file_extension, n)

  file_names <- paste0(names(x),".",extn)
  path_names <- paste(save_dir, file_names, sep = "/")
  path_names <- gsub("//", "/", path_names)

  existing_files <-list.files(save_dir)
  existing_files_paths <- paste(save_dir, existing_files, sep = "/")
  existing_files_paths <- gsub("//", "/", existing_files_paths)

  conflicts <- path_names %in% existing_files_paths
  n_conflicts <- sum(conflicts)
  if(n_conflicts > 0){
    choice_conflicts <- utils::menu(
      choices = c("overwrite", "skip", "exit"),
      title = paste0("There are ",
                     n_conflicts,
                     " files in your directory with the same name. Select action for all files."))
    if(choice_conflicts == 2){
      x <- x[!conflicts]
      path_names <- path_names[!conflicts]
      n <- length(x)
    }
    if(choice_conflicts == 3){
      message("\nExit download")
      return(invisible(NULL))
    }

  }

  for (i in seq_len(n)){
    try(utils::download.file(x[i], path_names[i], mode = "wb", quiet = TRUE))
    cat("Downloading", i, "out of", n,"                                \r")
  }
  message("\nDownload complete!")
}

#' @title Add quotes around a string
#' @description Add single quotes around a string
#' @param x a character string
#' @return a character string
add_quote <- function(x){
  paste0("'",x,"'")
}

#' @title Move files between directories
#' @description
#' Basically a wrapper for \code{file.rename()}
#' @param from_dir the source directory
#' @param to_dir the destination directory
#' @param file_names a character vector of the file names of the files to be moved
#' @return \code{NULL}
move_files <- function(from_dir, to_dir, file_names){
  for (i in seq_along(file_names)){
    file.rename(
      paste(from_dir,file_names[i], sep = "/"),
      paste(to_dir,file_names[i], sep = "/")
    )
  }
}

#' @title Delete image files from a directory
#' @description
#' Basically a wrapper for \code{file.remove()}
#' @param dir the directory for which image files are searched
#' @param ask if \code{TRUE} (default), user is prompted with a double check menu
#' @param recursive if \code{TRUE}, all sub-directories of the supplied directories will also be searched
#' @param file_types a vector of the file types for which the function finds and deletes. Defaults are \code{c("jpg","jpeg","png","tif")}
#' @return \code{NULL}
sweep_images <- function(dir,
                         ask = TRUE,
                         recursive = FALSE,
                         file_types = c("jpg", "jpeg", "png", "tif")){
  file_types <- paste0(paste0(".",file_types),collapse = "|")
  files <- list.files(dir, pattern = file_types, full.names = TRUE, recursive = recursive)
  nfiles <- length(files)

  if(nfiles < 1){
    message(paste0("\n",add_quote(dir), " is empty!"))
  } else {
    choice <- 1
    if(ask){
      choice <- utils::menu(title = paste0("You are attempting to delete ",
                                    nfiles , " images in ", add_quote(dir),
                                    "\nContinue?"),
                     choices = c("Yes", "No")
      )
    }
  }

  if(choice == 1){
    for (i in seq_along(files)){
      file.remove(files[i])
    }
  }
  message("\nSweep completed!")
  return(invisible(NULL))
}

