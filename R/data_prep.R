

#' @export
rep_data.frame<-function(x,n){
  if(inherits(x,"data.frame")){
    if(n>1){
      x.og<-x
      for (i in seq_len(n-1)){
        x<-rbind(x,x.og)
      }
    }
  }
  return(x)
}


#' @export
.id_screener<-function(data,variable = c("percHerbPlant",
                                        "percLf",
                                        "leaf",
                                        "plant"),
                      cond=TRUE,
                      min.leaves = 1,
                      min.plants = 15,
                      min.herb = 0.005){
  variable <- match.arg(variable)
  variable <- switch(variable,
                     "leaf" = "leaf",
                     "plant" = "plant",
                     "percLf" = "leaf",
                     "percHerbPlant" = "plant")

  if(is.list(data)){
    if(variable=="leaf"){
      data<-data[["long"]]
    } else {
      data<-data[["short"]]
    }
  }

  if(variable=="leaf"){
    ids<-data %>%
      dplyr::filter(eval(parse(text = cond))) %>%
      dplyr::filter(!is.na(percLf)) %>%
      dplyr::group_by(plantID3.long) %>%
      dplyr::mutate(n.leaves.on.plant= length(unique(leafID))) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(surveyID) %>%
      dplyr::mutate(n.plants.in.survey = length(unique(plantID3.long))) %>%
      dplyr::ungroup() %>%
      dplyr::filter(n.leaves.on.plant>=min.leaves) %>%
      dplyr::filter(n.plants.in.survey>=min.plants) %>%
      dplyr::mutate(percLf = ifelse(percLf < min.herb,
                    0,
                    percLf)) %>%
      dplyr::group_by(plantID3.long) %>%
      dplyr::mutate(all_zero_check = all(percLf < min.herb)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!all_zero_check) %>%
      dplyr::select(plantID3.long) %>%
      unique() %>%
      unlist()
  } else {
    if(variable=="plant"){
      ids<-data %>%
        dplyr::filter(!is.na(percHerbPlant)) %>%
        dplyr::filter(eval(parse(text = cond))) %>%
        dplyr::filter(n.avg >= min.leaves) %>%
        dplyr::group_by(surveyID) %>%
        dplyr::mutate(n.plants.in.survey = length(unique(plantID3.long)),
                      n.avg.med = median(n.avg)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(n.plants.in.survey >=
                        min.plants) %>%
        dplyr::mutate(percHerbPlant = ifelse(percHerbPlant < min.herb * 1/n.avg.med,
                                              0,
                                              percHerbPlant)) %>%
        dplyr::group_by(surveyID) %>%
        dplyr::mutate(all_zero_check = all(percHerbPlant == 0)) %>%
        dplyr::filter(!all_zero_check) %>%
        dplyr::select(surveyID) %>%
        unique() %>%
        unlist()
    }
  }

  return(unname(ids))
}



#' @export
.get_data_list<-function(data,
                         variable = c("percHerbPlant",
                                      "percLf",
                                      "leaf",
                                      "plant"),
                         id.name,
                         ok.ids,
                         group.n = "n.avg"){
  variable <- match.arg(variable)
  variable <- switch(variable,
                     "leaf" = "percLf",
                     "plant" = "percHerbPlant",
                     "percLf" = "percLf",
                     "percHerbPlant" = "percHerbPlant")
  if(is.list(data)){
    if(variable == "percLf"){
      data<-data[["long"]]
    } else {
      data<-data[["short"]]
    }
  }

  data$ids <- data[,id.name,drop=TRUE]
  data <- data %>% dplyr::filter(ids%in%ok.ids)
  id.vector <- unique(ok.ids)

  herb.obs.list <- lapply(id.vector, function(x, data, variable){
    data[data$ids == x, variable, drop = TRUE]
  }, data = data, variable = variable)

  names(herb.obs.list) <- id.vector

  return(herb.obs.list)
}



.randomize_leaves_engin <- function(data_list, meta_data, summarise_plant){
  #first column of meta_data is the surveyID
  #second column of meta_data is the plant id
  #third column of the meta_data is the number of leaves on that plant

  survey.ids <- names(data_list)
  rand.list <- lapply(seq_along(data_list),
                      function(i){
                        data<-meta_data[meta_data[,1]==survey.ids[i],]
                        out <- data.frame(
                          "herb" = data_list[[i]],
                          "group"= sample(rep(data[,2,drop=TRUE],data[,3,drop=TRUE]),
                                          size = length(data_list[[i]]),
                                          replace = FALSE)
                        )
                        return(out)
                        })

  if(summarise_plant){
    rand.list <- lapply(rand.list,
                        function(x){
                          unname(tapply(x$herb,
                                        (x$group),
                                        function(z) {
                                          mean(z, na.rm = TRUE)
                                        })
                          )
                        })

    names(rand.list) <- survey.ids
  } else {
    rand.list <- lapply(rand.list,
                                 function(x){
                                   plant.ids <- unique(x$group)
                                   out <- lapply(plant.ids, function(group, data){
                                     data[data$group == group, "herb"]
                                   }, data = x)
                                   names(out) <- plant.ids
                                   return(out)
                                 }) %>%
      unlist(recursive = FALSE)
  }
  return(rand.list)
}

#' @title Randomize Leaf Herbivory Observations
#' @description Randomize leaf herbivory observations using provided meta data. Used to generate null distribution where leaves on individual plants are random draws from all leaves within a population or survey.
#' @param data_list a named list of data to be randomized. The name should match the data set ID provided with the meta_data argument.
#' @param meta_data a data.frame of meta data to be used in the randomization procedure. First column of meta_data is the data set ID. Second column of meta_data is the plant ID. Third column of the meta_data is the number of leaves on that plant that has herbivory data.
#' @param nboot the number of times to repeat the randomization procedure. Default is one.
#' @param summarise_plant if \code{TRUE} (default), the randomized leaf herbivory will be averaged by plant.
#' @return a list of list of numeric vectors
#' @export
randomize_leaves<-function(data_list, meta_data, nboot = 1, summarise_plant = TRUE){
  out <- lapply(
    seq_len(nboot),
    FUN = function(i, data_list, meta_data, summarise_plant){
      rout <- .randomize_leaves_engin(data_list,
                             meta_data = meta_data,
                             summarise_plant = summarise_plant)
      cat(i, "/", nboot, "                    \r")
      return(rout)
    },
    data_list = data_list,
    summarise_plant = summarise_plant,
    meta_data = meta_data
  )
  names(out) <- paste0("randomized",seq_len(nboot))

  return(out)
}
