
#' @export
phi.min.obs.test<-function(mean.phi.T,n,min.phi=0.005,max.phi=1,a=14/9){
  lambda<-mean.phi.T*(2-a)/(1-a)*((max.phi^(1-a)-min.phi^(1-a))/(max.phi^(2-a)-min.phi^(2-a))) # lambda = mean phi.T / mean phi
  prob<-(1-dpois(1,lambda = lambda))^n
  return(prob)
} #Test probability of not observing the minimum phi in n samples of leaves
#Prob of seeing phi min = P(k=1|lambda=mean.phi.T/mean.phi)

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
id.screener<-function(data,variable = c("percHerbPlant",
                                        "percLf",
                                        "leaf",
                                        "plant"),cond=TRUE,min.leaves = 1,min.plants = 15,
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
        dplyr::filter(!is.na(percHerbPlant2)) %>%
        dplyr::filter(eval(parse(text = cond))) %>%
        dplyr::filter(n.avg >= min.leaves) %>%
        dplyr::group_by(surveyID) %>%
        dplyr::mutate(n.plants.in.survey = length(unique(plantID3.long)),
                      n.avg.med = median(n.avg)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(n.plants.in.survey >=
                        min.plants) %>%
        dplyr::mutate(percHerbPlant2 = ifelse(percHerbPlant2 < min.herb * 1/n.avg.med,
                                              0,
                                              percHerbPlant2)) %>%
        dplyr::group_by(surveyID) %>%
        dplyr::mutate(all_zero_check = all(percHerbPlant2 == 0)) %>%
        dplyr::filter(!all_zero_check) %>%
        dplyr::select(surveyID) %>%
        unique() %>%
        unlist()
    }
  }

  return(unname(ids))
}



#' @export
get.data.lists<-function(data,variable,id.name,ok.ids,n.sim=NA,
                         min.phi.prop=FALSE,max.phi.prop=FALSE,min.phi=0.005,max.phi=1,
                         correct.min.phi=FALSE,try=FALSE,group.n="n.avg",sim.pred=TRUE,...){
  variable<-ifelse(variable=="l","percLf",
                   ifelse(variable=="p1","percHerbPlant",
                          ifelse(variable=="p2","percHerbPlant2",
                                 variable)))
  if(is.list(data)){
    if(variable=="percLf"){
      data<-data[["long"]]
    } else {
      data<-data[["short"]]
    }
  }
  data$ids<-data[,id.name,drop=T]
  data<-data %>% dplyr::filter(ids%in%ok.ids)

  n.ids<-length(unique(data$ids))
  herb.obs.list<-vector(mode="list",length=n.ids)
  herb.pred.list<-vector(mode="list",length=n.ids)
  id.vector<-unique(data$ids)

  for (i in seq_len(n.ids)){
    herb.obs.list[[i]]<-data[data$ids==id.vector[i],variable,drop=T]
    if(sim.pred){
      if(id.name=="surveyID"){
        leaf.prop.of.whole<-1/mean(data[data$ids==id.vector[i],
                                        group.n,drop=T],na.rm = T)
      }
      min.phi.corrected<-ifelse(correct.min.phi,
                                mean(data[data$ids==id.vector[i],
                                          "phi.min.corrected",drop=T],na.rm = T),
                                min.phi)
      herb.pred.list[[i]]<-tryCatch(
        allometry.herb.quasi.sim(
          mean.phi.T = mean(herb.obs.list[[i]]),
          n.sim = ifelse(is.na(n.sim),
                         length(herb.obs.list[[i]]),
                         n.sim),
          min.phi = ifelse(min.phi.prop&id.name=="surveyID",
                           min.phi.corrected*leaf.prop.of.whole,
                           min.phi.corrected),
          max.phi = ifelse(max.phi.prop&id.name=="surveyID",
                           max.phi*leaf.prop.of.whole,
                           max.phi),
          ...),
        error=function(e){
          message(e,"\n")
          message("error id = ",i,"\n")
          if(try){
            return(NA)
          } else {
            stop(id.vector[i])
          }
        }
      )
    } else {
      herb.pred.list[[i]]<-NA
    }

    cat(signif(i/n.ids,2)*100,"%     ",
        id.vector[i],"                                ",
        "\r")
  }
  names(herb.obs.list)<-id.vector
  names(herb.pred.list)<-id.vector

  out<-list(
    "herb.obs.list"=herb.obs.list,
    "herb.pred.list"=herb.pred.list
  )
  return(out)
}



.randomize_leaves_engin <- function(data.list, meta.data, summarise.plant){
  #first column of meta.data is the surveyID
  #second column of meta.data is the plant id
  #third column of the meta.data is the number of leaves on that plant
  data.list.length <- length(data.list)
  rand.list <- vector(mode = "list", length = data.list.length)
  survey.ids <- names(data.list)
  for (i in seq_len(data.list.length)){
    n.leaves<-length(data.list[[i]])
    data<-meta.data[meta.data[,1]==survey.ids[i],]
    rand.list[[i]]<-data.frame(
      "herb"=data.list[[i]],
      "group"=rep(data[,2,drop=TRUE],data[,3,drop=TRUE]) %>%
        sample(size = n.leaves,replace = FALSE)
    )
  }

  if(summarise.plant){
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

#' @export
randomize_leaves<-function(data.list.list, meta.data, summarise.plant = TRUE, bind.lists = NULL){
  out <- lapply(
    data.list.list,
    FUN = function(x){
      .randomize_leaves_engin(data.list = x,
                             meta.data = meta.data[,c("surveyID","plantID3.long","n.avg")],
                             summarise.plant = summarise.plant)
    }
  )
  names(out) <- paste0("rand.",names(data.list.list))


  if(!is.null(bind.lists)){
    for (i in seq_along(bind.lists)){
      bind.lists[[i]]<-bind.lists[[i]][names(out[[1]])]
    }
    out<-c(out,bind.lists)
  }
  return(out)
}
