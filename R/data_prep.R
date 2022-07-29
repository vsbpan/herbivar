
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
id.screener<-function(data,variable,cond=T,min.leaves=2,min.plants=2,
                      min.herb=0.00000001){
  variable<-ifelse(variable=="l","percLf",
                   ifelse(variable=="p1","percHerbPlant",
                          ifelse(variable=="p2","percHerbPlant2",
                                 variable)))
  if(!any(
    c("percLf","percHerbPlant","percHerbPlant2")%in%variable
  )){
    stop("Variable not supported")
  }
  if(is.list(data)){
    if(variable=="percLf"){
      data<-data[["long"]]
    } else {
      data<-data[["short"]]
    }
  }

  if(variable=="percLf"){
    ids<-data %>%
      dplyr::filter(eval(parse(text = cond))) %>%
      dplyr::filter(!is.na(percLf)) %>%
      group_by(plantID3.long) %>%
      mutate(n.leaves.on.plant= length(unique(leafID))) %>%
      ungroup() %>%
      group_by(surveyID) %>%
      mutate(n.plants.in.survey = length(unique(plantID3.long))) %>%
      ungroup() %>%
      dplyr::filter(n.leaves.on.plant>=min.leaves) %>%
      dplyr::filter(n.plants.in.survey>=min.plants) %>%
      group_by(plantID3.long) %>%
      mutate(avgHerb=mean(percLf)) %>%
      dplyr::filter(avgHerb>min.herb) %>%
      ungroup() %>%
      dplyr::select(plantID3.long) %>%
      unique() %>%
      unlist()
  } else {
    if(variable=="percHerbPlant"){
      data<-data %>%
        dplyr::filter(!is.na(percHerbPlant)) %>%
        mutate(percHerbPlant2=percHerbPlant)
    }
    ids<-data %>%
      dplyr::filter(eval(parse(text = cond))) %>%
      dplyr::filter(n.avg>=min.leaves) %>%
      group_by(surveyID) %>%
      mutate(n.plants.in.survey = length(unique(plantID3.long))) %>%
      ungroup() %>%
      dplyr::filter(n.plants.in.survey>=min.plants) %>%
      group_by(surveyID) %>%
      mutate(avgHerb=mean(percHerbPlant2)) %>%
      dplyr::filter(avgHerb>=min.herb) %>%
      ungroup() %>%
      dplyr::select(surveyID) %>%
      unique() %>%
      unlist()
  }
  ids<-ids %>% as.vector()
  return(ids)
}



#' @export
get.data.lists<-function(data,variable,id.name,ok.ids,n.sim=NA,
                         min.phi.prop=F,max.phi.prop=F,min.phi=0.005,max.phi=1,
                         correct.min.phi=T,try=F,group.n="n.avg",sim.pred=T,...){
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
          if(try==T){
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


#' @export
randomize.leaves<-function(data.list,n.sim=NA,try=F,
                           meta.data,bind.lists=NULL,...){
  herb.rand.obs.list<-vector(mode="list",length=length(data.list[[1]]))
  herb.rand.pred.list<-vector(mode="list",length=length(data.list[[1]]))

  for (i in seq_along(data.list[[1]])){
    n.leaves<-length(data.list[[1]][[i]])
    data<-meta.data %>% dplyr::filter(surveyID==names(data.list[[1]])[i])
    herb.rand.obs.list[[i]]<-data.frame(
      "herb"=data.list[[1]][[i]],
      "group"=rep(data$plantID3.long,data$n.avg) %>%
        sample(size = n.leaves,replace = F)
    ) %>%
      group_by(group) %>%
      summarise(percHerbPlant3=mean(herb)) %>%
      ungroup() %>%
      select(percHerbPlant3) %>%
      unlist() %>%
      unname()
  }
  names(herb.rand.obs.list)<-names(data.list[[1]])

  for (i in seq_along(data.list[[1]])){
    n.leaves<-length(data.list[[1]][[i]])
    data<-meta.data %>% dplyr::filter(surveyID==names(data.list[[1]])[i])
    herb.rand.pred.list[[i]]<-data.frame(
      "herb"=data.list[[2]][[i]],
      "group"=rep(data$plantID3.long,data$n.avg) %>%
        sample(size = n.leaves,replace = F)
    ) %>%
      group_by(group) %>%
      summarise(percHerbPlant3=mean(herb)) %>%
      ungroup() %>%
      select(percHerbPlant3) %>%
      unlist() %>%
      unname()
  }
  names(herb.rand.pred.list)<-names(data.list[[2]])

  out<-list(
    "herb.rand.obs.list"=herb.rand.obs.list,
    "herb.rand.pred.list"=herb.rand.pred.list
  )
  if(!is.null(bind.lists)){
    for (i in seq_along(bind.lists)){
      bind.lists[[i]]<-bind.lists[[i]][names(out[[1]])]
    }
    out<-c(out,bind.lists)
  }
  return(out)
}
