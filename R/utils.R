#### Convenience functions ####

#' @export
slide.show<-function(data.list,expression,sleep.time=1,length=1){
  for (i in seq_len(length)){
    eval(parse(text=expression))
    Sys.sleep(sleep.time)
    cat(i,"/",length,      "\r")
  }
}



#' @export
combine.data.lists<-function(data.list,data.list2){
  for (i in seq_along(data.list2)){
    data.list2[[i]]<-data.list2[[i]][names(data.list[[1]])]
  }
  out<-c(data.list,data.list2)
  return(out)
}













