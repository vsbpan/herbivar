#### BSL Wrangling ####
#' @export
merge_BSL_lists<-function(BSL.lists3){
  #Need BSL.list3 to be list of list of list of BSL
  #Outer list is just for merging-- returns list of list of BSL
  #Middle list is a list of surveys or plants (i.e. unique data set)
  #Inner list is a list of chains of that data set
  #Middle list needs to be named!
  #Inner list can pass on chain id name
  keys<-unique(unlist(lapply(BSL.lists3, names)))
  temp.list<-vector(mode = "list",length = length(keys))
  names(temp.list)<-keys
  for (i in seq_along(keys)){
    temp.list[keys[i]]<-list(unlist(purrr::map(BSL.lists3,keys[i])))
  }
  return(temp.list)
}


#' @export
check.BSL.convergence<-function(BSL.fit.list,trace.plot=T,burn=100){
  burn<-burn+1
  out<-data.frame("theta"=NA,"rhat"=NA,
                  "ess_bulk"=NA,"ess_tail"=NA,
                  "n_chains"=NA,"iter"=NA)
  chains<-list()
  for (i in seq_len(ncol(BSL.fit.list[[1]]@theta))){
    chains[[i]]<-lapply(BSL.fit.list,function(x){
      x@theta[,i]
    }) %>% do.call("cbind",.)
    chains[[i]]<-array(chains[[i]],dim = c(nrow(chains[[i]]),ncol(chains[[i]])))
    out[i,1] <- paste0("theta",i)
    out[i,2] <- Rhat(chains[[i]][burn:nrow(chains[[i]]),])
    out[i,3] <- ess_bulk(chains[[i]][burn:nrow(chains[[i]]),])
    out[i,4] <- ess_tail(chains[[i]][burn:nrow(chains[[i]]),])
    out[i,5] <- ncol(chains[[i]])
    out[i,6] <- nrow(chains[[i]])
  }
  if(trace.plot){
    opar <- par("mfrow","mar")
    par(mfrow=c(length(chains),2),mar=c(4,4,0.5,0.5))
    for (i in seq_along(chains)){
      for(j in seq_len(ncol(chains[[i]])))
        if(j==1){
          plot(x = burn:nrow(chains[[i]]),
               y = chains[[i]][burn:nrow(chains[[i]]),j],
               type="l",
               col=j,
               lwd=1.5,
               xlab="iter",
               ylab=paste0("theta",i))
        } else {
          lines(x = burn:nrow(chains[[i]]),
                y = chains[[i]][burn:nrow(chains[[i]]),j],
                lwd=1.5,
                col=j)
        }
      as.vector(chains[[i]][burn:nrow(chains[[i]]),]) %>%
        density() %>%
        plot(lwd=2,
             main="")
    }
    par(opar)
  }
  out$burn<-burn-1
  return(out)
}

#' @export
BSL.summary<-function(BSL.fit.list,burn=100,sigfig=3,check=T,trace.plot=F,...){
  out<-data.frame("theta"=NA,"Estimate"=NA,"Std"=NA,"lower"=NA,"upper"=NA)
  for (i in seq_len(ncol(BSL.fit.list[[1]]@theta))){
    draws<-lapply(BSL.fit.list,function(x){
      x@theta[(burn+1):nrow(x@theta),i]
    }) %>%
      do.call("cbind",.) %>%
      as.vector()
    out[i,1] <- paste0("theta",i)
    out[i,2] <- mean(draws)
    out[i,3] <- sd(draws)
    out[i,4] <- quantile(draws,prob=0.025)
    out[i,5] <- quantile(draws,prob=0.975)
  }
  if(check){
    out2<-check.BSL.convergence(BSL.fit.list,
                                burn = burn,
                                trace.plot = trace.plot,
                                ...)
    out<-left_join(out,out2,by="theta")
  }
  out<-out %>% mutate_if(is.numeric,function(x){
    signif(x,digits = sigfig)
  })
  return(out)
}

#' @export
BSL.posterior.draws<-function(BSL.fit.list,n.draws=1000,burn=100){
  #BSL.fit.list is a list of chains
  draws.list<-list()
  for (i in seq_len(ncol(BSL.fit.list[[1]]@theta))){
    draws.list[[i]]<-lapply(BSL.fit.list,function(x){
      x@theta[(burn+1):nrow(x@theta),i]
    }) %>%
      do.call("cbind",.) %>%
      as.vector()
  }
  names(draws.list)<-paste0("theta",1:ncol(BSL.fit.list[[1]]@theta))
  draws<-do.call("cbind",draws.list)
  draws<-draws[sample(seq_len(nrow(draws)),size = n.draws,replace = F),]
  return(draws)
}

#' @export
predict_draws_allo_herb<-function(draws,n.sim=1000,summarise=T,
                                  draws.name=NULL,
                                  param=c("mean.phi.T"=0.1,
                                          "max.phi"=1,
                                          "min.phi"=0.005,
                                          "a"=14/9,
                                          "truncate"=T)){
  theta.names<-c("mean.phi.T","max.phi","min.phi","a","truncate")

  if(!is.null(draws.name)){
    colnames(draws)<-draws.name #Rename supplied draws matrix if names provided
  }

  if(any(!colnames(draws)%in%theta.names)){
    stop("Input parameter not supported","\n",
         "Acceptable parameters are: mean.phi.T, max.phi, min.phi, a, truncate")
  }
  draws.full<-matrix(ncol=length(theta.names),
                     nrow=nrow(draws),
                     dimnames=list(NULL, theta.names))

  for(i in seq_along(theta.names)){
    draws.full[,theta.names[i]]<-ifelse(
      theta.names[i]%in%colnames(draws),# Check if parameter already existed
      draws[,theta.names[i]], # Stays the same
      param[theta.names[i]] # Index default if none exist
    )
  }

  out<-mapply(FUN = function(mean.phi.T, max.phi, min.phi, a, truncate){
    allometry.herb.quasi.sim(mean.phi.T = mean.phi.T,
                             max.phi = max.phi,
                             min.phi = min.phi,
                             a = a,
                             truncate = truncate,
                             n.sim = n.sim)
  },
  draws.full[,"mean.phi.T"],
  draws.full[,"max.phi"],
  draws.full[,"min.phi"],
  draws.full[,"a"],
  draws.full[,"truncate"])

  if(summarise){
    out<-apply(out,2,probe_distribution)
  }

  out<-t(out)
  return(out) #Each row is a draw
}
