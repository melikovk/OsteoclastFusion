# Tools to work with syncytia data

#' Function to read data in xls/xlsx format with each experimental condition
#' on a separate sheet, named accordingly and different field of vew in 
#' different columns. Each measurement is in separate cell. Different types
#' of measurements are next to each other.
#' Area | Nuclei | Area | Nuclei | Area | Nuclei | ...
#' a11  | n11    | a12  | n12    | a13  | n13
#' a21  | n21    | a22  | n22    | a23  | n23
#' ...  | ...    | ...  | ...    | ...  | ...
#' @param fname 
#'
#' @return tidy dataframe
#' @export
#'
#' @examples

require(tibble)
require(dplyr)
# require(wrswoR)
require(parallel)
require(doParallel)
require(foreach)
require(iterators)
no_of_cores = detectCores()

pad_zeros<-function(vec, lmax) {
  return(c(vec,rep(0, lmax-length(vec))))
}

read.excel<-function(fname) {
  require(readxl)
  require(dplyr)
  exp_names=excel_sheets(fname)
  wb<-lapply(exp_names, read_excel, path=fname)
  df <- bind_rows(mapply(reshape.sheet, wb, exp_names, SIMPLIFY=F))
  df$fname = fname
  return(df)
}

reshape.sheet <- function(data, exp_name) {
  require(dplyr)
  require(tibble)
  measures<-unique(colnames(data))
  mnum=length(measures)
  colnames(data)<-make.unique(colnames(data))
  fields_num<-ncol(data)/mnum
  select.field<-function(fnum) {
    data %>% select((1:mnum)+mnum*(fnum-1)) %>% na.omit() -> df 
    colnames(df)<-measures
    df<-add_column(df, field=rep(fnum,nrow(df)))
    return(df)
  }
  add.expname <- function(x) add_column(x, Exp=rep(exp_name,nrow(x)))
  bind_rows(lapply(1:fields_num, select.field)) %>% add.expname
}

#' Function to model different measures of multinucleation
#'
#' @param N 
#' @param nfus 
#' @param B 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
fusion.indexes.mod<-function(N, nfus=round(0.3*N), B=1) {
  require(dplyr)
  require(tibble)
  require(tidyr)
  nuc_num<-matrix(rep(0, B*nfus), nrow=B)
  syn_num<-matrix(rep(0, B*nfus), nrow=B)
  syn_size<-matrix(rep(0, B*nfus), nrow=B)
  for (i in 1:B) {
    dist_mat<-diag(rep(0,nfus+1))
    dist_mat[1,1]<-N
    for (j in 2:(nfus+1)) {
      dist_mat[j,]=dist_mat[j-1,]
      c1<-sample((nfus+1), 1, prob=dist_mat[j,])
      dist_mat[j,c1]<-dist_mat[j,c1]-1
      c2<-sample((nfus+1), 1, prob=dist_mat[j,])
      dist_mat[j,c2]<-dist_mat[j,c2]-1
      dist_mat[j,c1+c2]<-dist_mat[j,c1+c2]+1
    }
    nuc_num[i,]<-apply(dist_mat[-1,-1], 1, function(x) sum(x*(2:(length(x)+1))))
    syn_num[i,]<-apply(dist_mat[-1,-1], 1, sum)
    syn_size[i,]<-apply(dist_mat[-1,-1], 1, function(x) sum(x*(2:(length(x)+1)))/sum(x))
  }
  add.expname <- function(x, exp_name) add_column(x, Exp=rep(exp_name,nrow(x)))
  nuc_df<-as_tibble(nuc_num)  %>% (function(x) add.expname(x, "Nuclei count"))
  syn_df<-as_tibble(syn_num)  %>% (function(x) add.expname(x, "Syncitia count"))
  size_df<-as_tibble(syn_size)  %>% (function(x) add.expname(x, "Syncitia size"))
  bind_rows(nuc_df,syn_df,size_df) %>% gather(key=T, value=Val, -Exp) %>% 
    mutate(Step=as.integer(substr(T,2,100))) %>% select(-T)
}


#' Function to simulate random cell fusions
#'
#' @param dist 
#' @param fnum 
#' @param weights 
#'
#' @return
#' @export
#'
#' @examples
#' 

simulate.random<-function(dist, fnum, weights=identity, B = 1) {
  out_df<-tibble(Nuclei=integer(), N=integer(), fusNum=integer())
  if (!is.vector(dist)) {
    
  }
  out_dists <- lapply(1:B, function(x) dist)
  fnum <- c(0, sort(fnum))
  fnum_len <- length(fnum)
  for (tpoint in 2:fnum_len) {
    for (f in fnum[tpoint-1]:(fnum[tpoint]-1)) {
      # next fusion
      for (k in 1:B) {
        # Advance each replicate 1 step forward
        max_nuc<-length(out_dists[[k]])
        c1<-sample(max_nuc, 1, replace = TRUE, prob=weights(out_dists[[k]]))
        out_dists[[k]][c1]<-out_dists[[k]][c1]-1
        c2<-sample(max_nuc, 1, replace = TRUE, prob=weights(out_dists[[k]]))
        out_dists[[k]][c2]<-out_dists[[k]][c2]-1
        if (c1+c2 > max_nuc) {
          out_dists[[k]]<-c(out_dists[[k]], rep(0,max_nuc))
        }
        out_dists[[k]][c1+c2]<-out_dists[[k]][c1+c2]+1
      }
    }
    # average replicates and record to tibble
    max_dist<-max(sapply(out_dists, length))
    out_matrix<-matrix(unlist(lapply(out_dists, function(x) pad_zeros(x,max_dist))), 
                       nrow=max_dist,ncol=B)
    out_mean<-.rowMeans(out_matrix, max_dist, B)
    out_df<-bind_rows(out_df, tibble(Nuclei=1:max_dist, N=out_mean, 
                                     fusNum=rep(fnum[tpoint], max_dist)) %>% 
                        filter(!(Nuclei>1 & N<1))) 
  }
  return(out_df)
}

psimulate.random<-function(dist, fnum, weights=identity, B = 1) {
  # cl = makeCluster(no_of_cores, type = "SOCK")
  registerDoParallel(no_of_cores)
  # on.exit(stopCluster(cl))
  out_df<-tibble(Nuclei=integer(), N=integer(), fusNum=integer())
  if (!is.vector(dist)) {
    
  }
  out_dists <- lapply(1:B, function(x) dist)
  fnum <- c(0, sort(fnum))
  fnum_len <- length(fnum)
  for (tpoint in 2:fnum_len) {
    # for (k in 1:B) 
    out_dists<-foreach(next_dist=iter(out_dists)) %dopar% {
      # next simulation
      for (f in fnum[tpoint-1]:(fnum[tpoint]-1)) {
        # next fusion
        max_nuc<-length(next_dist)
        c1<-sample(max_nuc, 1, replace = TRUE, prob=weights(next_dist))
        next_dist[c1]<-next_dist[c1]-1
        c2<-sample(max_nuc, 1, replace = TRUE, prob=weights(next_dist))
        next_dist[c2]<-next_dist[c2]-1
        if (c1+c2 > max_nuc) {
          next_dist<-c(next_dist, rep(0,max_nuc))
        }
        next_dist[c1+c2]<-next_dist[c1+c2]+1
      }
      next_dist
    }
    # average replicates and record to tibble
    max_dist<-max(foreach(dist=out_dists, .combine=c, .multicombine = T) %dopar% length(dist))
    #out_mean<-(foreach(dist=out_dists, .combine='+', .export=c("pad_zeros")) 
    #         %dopar% pad_zeros(dist,max_dist))/B
    out_matrix<-matrix(foreach(dist=out_dists, .combine=c, .multicombine=T, .export=c("pad_zeros")) %dopar% 
     pad_zeros(dist,max_dist),max_dist, B)
    # out_matrix<-matrix(unlist(lapply(out_dists, function(x) pad_zeros(x,max_dist))), 
    #                   nrow=max_dist,ncol=B)
    out_mean<-.rowMeans(out_matrix, max_dist, B)
    out_df<-bind_rows(out_df, tibble(Nuclei=1:max_dist, N=out_mean, 
                                     fusNum=rep(fnum[tpoint], max_dist)) %>% 
                        filter(!(Nuclei>1 & N<1))) 
  }
  return(out_df)
}

psimulate.random_c<-compiler::cmpfun(psimulate.random)
simulate.random_c<-compiler::cmpfun(simulate.random)

#' Function to simulate "founder" fusion model
#'
#' @param dist 
#' @param fnum 
#' @param weights 
#'
#' @return
#' @export
#'
#' @examples
#'
#'

simulate.founder<-function(dist, fnum, weights=identity, B=1) {
  out_df<-tibble(Nuclei=integer(), N=integer(), fusNum=integer())
  if (!is.vector(dist)) {
    
  }
  out_dists <- lapply(1:B, function(x) dist)
  fnum <- c(0, sort(fnum))
  fnum_len <- length(fnum)
  for (tpoint in 2:fnum_len) {
    for (f in fnum[tpoint-1]:(fnum[tpoint]-1)) {
      # next fusion
      for (k in 1:B) {
        # Advance each replicate 1 step forward
        max_nuc<-length(out_dists[[k]])
        c1<-sample(max_nuc, 1, replace = TRUE, prob=weights(out_dists[[k]]))
        out_dists[[k]][c1]<-out_dists[[k]][c1]-1
        if (c1+1 > max_nuc) {
          out_dists[[k]]<-c(out_dists[[k]], rep(0,max_nuc))
        }
        out_dists[[k]][c1+1]<-out_dists[[k]][c1+1]+1
      }
    }
    # average replicates and record to tibble
    max_dist<-max(sapply(out_dists, length))
    out_matrix<-matrix(unlist(lapply(out_dists, function(x) pad_zeros(x,max_dist))), 
                       nrow=max_dist,ncol=B)
    out_mean<-.rowMeans(out_matrix, max_dist, B)
    out_df<-bind_rows(out_df, tibble(Nuclei=1:max_dist, N=out_mean, 
                                     fusNum=rep(fnum[tpoint], max_dist)) %>% 
                        filter(!(Nuclei>1 & N<1))) 
  }
  return(out_df)
}

psimulate.founder<-function(dist, fnum, weights=identity, B = 1) {
  # cl = makeCluster(no_of_cores, type = "SOCK")
  registerDoParallel(no_of_cores)
  # on.exit(stopCluster(cl))
  out_df<-tibble(Nuclei=integer(), N=integer(), fusNum=integer())
  if (!is.vector(dist)) {
    
  }
  out_dists <- lapply(1:B, function(x) dist)
  fnum <- c(0, sort(fnum))
  fnum_len <- length(fnum)
  for (tpoint in 2:fnum_len) {
    # for (k in 1:B) 
    out_dists<-foreach(next_dist=iter(out_dists)) %dopar% {
      # next simulation
      for (f in fnum[tpoint-1]:(fnum[tpoint]-1)) {
        # next fusion
        max_nuc<-length(next_dist)
        c1<-sample(max_nuc, 1, replace = TRUE, prob=weights(next_dist))
        next_dist[c1]<-next_dist[c1]-1
        if (c1+1 > max_nuc) {
          next_dist<-c(next_dist, rep(0,max_nuc))
        }
        next_dist[c1+1]<-next_dist[c1+1]+1
      }
      next_dist
    }
    # average replicates and record to tibble
    max_dist<-max(foreach(dist=out_dists, .combine=c, .multicombine = T) %do% length(dist))
    #out_mean<-(foreach(dist=out_dists, .combine='+', .export=c("pad_zeros")) 
    #         %dopar% pad_zeros(dist,max_dist))/B
    out_matrix<-matrix(foreach(dist=out_dists, .combine=c, .multicombine=T, .export=c("pad_zeros")) %dopar% 
                         pad_zeros(dist,max_dist),max_dist, B)
    # out_matrix<-matrix(unlist(lapply(out_dists, function(x) pad_zeros(x,max_dist))), 
    #                   nrow=max_dist,ncol=B)
    out_mean<-.rowMeans(out_matrix, max_dist, B)
    out_df<-bind_rows(out_df, tibble(Nuclei=1:max_dist, N=out_mean, 
                                     fusNum=rep(fnum[tpoint], max_dist)) %>% 
                        filter(!(Nuclei>1 & N<1))) 
  }
  return(out_df)
}

simulate.founder_c<-compiler::cmpfun(simulate.founder)


simulate.founder2<-function(dist, mono_num, fnum, weights=identity, weights_all=identity, B=1) {
  out_df<-tibble(Nuclei=integer(), N=integer(), fusNum=integer())
  if (!is.vector(dist)) {
    
  }
  out_dists <- lapply(1:B, function(x) c(mono_num, dist))
  fnum <- c(0, sort(fnum))
  fnum_len <- length(fnum)
  for (tpoint in 2:fnum_len) {
    for (f in fnum[tpoint-1]:(fnum[tpoint]-1)) {
      # next fusion
      for (k in 1:B) {
        # Advance each replicate 1 step forward
        max_nuc<-length(out_dists[[k]])-1
        # choose first cell
        c1<-sample(max_nuc, 1, replace = TRUE, prob=weights(out_dists[[k]][-1]))
        out_dists[[k]][c1+1]<-out_dists[[k]][c1+1]-1
        # choose 2nd cell
        c2<-sample(max_nuc, 1, replace = TRUE, prob=weights_all(out_dists[[k]][-2]))
        if (c2>1) 
          out_dists[[k]][c2+1]<-out_dists[[k]][c2+1]-1
        else
          out_dists[[k]][c2]<-out_dists[[k]][c2]-1
        # add fusion product
        if (c1+c2 > max_nuc) {
          out_dists[[k]]<-c(out_dists[[k]], rep(0,max_nuc))
        }
        out_dists[[k]][c1+c2+1]<-out_dists[[k]][c1+c2+1]+1
      }
    }
    # average replicates and record to tibble
    max_dist<-max(sapply(out_dists, length))-1
    out_matrix<-matrix(unlist(lapply(out_dists, function(x) pad_zeros(x[-1],max_dist))), 
                       nrow=max_dist,ncol=B)
    out_mean<-.rowMeans(out_matrix, max_dist, B)
    out_df<-bind_rows(out_df, tibble(Nuclei=1:max_dist, N=out_mean, 
                                     fusNum=rep(fnum[tpoint], max_dist)) %>% 
                        filter(!(Nuclei>1 & N<1))) 
  }
  return(out_df)
}

psimulate.founder2<-function(dist, mono_num, fnum, weights=identity, weights_all=identity, B=1, ret.matrix=FALSE) {
  registerDoParallel(no_of_cores)
  #out_df<-tibble(Nuclei=integer(), N=integer(), fusNum=integer())
  if (!is.vector(dist)) {
    
  }
  max_nuc <- length(dist)
  out_dists <- lapply(1:B, function(x) c(mono_num, dist, rep.int(0, max_nuc)))
  fnum <- c(0, sort(fnum))
  fnum_len <- length(fnum)
  out_list<- vector('list', fnum_len-1)
  names(out_list)<-as.character(fnum[-1])
  for (tpoint in 2:fnum_len) {
    out_dists<-foreach(next_dist=iter(out_dists)) %dopar% {
      # next fusion
      for (f in fnum[tpoint-1]:(fnum[tpoint]-1)) {
        # Advance each replicate 1 step forward
        max_nuc<-length(next_dist)-1
        # choose first cell
        c1<-sample(max_nuc, 1, replace = TRUE, prob=weights(next_dist[-1]))
        next_dist[c1+1]<-next_dist[c1+1]-1
        # choose 2nd cell
        c2<-sample(max_nuc, 1, replace = TRUE, prob=weights_all(next_dist[-2]))
        if (c2 > 1) 
          next_dist[c2+1]<-next_dist[c2+1]-1
        else 
          next_dist[c2]<-next_dist[c2]-1
        # add fusion product
        if (c1+c2 > max_nuc) {
          next_dist<-c(next_dist, rep(0,max_nuc))
        }
        next_dist[c1+c2+1]<-next_dist[c1+c2+1]+1
      }
      next_dist
    }
    # average replicates and record to tibble
    max_dist<-max(sapply(out_dists, length))-1
    out_matrix<-matrix(unlist(lapply(out_dists, function(x) pad_zeros(x[-1],max_dist))), 
                       nrow=max_dist,ncol=B)
    out_mean<-.rowMeans(out_matrix, max_dist, B)
    max_syn <- Position(function(x) x>0, out_mean, right=TRUE)
    # out_mean<-out_mean/sum(out_mean)
    out_list[[as.character(fnum[tpoint])]]<-out_mean[1:max_syn]
  }
  if (ret.matrix) {
    max_dist<-max(sapply(out_list, length))
    out_data <- matrix(unlist(lapply(out_list, function(x) pad_zeros(x,max_dist))), 
                       nrow=max_dist,ncol=fnum_len-1) 
  } else {
    out_data<-foreach(odist=iter(out_list), fus_num=fnum[-1], .combine = bind_rows, .multicombine = TRUE) %do% {
      tibble(Nuclei=seq_along(odist), N=odist, fusNum=rep.int(fus_num, length(odist))) %>% 
        filter(N>0 | Nuclei==1)
    }
  }
  return(out_data)
}

simulate.founder2_c<-compiler::cmpfun(simulate.founder2)
psimulate.founder2_c<-compiler::cmpfun(psimulate.founder2)

fit.founder2<-function(data, mono, control.name, max.step=100) {
  control<- data %>% filter(Exp==control.name) %>% summarise(fuse.num = sum(Nuclei-1), syn.num=n())
  fnums <- data %>% filter(Exp!=control.name) %>% group_by(Exp) %>% summarise(fuse.num = sum(Nuclei-1), syn.num=n())
  fnums$fuse.num<-fnums$fuse.num-control$fuse.num
  fnums$syn.num<-fnums$syn.num-control$syn.num
  # create data matrix
  exp_matrix<-foreach(name=fnums$Exp[order(fnums$fuse.num)]) %do% {
    initial_dist(data %>% filter(Exp==name),0)[-1]
  }
  max_dist<-max(sapply(exp_matrix, length))
  exp_matrix <- matrix(unlist(lapply(exp_matrix, function(x) pad_zeros(x,max_dist))), 
                     nrow=max_dist,ncol=nrow(fnums))
  # create initial dist
  dist0<-initial_dist(data %>% filter(Exp==control.name),0)
  minErr <- Inf
  bestPar <- 0
  for (fit.num in seq(from=max(fnums$syn.num), length.out = max.step)) {
      dist0[1]<-fit.num
      # cat(fit.num, "\n")
      model <- psimulate.founder2(dist0, mono, fnums$fuse.num, 
                                  weights = function(x) x*(1:length(x)), B=100, ret.matrix=T)
      if (nrow(model)-1 < max_dist) {
        model<-apply(model, 2, function(x) pad_zeros(x, max_dist+1))
      } 
      model<-model[2:(max_dist+1),]  
      error<-sum((model-exp_matrix)^2)
      if (error<minErr) {
        minErr <-error
        bestPar <- fit.num
        cat("Current best fit is ", bestPar, " with error of ", minErr, "\n")
      }
      if (fit.num %% 500 == 0) {
        cat(fit.num, " Current best fit is ", bestPar, " with error of ", minErr, "\n")
      }
  }
}

initial_dist <- function(data, mono_num) {
  counts <- data %>% group_by(Nuclei) %>% summarise(N = n())
  maxn <- max(counts$Nuclei)
  dist <- c(round(mono_num), rep(0, maxn-1))
  for (i in seq_len(nrow(counts))) {
    dist[counts$Nuclei[i]]<-counts$N[i]
  }
  return(dist)
}

dist.to.vec<-function(data) {
  maxn <- max(data$Nuclei)
  dist <- rep(0, maxn)
  for (i in seq_along(data$Nuclei)) {
    dist[data$Nuclei[i]]<-data$N[i]
  }
  return(dist)
}

subtruct.dist<-function(data, dist) {
  for (i in seq_along(dist)) {
    data$N[data$Nuclei==i]<-data$N[data$Nuclei==i]-dist[i]
  }
  data$N[data$N<0]<-0
  return(data)
}


dist.to.cdf<-function(dist, norm=TRUE) {
  ret <- cumsum(dist)
  if (norm) ret<-ret/sum(dist)
  return(ret)
}


