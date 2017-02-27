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
require(wrswoR)


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
simulate.random<-function(dist, fnum, weights=rep(1, length(dist)), B = 1) {
  out_df<-tibble(Nuclei=integer(), N=integer(), fusNum=integer())
  smax<-length(dist)
  init_dist <- dist
  out_dist <- rep_len(0, length(dist))
  # print(c(length(dist), length(weights), fnum))
  for (k in 1:B) { 
    dist <- init_dist
    # probs <- weights*dist
    for (i in 1:fnum) {
      c1<-sample(smax, 1, prob=weights*dist)
      dist[c1]<-dist[c1]-1
      c2<-sample(smax, 1, prob=weights*dist)
      dist[c2]<-dist[c2]-1
      dist[c1+c2]<-dist[c1+c2]+1
    }
    out_dist<-out_dist+dist
  }
  out_dist<-out_dist/B
  out_df<-tibble(Nuclei=1:smax, N=out_dist, fusNum=rep(fnum, smax)) %>% 
    filter(!(Nuclei>1 & N<1))
  # print(c(smax, max(out_df$Nuclei)))
  return(out_df)
}

simulate.random_c<-compiler::cmpfun(simulate.random)

simulate.random1<-function(dist, fnum, weights=identity, B = 1) {
  out_df<-tibble(Nuclei=integer(), N=integer(), fusNum=integer())
  init_dist <- dist
  out_dist <- rep_len(0, length(dist))
  # print(c(length(dist), length(weights), fnum))
  for (k in 1:B) { 
    dist <- init_dist
    # probs <- weights*dist
    for (i in 1:fnum) {
      smax<-length(dist)
      c1<-sample(smax, 1, prob=weights(dist))
      dist[c1]<-dist[c1]-1
      c2<-sample(smax, 1, prob=weights(dist))
      dist[c2]<-dist[c2]-1
      if (c1+c2 > smax) {
        dist<-c(dist, rep(0,smax))
      }
      dist[c1+c2]<-dist[c1+c2]+1
    }
    dl <- length(out_dist)-length(dist)
    if (dl>0) {
      out_dist<-out_dist+c(dist, rep(0,dl))
    } else {
      out_dist<-c(out_dist, rep(0,-dl))+dist
    }
  }
  out_dist<-out_dist/B
  out_df<-tibble(Nuclei=1:length(out_dist), N=out_dist, fusNum=rep(fnum, length(out_dist))) %>% 
    filter(!(Nuclei>1 & N<1))
  # print(c(smax, max(out_df$Nuclei)))
  return(out_df)
}

simulate.random1_c<-compiler::cmpfun(simulate.random1)



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
simulate.founder<-function(dist, fnum, weights=rep(1, length(dist))) {
  require(tibble)
  require(dplyr)
  smax<-length(dist)
  out_df<-tibble(Nuclei=integer(), N=integer(), fusNum=integer())
  for (i in 1:fnum) {
    c1<-sample(smax, 1, prob=dist*weights)
    dist[c1]<-dist[c1]-1
    dist[c1+1]<-dist[c1+1]+1
  }
  out_df<-tibble(Nuclei=1:smax, N=dist, fusNum=rep(fnum, smax)) %>% 
    filter(!(Nuclei>1 & N<1))
  return(out_df)
}

simulate.founder2<-function(dist, dist_all, fnum, weights=rep(1, length(dist))) {
  require(tibble)
  require(dplyr)
  smax<-length(dist)
  out_df<-tibble(Nuclei=integer(), N=integer(), fusNum=integer())
  for (i in 1:fnum) {
    c1<-sample(smax, 1, prob=dist*weights)
    dist[c1]<-dist[c1]-1
    if (c1>1) dist_all[c1]<-dist_all[c1]-1
    c2<-sample(smax, 1, prob=dist_all*weights)
    dist_all[c2]<-dist_all[c2]-1
    if (c2>1) dist[c2]<-dist[c2]-1
    dist[c1+c2]<-dist[c1+c2]+1
    dist_all[c1+c2]<-dist_all[c1+c2]+1
    
  }
  out_df<-tibble(Nuclei=1:smax, N=dist, fusNum=rep(fnum, smax)) %>% 
    filter(!(Nuclei>1 & N<1))
  return(out_df)
}



initial_dist <- function(data, fnum, mono_num) {
  counts <- data %>% group_by(Nuclei) %>% summarise(N = n())
  maxn <- sum(counts$Nuclei*counts$N)+fnum
  dist <- c(mono_num, rep(0, maxn))
  for (i in seq_len(nrow(counts))) {
    dist[counts$Nuclei[i]]<-counts$N[i]
  }
  return(dist)
}

initial_dist1 <- function(data, fnum, mono_num) {
  counts <- data %>% group_by(Nuclei) %>% summarise(N = n())
  maxn <- max(counts$Nuclei)*2
  dist <- c(mono_num, rep(0, maxn))
  for (i in seq_len(nrow(counts))) {
    dist[counts$Nuclei[i]]<-counts$N[i]
  }
  return(dist)
}




dist.to.cdf<-function(dist) {
  cumsum(dist)/sum(dist)
}


