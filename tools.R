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
read.excel<-function(fname) {
  require(readxl)
  require(dplyr)
  exp_names=excel_sheets(fname)
  wb<-lapply(exp_names, read_excel, path=fname)
  bind_rows(mapply(reshape.sheet, wb, exp_names, SIMPLIFY=F))
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


simulate.random<-function(dist, fnum, weights=rep(1, length(dist))) {
  smax<-length(dist)
  for (i in 1:fnum) {
    c1<-sample(smax, 1, prob=weights*dist)
    dist[c1]<-dist[c1]-1
    c2<-sample(smax, 1, prob=weights*dist)
    dist[c2]<-dist[c2]-1
    dist[c1+c2]<-dist[c1+c2]+1
  }
  return(dist)
}

simulate.founder<-function(dist, fnum, weights=rep(1, length(dist))) {
  smax<-length(dist)
  for (i in 1:fnum) {
    c1<-sample(smax, 1, prob=dist*weights)
    dist[c1]<-dist[c1]-1
    dist[c1+1]<-dist[c1+1]+1
  }
  return(dist)
}


