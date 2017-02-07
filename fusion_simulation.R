require(plyr)
require(ggplot2)

fusion.sim<-function(dist, fnum, weights=rep(1, length(dist))) {
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

fusion.sim.kinetics<-function(dist, fnum, weights=rep(1, length(dist)), points=50) {
  smax<-length(dist)
  kin.dist<-matrix(0, nrow=length(dist), ncol=fnum+1)
  kin.dist[,1]<-dist
  for (i in 1:fnum) {
    c1<-sample(smax, 1, prob=dist*weights)
    dist[c1]<-dist[c1]-1
    c2<-sample(smax, 1, prob=dist*weights)
    dist[c2]<-dist[c2]-1
    dist[c1+c2]<-dist[c1+c2]+1
    kin.dist[,i+1]<-dist
  }
  tpoints<-round(seq(from=1, to=fnum+1, length.out=points))
  kin.dist<-kin.dist[,tpoints]
  dimnames(kin.dist)<-list(NULL, tpoints-1)
  return(kin.dist)
}

fusion.reservouir.sim.kinetics<-function(dist, fnum, weights=rep(1, length(dist)), points=50) {
  smax<-length(dist)
  kin.dist<-matrix(0, nrow=length(dist), ncol=fnum+1)
  kin.dist[,1]<-dist
  for (i in 1:fnum) {
    c1<-sample(smax, 1, prob=dist*weights)
    dist[c1]<-dist[c1]-1
    dist[c1+1]<-dist[c1+1]+1
    kin.dist[,i+1]<-dist
  }
  tpoints<-round(seq(from=1, to=fnum+1, length.out=points))
  kin.dist<-kin.dist[,tpoints]
  dimnames(kin.dist)<-list(NULL, tpoints-1)
  return(kin.dist)
}

fusion.reservouir.sim<-function(dist, fnum, weights=rep(1, length(dist))) {
  smax<-length(dist)
  for (i in 1:fnum) {
    c1<-sample(smax, 1, prob=dist*weights)
    dist[c1]<-dist[c1]-1
    dist[c1+1]<-dist[c1+1]+1
  }
  return(dist)
}

fusion.founder.sim<-function(dist_found, fnum, dist_foll, weights=NA) {
  smax<-length(dist_found)
  for (i in 1:fnum) {
    c1<-sample(smax, 1, prob=dist_found)
    dist_found[c1]<-dist_found[c1]-1
    dist_foll[c1]<-dist_foll[c1]-1
    c2<-sample(smax, 1, prob=dist_foll)
    dist_foll[c2]<-dist_foll[c2]-1
    if (c2!=1) 
      dist_found[c2]<-dist_found[c2]-1
    dist_found[c1+c2]<-dist_found[c1+c2]+1
    dist_foll[c1+c2]<-dist_foll[c1+c2]+1
  }
  return(dist_found)
}

make.cdf <- function(mydata, xname, yname){
  if (! is.element(xname, colnames(mydata)))
    stop("ColName does not exist")
  if (! is.element(yname, colnames(mydata)))
    stop("ColName does not exist")
  explevels<-levels(mydata[,yname])
  xmax<-max(mydata[,xname])
  finaldata<-read.csv(text="x, freq, y, Exp")
  for (expname in explevels) {
    tmpdata<-count(mydata[mydata[,yname]==expname,xname])
    tmpdata$y<-cumsum(tmpdata$freq)/sum(tmpdata$freq)
    tmpdata<-rbind(list(x=0, freq=0, y=0), tmpdata)
    tmpdata<-rbind(tmpdata, list(x=xmax, freq=0, y=1.00001))
    tmpdata$Exp<-expname
    finaldata<-rbind(finaldata, tmpdata)
  }
  return(finaldata)
}

get.dist<- function(cdfdata, expname, max.nuc) {
  mydist<-rep(0, max.nuc)
  mydata<-cdfdata[cdfdata$Exp==expname, c(1,2)]
  for (i in 2:nrow(mydata)-1) {
    mydist[mydata$x[i]]<-mydata$freq[i]
  }
  return(mydist)
}

plot.kinetics<-function(dist, lim=NULL) {
  rows<-nrow(dist)
  plot.data<-data.frame(x=numeric(0), y=numeric(0), Exp=character(0))
  for (i in 1:ncol(dist)) {
    tdata<-data.frame(x=2:rows, y=dist[-1,i])
    tdata<-tdata[tdata$y>0,]
    tdata$y<-cumsum(tdata$y)/sum(tdata$y)
    tdata<-rbind(list(x=0,y=0), tdata)
    tdata$Exp<-colnames(dist)[i]
    plot.data<-rbind(plot.data, tdata)
  }
  ggplot(plot.data, aes(x=x, y=y, color=Exp))+geom_step()+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+
    coord_cartesian(xlim=lim)
}

plot.dists<-function(exp.dist, sim.dist, lim=NULL) {
  exp.data<-data.frame(x=2:length(exp.dist), freq=exp.dist[-1])
  exp.data<-exp.data[exp.data$freq>0,]
  exp.data$y<-cumsum(exp.data$freq)/sum(exp.data$freq)
  exp.data<-rbind(list(x=0, freq=0, y=0), exp.data)
  exp.data$Exp<-"Experiment"
  sim.data<-data.frame(x=2:length(sim.dist), freq=sim.dist[-1])
  sim.data<-sim.data[sim.data$freq>0,]
  sim.data$y<-cumsum(sim.data$freq)/sum(sim.data$freq)
  sim.data<-rbind(list(x=0, freq=0, y=0), sim.data)
  sim.data$Exp<-"Simulation"
  plot.data<-rbind(exp.data,sim.data)
  ggplot(plot.data, aes(x=x, y=y, color=Exp))+geom_step()+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+
    coord_cartesian(xlim=lim)
}

# test data
#homedir="~/Documents/MacrophagePaper/CDF" #Linux
homedir="C:/Users/melikovk/Documents/MacrophagePaper/CDF"
syn_ab<-read.csv(paste(homedir,"anti Syn1 on human monocyte CDF.csv",sep='/'))
levels(syn_ab$Exp)<-list("+LPC"="LPC on", 
                         "+LPC/-LPC"="LPC wash", 
                         "+LPC/-LPC Syn Ab"="anti Syn1 Santa cruz H280", 
                         "+LPC/-LPC IgG"="anti IgG")
syn_ab_cdf<-make.cdf(syn_ab, "Nuclei", "Exp")
syn_ab_dist1<-get.dist(syn_ab_cdf, "+LPC", 500)
syn_ab_dist2<-get.dist(syn_ab_cdf, "+LPC/-LPC", 500)
fuse.num<-sum(syn_ab$Nuclei[syn_ab$Exp=='+LPC/-LPC'])-length(syn_ab$Nuclei[syn_ab$Exp=='+LPC/-LPC'])-
  sum(syn_ab$Nuclei[syn_ab$Exp=='+LPC'])+length(syn_ab$Nuclei[syn_ab$Exp=='+LPC'])
mono<-sum(syn_ab$Nuclei[syn_ab$Exp=='+LPC/-LPC'])*4+sum(syn_ab$Nuclei[syn_ab$Exp=='+LPC/-LPC'])-sum(syn_ab$Nuclei[syn_ab$Exp=='+LPC'])
syn_ab_dist1[1]<-mono
dist.sim<-fusion.sim(syn_ab_dist1, fuse.num)
plot.dists(syn_ab_dist2, dist.sim)
syn_ab_dist1b[1]<-sum(syn_ab$Nuclei[syn_ab$Exp=='+LPC/-LPC'])-sum(syn_ab$Nuclei[syn_ab$Exp=='+LPC'])

# Lact experiment siimulation
lact<-read.csv(paste(homedir,"Lactadherin on synchronized human monocytes CDF.csv",sep='/'))
levels(lact$Exp)<-list("+LPC"="LPC on", 
                       "+LPC/-LPC"="LPC wash", 
                       "+LPC/-LPC 50 ul Lact"="Lact_50ul per ml",
                       "+LPC/-LPC 100 ul Lact"="Lact_100ul per ml")
lact_cdf<-make.cdf(lact, "Nuclei", "Exp")
lact_dist1<-get.dist(lact_cdf, "+LPC", 500)
lact_dist2<-get.dist(lact_cdf, "+LPC/-LPC", 500)
lact_dist3<-get.dist(lact_cdf, "+LPC/-LPC 100 ul Lact", 500)
lact.fnum1<-sum(lact$Nuclei[lact$Exp=='+LPC/-LPC'])-length(lact$Nuclei[lact$Exp=='+LPC/-LPC'])-
  sum(lact$Nuclei[lact$Exp=='+LPC'])+length(lact$Nuclei[lact$Exp=='+LPC'])
lact.mono1<-sum(lact$Nuclei[lact$Exp=='+LPC/-LPC'])-sum(lact$Nuclei[lact$Exp=='+LPC'])
lact_dist1[1]<-round(lact.mono1*1.01)
dist.sim<-fusion.sim(lact_dist1, lact.fnum1)
plot.dists(lact_dist2, dist.sim)
lact.fnum2<-sum(lact$Nuclei[lact$Exp=='+LPC/-LPC 100 ul Lact'])-length(lact$Nuclei[lact$Exp=='+LPC/-LPC 100 ul Lact'])-
  sum(lact$Nuclei[lact$Exp=='+LPC'])+length(lact$Nuclei[lact$Exp=='+LPC'])
lact.mono2<-sum(lact$Nuclei[lact$Exp=='+LPC/-LPC 100 ul Lact'])-sum(lact$Nuclei[lact$Exp=='+LPC'])
lact_dist1[1]<-round(lact.mono2*1.15)
dist.sim<-fusion.sim(lact_dist1, lact.fnum2)
plot.dists(lact_dist3, dist.sim)

# DC stamp simulation
dcstamp<-read.csv(paste(homedir,"DC STAMP RAW cells inhibition.csv",sep='/'))
levels(dcstamp$Exp)<-list("+LPC"="LPC on",
                          "+LPC/-LPC"="LPC wash",
                          "+LPC/-LPC DC-STAMP Ab"="anti DC-STAMP ", 
                          "+LPC/-LPC IgG"="anti IgG")
dcstamp_cdf<-make.cdf(dcstamp, "Nuclei", "Exp")
dcstamp_dist1<-get.dist(dcstamp_cdf, "+LPC", 5000)
dcstamp_dist2<-get.dist(dcstamp_cdf, "+LPC/-LPC", 5000)
dcstamp_dist3<-get.dist(dcstamp_cdf, "+LPC/-LPC DC-STAMP Ab", 5000)
dcstamp.fnum1<-sum(dcstamp$Nuclei[dcstamp$Exp=='+LPC/-LPC'])-length(dcstamp$Nuclei[dcstamp$Exp=='+LPC/-LPC'])-
  sum(dcstamp$Nuclei[dcstamp$Exp=='+LPC'])+length(dcstamp$Nuclei[dcstamp$Exp=='+LPC'])
dcstamp.mono1<-sum(dcstamp$Nuclei[dcstamp$Exp=='+LPC/-LPC'])-sum(dcstamp$Nuclei[dcstamp$Exp=='+LPC'])
dcstamp_dist1[1]<-dcstamp.mono1*1.015
dist.sim<-fusion.sim(dcstamp_dist1, dcstamp.fnum1)
plot.dists(dcstamp_dist2, dist.sim)
dcstamp.fnum2<-sum(dcstamp$Nuclei[dcstamp$Exp=='+LPC/-LPC DC-STAMP Ab'])-length(dcstamp$Nuclei[dcstamp$Exp=='+LPC/-LPC DC-STAMP Ab'])-
  sum(dcstamp$Nuclei[dcstamp$Exp=='+LPC'])+length(dcstamp$Nuclei[dcstamp$Exp=='+LPC'])
dcstamp.mono2<-sum(dcstamp$Nuclei[dcstamp$Exp=='+LPC/-LPC DC-STAMP Ab'])-sum(dcstamp$Nuclei[dcstamp$Exp=='+LPC'])
dcstamp_dist1[1]<-dcstamp.mono2*1.15
dist.sim<-fusion.sim(dcstamp_dist1, dcstamp.fnum2)
plot.dists(dcstamp_dist3, dist.sim)

# Fusion kinetics simulation
fusion_kinetics<-read.csv(paste(homedir,"Fusion kinetics RAW cells.csv",sep='/'))
levels(fusion_kinetics$Exp)<-list("No LPC"="No LPC",
                                  "+LPC"="LPC on", 
                                  "+LPC/-LPC 15 min"="15 min", 
                                  "+LPC/-LPC 30 min"="30 min",
                                  "+LPC/-LPC 60 min"="60 min",
                                  "+LPC/-LPC 90 min"="90 min")
fusion_kinetics_cdf<-make.cdf(fusion_kinetics, "Nuclei", "Exp")
fusion_kinetics_dist1<-get.dist(fusion_kinetics_cdf, "+LPC", 5000)
fusion_kinetics_dist2<-get.dist(fusion_kinetics_cdf, "+LPC/-LPC 90 min", 5000)
fusion_kinetics_dist3<-get.dist(fusion_kinetics_cdf, "+LPC/-LPC 30 min", 5000)
fusion_kinetics_dist4<-get.dist(fusion_kinetics_cdf, "+LPC/-LPC 60 min", 5000)
fusion_kinetics.fnum1<-sum(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC/-LPC 90 min'])-length(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC/-LPC 90 min'])-
  sum(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC'])+length(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC'])
fusion_kinetics.mono1<-sum(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC/-LPC 90 min'])-sum(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC'])
fusion_kinetics_dist1[1]<-fusion_kinetics.mono1*1.03
dist.sim<-fusion.sim(fusion_kinetics_dist1, fusion_kinetics.fnum1)
plot.dists(fusion_kinetics_dist2, dist.sim, c(0,200))
fusion_kinetics.fnum2<-sum(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC/-LPC 30 min'])-length(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC/-LPC 30 min'])-
  sum(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC'])+length(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC'])
fusion_kinetics.mono2<-sum(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC/-LPC 30 min'])-sum(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC'])
fusion_kinetics_dist1[1]<-fusion_kinetics.mono2*1.2
dist.sim<-fusion.sim(fusion_kinetics_dist1, fusion_kinetics.fnum2)
plot.dists(fusion_kinetics_dist3, dist.sim)
fusion_kinetics.fnum3<-sum(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC/-LPC 60 min'])-length(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC/-LPC 60 min'])-
  sum(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC'])+length(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC'])
fusion_kinetics.mono3<-sum(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC/-LPC 60 min'])-sum(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC'])
fusion_kinetics_dist1[1]<-fusion_kinetics.mono3*1.1
dist.sim<-fusion.sim(fusion_kinetics_dist1, fusion_kinetics.fnum3)
plot.dists(fusion_kinetics_dist4, dist.sim)
fusion_kinetics.fnum4<-sum(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC/-LPC 90 min'])-length(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC/-LPC 90 min'])-
  sum(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC/-LPC 60 min'])+length(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC/-LPC 60 min'])
fusion_kinetics.mono4<-sum(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC/-LPC 90 min'])-sum(fusion_kinetics$Nuclei[fusion_kinetics$Exp=='+LPC/-LPC 60 min'])
fusion_kinetics_dist3[1]<-fusion_kinetics.mono4*1.05
dist.sim<-fusion.sim(fusion_kinetics_dist3, fusion_kinetics.fnum4)
plot.dists(fusion_kinetics_dist2, dist.sim, c(0,200))

# test

lact.fnum1<-sum(lact$Nuclei[lact$Exp=='+LPC/-LPC'])-length(lact$Nuclei[lact$Exp=='+LPC/-LPC'])-
  sum(lact$Nuclei[lact$Exp=='+LPC'])+length(lact$Nuclei[lact$Exp=='+LPC'])
lact.mono1<-sum(lact$Nuclei[lact$Exp=='+LPC/-LPC'])-sum(lact$Nuclei[lact$Exp=='+LPC'])
lact_dist1[1]<-0
lact_dist1[2]<-lact_dist1[2]+lact.mono1/2
dist.sim<-fusion.sim(lact_dist1, lact.fnum1-lact.mono1/2)
plot.dists(lact_dist2, dist.sim)
lact.syn1<-length(lact$Nuclei[lact$Exp=='+LPC/-LPC'])
lact.fnumt<-sum(lact$Nuclei[lact$Exp=='+LPC/-LPC'])-length(lact$Nuclei[lact$Exp=='+LPC/-LPC'])
dist.sim<-fusion.reservouir.sim(c(lact.syn1-sum(lact_dist1[-1]), lact_dist1[-1]), lact.fnum1, c(1:100, rep(100,400)))
plot.dists(lact_dist2, dist.sim)
dist.kin<-fusion.sim.kinetics(lact_dist1, lact.fnum1, points=10)
plot.kinetics(dist.kin)
dist.kin2<-fusion.reservouir.sim.kinetics(c(lact.syn1-sum(lact_dist1[-1]), lact_dist1[-1]), lact.fnum1, c(1:100, rep(100,400)), points=10)
plot.kinetics(dist.kin2)
