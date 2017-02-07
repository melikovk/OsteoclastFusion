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


# Read Data
homedir="~/Documents/Osteoclast"
#homedir="C:/Users/melikovk/Documents/MacrophagePaper/CDF"
# Mice csv
mice<-read.csv(paste(homedir,"11 wk old mice CDF.csv",sep='/'))
levels(mice$Exp2)<-list("w.t."="Wt", "AnxA1 -/-"="AnxA1", "AnxA5 -/-"="AnxA5")
# Syn Antibodies csv
syn_ab<-read.csv(paste(homedir,"anti Syn1 on human monocyte CDF.csv",sep='/'))
levels(syn_ab$Exp)<-list("+LPC"="LPC on", 
                         "+LPC/-LPC"="LPC wash", 
                         "+LPC/-LPC Syn Ab"="anti Syn1 Santa cruz H280", 
                         "+LPC/-LPC IgG"="anti IgG")
# Syn peptide csv
syn_pept<-read.csv(paste(homedir,"Syn1 peptide on human monocyte CDF.csv",sep='/'))
levels(syn_pept$Exp)<-list("+LPC"="LPC on", 
                           "+LPC/-LPC"="LPC wash", 
                           "+LPC/-LPC Syn1 peptide"="Syn pep_20ug per ml", 
                           "+LPC/-LPC scr Syn1 peptide"="Syn pep sc_20ug per ml")
# Annexin peptide csv
anx_pept<-read.csv(paste(homedir,"Anx peptide on human monocytes CDF.csv",sep='/'))
levels(anx_pept$Exp)<-list("+LPC"="LPC on", 
                           "+LPC/-LPC"="LPC wash", 
                           "+LPC/-LPC Anx1 peptide"="Anx1 pep", 
                           "+LPC/-LPC scr Anx1 peptide"="Anx1 sc pep",
                           "+LPC/-LPC Anx5 peptide"="Anx5 pep", 
                           "+LPC/-LPC scr Anx5 peptide"="Anx5 sc pep")
# Annexin Ab csv
anx_ab<-read.csv(paste(homedir,"Anx antibody on human monocytes CDF.csv",sep='/'))
levels(anx_ab$Exp)<-list("+LPC"="LPC on", 
                         "+LPC/-LPC"="LPC wash", 
                         "+LPC/-LPC Anx1 Ab"="anti Anx 1_35ug_Wash", 
                         "+LPC/-LPC Anx5 Ab"="anti Anx 5_20ug_Wash", 
                         "+LPC/-LPC IgG"="anti IgG_20ug_Wash")
# a01 csv
a01<-read.csv(paste(homedir,"A01 on synchronized human monocytes CDF.csv",sep='/'))
levels(a01$Exp)<-list("+LPC"="LPC on", 
                      "+LPC/-LPC"="LPC wash", 
                      "+LPC/-LPC A01"="A01_120uM")
# Lact csv
lact<-read.csv(paste(homedir,"Lactadherin on synchronized human monocytes CDF.csv",sep='/'))
levels(lact$Exp)<-list("+LPC"="LPC on", 
                       "+LPC/-LPC"="LPC wash", 
                       "+LPC/-LPC 50 ul Lact"="Lact_50ul per ml",
                       "+LPC/-LPC 100 ul Lact"="Lact_100ul per ml")
# a01 raw csv
a01_raw<-read.csv(paste(homedir,"A01 on RAW cells.csv",sep='/'))
levels(a01_raw$Exp)<-list("+LPC"="LPC on", 
                          "+LPC/-LPC"="LPC wash", 
                          "+LPC/-LPC A01"="A01_120uM")
# Fusion kinetics csv
fusion_kinetics<-read.csv(paste(homedir,"Fusion kinetics RAW cells.csv",sep='/'))
levels(fusion_kinetics$Exp)<-list("No LPC"="No LPC",
                                  "+LPC"="LPC on", 
                                  "+LPC/-LPC 15 min"="15 min", 
                                  "+LPC/-LPC 30 min"="30 min",
                                  "+LPC/-LPC 60 min"="60 min",
                                  "+LPC/-LPC 90 min"="90 min")

# DC stamp csv
dcstamp<-read.csv(paste(homedir,"DC STAMP RAW cells inhibition.csv",sep='/'))
levels(dcstamp$Exp)<-list("+LPC"="LPC on",
                          "+LPC/-LPC"="LPC wash",
                          "+LPC/-LPC DC-STAMP Ab"="anti DC-STAMP ", 
                          "+LPC/-LPC IgG"="anti IgG")
# Anx Ab on RAW cell
anx_ab_raw<-read.csv(paste(homedir,"Anx antibody RAW cells CDF.csv",sep='/'))
levels(anx_ab_raw$Exp)<-list("+LPC"="LPC on", 
                             "+LPC/-LPC"="LPC wash", 
                             "+LPC/-LPC Anx1 Ab"="anti AnxA1", 
                             "+LPC/-LPC Anx5 Ab"="anti AnxA5", 
                             "+LPC/-LPC IgG"="anti IgG")
# Anx peptide on RAW cells
anx_pept_raw<-read.csv(paste(homedir,"Anx peptide RAW cells CDF.csv",sep='/'))
levels(anx_pept_raw$Exp)<-list("+LPC"="LPC on", 
                               "+LPC/-LPC"="LPC wash", 
                               "+LPC/-LPC Anx1 peptide"="Anx1 peptide", 
                               "+LPC/-LPC scr Anx1 peptide"="Anx1 sc peptide",
                               "+LPC/-LPC Anx5 peptide"="Anx5 peptide", 
                               "+LPC/-LPC scr Anx5 peptide"="Anx5 sc peptide")

# Fig 1 Data
fig1data_raw<-read.csv(paste(homedir,"Data.csv",sep='/'))
levels(fig1data_raw$Exp)<-list("+LPC"="LPC on", 
                               "+LPC/-LPC"="LPC wash", 
                               "Control"="no LPC")
# Create CDF data frames
# mice_cdf<-make.cdf(mice, "Nuclei", "Exp2")
# syn_ab_cdf<-make.cdf(syn_ab, "Nuclei", "Exp")
# syn_pept_cdf<-make.cdf(syn_pept, "Nuclei", "Exp")
# anx_pept_cdf<-make.cdf(anx_pept, "Nuclei", "Exp")
# anx_ab_cdf<-make.cdf(anx_ab, "Nuclei", "Exp")
# a01_cdf<-make.cdf(a01, "Nuclei", "Exp")
# lact_cdf<-make.cdf(lact, "Nuclei", "Exp")
# a01_raw_cdf<-make.cdf(a01_raw, "Nuclei", "Exp")
# fusion_kinetics_cdf<-make.cdf(fusion_kinetics, "Nuclei", "Exp")
# dcstamp_cdf<-make.cdf(dcstamp, "Nuclei", "Exp")
# anx_ab_raw_cdf<-make.cdf(anx_ab_raw,"Nuclei","Exp")
# anx_pept_raw_cdf<-make.cdf(anx_pept_raw,"Nuclei","Exp")
# fig1data_raw_cdf<-make.cdf(fig1data_raw,"Nuclei","Exp")

#Plot graphs
fsize=36
lsize=3
myPalette1<-c("firebrick4","chartreuse4", "yellow4","navy","maroon4")
myPalette2<-c("black","red", "green","yellow","blue","magenta")
mytheme<-theme_bw(base_size=fsize)+theme(panel.border=element_rect(color="black", size=1.5), panel.grid=element_blank(),
                                         legend.key.size=unit(5,"char"), legend.key.height=unit(2,"lines"),
                                         axis.text=element_text(size=rel(0.6)), legend.text=element_text(size=rel(0.6)),
                                         legend.title=element_blank(), legend.position=c(.6,.5))
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# ggplot(mice_cdf, aes(x=x, y=y, color=Exp))+geom_step(size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,50))
ggplot(mice, aes(x=Nuclei, color=Exp2))+stat_ecdf(geom='step', size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,50))
ggsave("mice_cdf.png", path=homedir, width=9, height=6.5, dpi=100)
# ggplot(syn_ab_cdf, aes(x=x, y=y, color=Exp))+geom_step(size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,80))
ggplot(syn_ab, aes(x=Nuclei, color=Exp))+stat_ecdf(geom='step', size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,60))
ggsave("syn_ab_cdf.png", path=homedir, width=9, height=6.5, dpi=100)
# ggplot(syn_pept_cdf, aes(x=x, y=y, color=Exp))+geom_step(size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,80))
ggplot(syn_pept, aes(x=Nuclei, color=Exp))+stat_ecdf(geom='step', size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,80))
ggsave("syn_pept_cdf.png", path=homedir, width=9, height=6.5, dpi=100)
# ggplot(anx_ab_cdf, aes(x=x, y=y, color=Exp))+geom_step(size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,80))
ggplot(anx_ab, aes(x=Nuclei, color=Exp))+stat_ecdf(geom='step', size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,80))
ggsave("anx_ab_cdf.png", path=homedir, width=9, height=6.5, dpi=100)
# ggplot(anx_pept_cdf, aes(x=x, y=y, color=Exp))+geom_step(size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,80))
ggplot(anx_pept, aes(x=Nuclei, color=Exp))+stat_ecdf(geom='step', size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,80))
ggsave("anx_pept_cdf.png", path=homedir, width=9, height=6.5, dpi=100)
#ggplot(a01_cdf, aes(x=x, y=y, color=Exp))+geom_step(size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,60))
ggplot(a01, aes(x=Nuclei, color=Exp))+stat_ecdf(geom='step', size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,60))
ggsave("a01_cdf.png", path=homedir, width=9, height=6.5, dpi=100)
# ggplot(lact_cdf, aes(x=x, y=y, color=Exp))+geom_step(size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,60))
ggplot(lact, aes(x=Nuclei, color=Exp))+stat_ecdf(geom='step', size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,60))
ggsave("lact_cdf.png", path=homedir, width=9, height=6.5, dpi=100)
# ggplot(dcstamp_cdf, aes(x=x, y=y, color=Exp))+geom_step(size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,200))
ggplot(dcstamp, aes(x=Nuclei, color=Exp))+stat_ecdf(geom='step', size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,200))
ggsave("dcstamp_cdf.png", path=homedir, width=9, height=6.5, dpi=100)
# ggplot(fusion_kinetics_cdf, aes(x=x, y=y, color=Exp))+geom_step(size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,150))
ggplot(fusion_kinetics, aes(x=Nuclei, color=Exp))+stat_ecdf(geom='step', size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,150))
ggsave("fusion_kinetics_cdf.png", path=homedir, width=9, height=6.5, dpi=100)
# ggplot(a01_raw_cdf, aes(x=x, y=y, color=Exp))+geom_step(size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,300))
ggplot(a01_raw, aes(x=Nuclei, color=Exp))+stat_ecdf(geom='step', size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,300))
ggsave("a01_raw_cdf.png", path=homedir, width=9, height=6.5, dpi=100)
# ggplot(anx_ab_raw_cdf, aes(x=x, y=y, color=Exp))+geom_step(size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,150))
ggplot(anx_ab_raw, aes(x=Nuclei, color=Exp))+stat_ecdf(geom='step', size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,150))
ggsave("anx_ab_raw_cdf.png", path=homedir, width=9, height=6.5, dpi=100)
# ggplot(anx_pept_raw_cdf, aes(x=x, y=y, color=Exp))+geom_step(size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,300))
ggplot(anx_pept_raw, aes(x=Nuclei, color=Exp))+stat_ecdf(geom='step', size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,300))
ggsave("anx_pept_raw_cdf.png", path=homedir, width=9, height=6.5, dpi=100)
# ggplot(fig1data_raw_cdf, aes(x=x, y=y, color=Exp))+geom_step(size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,300))
ggplot(fig1data_raw, aes(x=Nuclei, color=Exp))+stat_ecdf(geom='step', size=lsize)+mytheme+scale_colour_manual(values=cbPalette)+ylab("Cumulative frequency")+scale_x_continuous(name="Nuclei/Syncytium")+coord_cartesian(xlim=c(0,50))
ggsave("fig1data_raw_cdf.png", path=homedir, width=9, height=6.5, dpi=100)