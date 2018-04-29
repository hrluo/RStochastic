setwd("/Users/xiaofeizhou/Documents/7620/Presentation")

# source("http://www.bioconductor.org/biocLite.R")
# biocLite("GEOquery")
library(Biobase)
library(GEOquery)

loadGeoFile <- function(geoFilename) {
  temp <- readLines(geoFilename)    # Load the file
  temp <- temp[grep("\t", temp)]    # Keep only lines with tabs
  temp <- gsub("\t$", "\tNA", temp) # Deal with NA
  temp <- strsplit(temp, "\t")      # Split the strings at each tab
  temp <- t(sapply(temp, unlist))   # Turn each line into a vector, transpose
  colnames(temp) <- temp[1, ]
  rownames(temp) <- temp[ ,1]
  #Remove the row/col names from the data, and return it.
  #Note that all the entries are strings/characters, not numeric!
  temp[-1,-1]
}

#temp <- getGEO(filename="/Users/xiaofeizhou/Downloads/GSE5418_series_matrix.txt.gz")

gse5418<-getGEO(filename='/Users/xiaofeizhou/Downloads/GSE5418_family.soft.gz',GSEMatrix = T)
Meta(gse5418)

gsmplatforms <- lapply(GSMList(gse5418),function(x) {Meta(x)$platform_id})
head(gsmplatforms)

gsmlist = Filter(function(gsm) {Meta(gsm)$platform_id=='GPL96'},GSMList(gse5418))
length(gsmlist)

Table(gsmlist[[71]])[1:5,]

n.ind<-22
ind<-vector("list",n.ind)
ind[[1]]<-cbind(Table(gsmlist$GSM123730),Table(gsmlist$GSM123732))
ind[[2]]<-cbind(Table(gsmlist$GSM123764),Table(gsmlist$GSM123765))
ind[[3]]<-cbind(Table(gsmlist$GSM123777),Table(gsmlist$GSM123778))
ind[[4]]<-cbind(Table(gsmlist$GSM123779),Table(gsmlist$GSM123780))
ind[[5]]<-cbind(Table(gsmlist$GSM123782),Table(gsmlist$GSM123784))
ind[[6]]<-cbind(Table(gsmlist$GSM123786),Table(gsmlist$GSM123787))
ind[[7]]<-cbind(Table(gsmlist$GSM123789),Table(gsmlist$GSM123791))
ind[[8]]<-cbind(Table(gsmlist$GSM123793),Table(gsmlist$GSM123795))
ind[[9]]<-cbind(Table(gsmlist$GSM123797),Table(gsmlist$GSM123799))
ind[[10]]<-cbind(Table(gsmlist$GSM123734),Table(gsmlist$GSM123736))
ind[[11]]<-cbind(Table(gsmlist$GSM123738),Table(gsmlist$GSM123740))
ind[[12]]<-cbind(Table(gsmlist$GSM123742),Table(gsmlist$GSM123744))
ind[[13]]<-cbind(Table(gsmlist$GSM123745),Table(gsmlist$GSM123746))
ind[[14]]<-cbind(Table(gsmlist$GSM123748),Table(gsmlist$GSM123750))
ind[[15]]<-cbind(Table(gsmlist$GSM123751),Table(gsmlist$GSM123752))
ind[[16]]<-cbind(Table(gsmlist$GSM123754),Table(gsmlist$GSM123756))
ind[[17]]<-cbind(Table(gsmlist$GSM123757),Table(gsmlist$GSM123758))
ind[[18]]<-cbind(Table(gsmlist$GSM123760),Table(gsmlist$GSM123761))
ind[[19]]<-cbind(Table(gsmlist$GSM123762),Table(gsmlist$GSM123763))
ind[[20]]<-cbind(Table(gsmlist$GSM123767),Table(gsmlist$GSM123769))
ind[[21]]<-cbind(Table(gsmlist$GSM123770),Table(gsmlist$GSM123771))
ind[[22]]<-cbind(Table(gsmlist$GSM123773),Table(gsmlist$GSM123774))




gse5418.ID.names<-as.character(Table(gsmlist[[71]])[,1])
log.ratio.mat<-matrix(0,nrow=length(gse5418.ID.names),ncol=n.ind)
avg.log.ratio<-numeric(length(gse5418.ID.names))
log.int.mat<-matrix(0,nrow=length(gse5418.ID.names),ncol=n.ind)
avg.log.int<-numeric(length(gse5418.ID.names))

min.abs<-2

for(i in 1:n.ind){
  min.abs<-min(c(abs(ind[[i]][,2]),abs(ind[[i]][,6])),min.abs)
  
  log.ratio.mat[,i]<-log(ind[[i]][,2]/ind[[i]][,6])
  log.int.mat[,i]<-0.5*(log(ind[[i]][,2])+log(ind[[i]][,6]))
}

min.abs

################################Get Genes Name################################
gpl96<-loadGeoFile('/Users/xiaofeizhou/Downloads/GPL96.annot')

gpl96.ID.names<-row.names(gpl96)

gse5418.gene.names<-gpl96[match(gse5418.ID.names,gpl96.ID.names),2]
gse5418.gene.names[1:5]
gse5418.gene.names<-unname(gse5418.gene.names,force = F)
gse5418.gene.names[1:5]
################################################################################

avg.log.ratio<-apply(log.ratio.mat,1,mean)
avg.log.int<-apply(log.int.mat,1,mean)



################################Find DE Genes ################################


library(DIME)

DIME.result<-DIME(avg.log.ratio,
                       avg=avg.log.int,
                       gng.K=4,inudge.K = 4,
                       gng.weights="lower",gng.max.iter = 3000,
                       inudge.weights = "lower", inudge.max.iter = 3000,
                       nudge.weights = "lower", nudge.max.iter = 3000)


data.out.sum<-data.frame(gse5418.gene.names,
                          DIME.result$best$fdr,
                          DIME.result$best$class,
                          avg.log.ratio,avg.log.int)
data.out.sig=data.out.sum[data.out.sum[,3]==1,]    #class needs to be 1
data.out.sig=data.out.sig[order(data.out.sig[,2]),]      #order it wrt fdr
colnames(data.out.sig)<-c("Gene","FDR","Class","Log Ratio","Log Intensity")
DIME.intermediate<-data.out.sig
DIME.final<-as.character(paste(data.out.sig[,1]))

save.image(file="Malaria22.RData")
################################################################################


#############################Obtain Enriched Gene Pathways#############################

# source("https://bioconductor.org/biocLite.R")
# biocLite("clusterProfiler")
# 
# source("https://bioconductor.org/biocLite.R")
# biocLite("ReactomePA")

setwd("/Users/xiaofeizhou/Documents/7620/Presentation")
load(file="Malaria22.RData")

library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(graphite)

DE.gene.names<-DIME.final[1:100]

malaria.new.name <- bitr(DE.gene.names, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
malaria.pathway <- enrichPathway(gene=malaria.new.name[,2],pvalueCutoff=0.05, readable=T)
summary(malaria.pathway)
barplot(malaria.pathway, showCategory=8,xlim=c(0,40))

#######################################################################################



#################################Get genes in pathways###################################
pathwayName<-summary(malaria.pathway)[,2]
gene.in.path<-vector("list",length(pathwayName))


for(i in 1:length(pathwayName)){
  p<-pathways("hsapiens", "reactome")[[pathwayName[i]]]
  p <- convertIdentifiers(p, "symbol")
  gene.in.path[[i]]<-unique(union(p@edges[,1],p@edges[,2]))
}
############################################################################################



#########################Construct Covariance from DE Genes#########################
sub.index<-match(DE.gene.names,gse5418.gene.names)
new.log.ratio.mat<-log.ratio.mat[sub.index,]

cov.mat<-cov(t(new.log.ratio.mat))
row.names(cov.mat)<-DE.gene.names
colnames(cov.mat)<-DE.gene.names
###########################################################################################



#############################Get Precision Matrix#############################
library(glasso)
library(Matrix)



nr<-100
rho<-seq(0.1,1,length=nr)
bic<-numeric(nr)
for(i in 1:nr){
  glasso.obj<-glasso(cov.mat,rho[i])
  pmat<-glasso.obj$wi
  diag(pmat)<-0
  p_off_d<-nnzero(pmat, na.counted = NA)/2
  bic[i]<--2*(glasso.obj$loglik)+p_off_d*log(n.ind)+4*p_off_d*log(nrow(pmat))
}

best<-which.min(bic)
best

pmat<-glasso(cov.mat,rho[best])$wi
diag(pmat)<-0
print(nnzero(pmat, na.counted = NA)/2)
########################################################################





################################Plot Networks with Pathways################################
source(file="ResultAnalysisFunctions.R")  #need to change

path.num<-3   #need to change
pathwayName<-pathwayName[1:path.num]


basePathNum=3    #need to change

colnames(pmat)<-DE.gene.names
row.names(pmat)<-DE.gene.names
new.pmat<-easy.pmat.fn(pmat,DE.gene.names,gene.in.path[[basePathNum]])

file.name.pt1<-pathwayName[basePathNum]
file.name.pt1<-gsub("/"," ",file.name.pt1)
file.name.pt1<-gsub(" ","",file.name.pt1)
png(filename = paste(file.name.pt1,".png",sep = ""),
    width = 960,height=960,res=150)
XFgraph.new.fn(new.pmat,gene.in.path,path.num,pmat,pathwayName,basePathNum)
title(main = "Malaria Networks")

dev.off()



########################################################################################




################################Plot Networks without Pathways################################
source(file="NetworkFnc_noPath.R")  #need to change

png(filename = paste("NoPathway",".png",sep=""),
    width = 960,height=960,res=150)
colnames(pmat)<-DE.gene.names
row.names(pmat)<-DE.gene.names

nrow.pmat<-nrow(pmat)
waste.index<-NULL
for(i in 1:nrow.pmat){
  if(nnzero(pmat[,i])==0){
    waste.index<-c(waste.index,i)
  }
}

new.pmat.index<-which(!(1:nrow.pmat) %in% waste.index)

new.pmat<-pmat[new.pmat.index,new.pmat.index]

XFgraph.noPath.fn(new.pmat)
title(main = "Malaria Networks")

dev.off()

########################################################################################