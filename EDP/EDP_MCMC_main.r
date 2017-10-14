##########
# Enriched Dirichlet Process in partitioning of high-dimensional data.
#(with general kernels K_Y_X,K_X replaceable, deafult to Gaussian/IG model)
# Author: Hengrui Luo (luo.619@osu.edu)
# Date: 2017-05-02/23
# Platform: R version 3.3.1 (2016-06-21)
# Length: 33,100 bytes, 881 lines.
##########
require(MASS)
require(Matrix)
require(gplots)
##########################################################################################
###Data input.
n<-24;
p<-20;
#We set up the number of observations and dimension of covariate X first.
#
y.data<-as.matrix(c(rnorm(n=12,14,1),rnorm(n=11,100,1),50) )
#One deviation, to see if EDP can detect it.
x.data<-rbind(
  mvrnorm(n=3,mu=rep(20,p),Sigma=diag(p) ),
  mvrnorm(n=9,mu=rep(30,p),Sigma=diag(p) ),
  mvrnorm(n=3,mu=rep(10,p),Sigma=diag(p) ),
  mvrnorm(n=3,mu=rep(20,p),Sigma=diag(p) ),
  mvrnorm(n=6,mu=rep(30,p),Sigma=diag(p) )
)
MHSteps<-1000;
#It should be 10000 or more in action.
totalSteps<-20000;
#You can increase it a bit if you have better computational power...it should be 10000 or more in action.
#The above data and sample size n, x dimension p are all you need to specify in order to use this code.
#Returning result contains a rho matrix showing the cluster partition structure proposed by [Wade et.al].
#[Wade et.al]Wade, Sara, et al. "Improving prediction from dirichlet process mixtures via enrichment." Journal of Machine Learning Research 15.1 (2014): 1041-1071.
y.cluster<-c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2);
#First level of DP is related with Y and theta.
x.cluster<-c(1,1,1,3,3,3,3,3,3,3,3,3,1,1,1,2,2,2,7,7,7,7,7,7);
#Second level of DP is related with X and psi; 
#The initial value is actually the correct cluster according to the generating procedure, we will see how EDP recovers this partition.

#since X has higher dimension than that of Y; a finer cluster structure is expected.
#In the observation, there are 7 clusters for covariates but only 4 clusters for responses.
#
##########################################################################################

##########################################################################################
require(ggplot2)
require(gridExtra)
matrix2list<-function(mat){
  nRow=dim(mat)[1];
  nCol=dim(mat)[2];
  if(nRow!=nCol){
    return(FALSE)
  }else{
    n=nRow;
    lis<-matrix(NA,nrow=n*n,ncol=3);
    colnames(lis)<-c("row","col","value");
    counter=1;
    for(i in 1:n){
      for(j in 1:n){
        lis[counter,]<-c(i,j,mat[i,j]);
        counter=counter+1;
      }
    }
    return(lis)
  }
}
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
analyze_plot<-function(rho,x.data,y.data){
  partitions<-paste(t(rho)[,1],t(rho)[,2])
  nums<-seq(1,dim(rho)[2],1)
  visdf<-as.data.frame(cbind(nums,x.data,y.data,partitions))
  require(ggplot2)
  x.data.astext<-round(x.data[,1],digits = 1)
  for(col_sel in 2:(dim(x.data)[2]) ){
    x.data.astext<-paste(x.data.astext,round(x.data[,col_sel],digits=1),sep=",")
  }
  x.data.astext<-paste("(",x.data.astext,")",sep="")
  for(i in 1:(dim(visdf)[2]-1)){
    visdf[,i]<-as.numeric(
      as.character(visdf[,i])
    )
  }
  g3<-ggplot(data=visdf)+
    geom_text(aes(x=nums,y=y.data),
              label=as.character(x.data.astext),size=5,alpha=.5, angle=45)+
    geom_point(aes(x=nums,y=y.data,colour=partitions),size=2)+
    xlab("# Observations")+
    scale_x_continuous()+
    theme_minimal()
  #require(multiplot)
  print(g3)
  
  plots <- list()
  for(xcol_sel in 1:p){
    g4<-ggplot(data=visdf)+
      geom_point(aes(x=as.numeric(as.character(visdf[,xcol_sel]))
                     ,y=y.data,colour=partitions),size=2)+
      ylab("Y") + xlab(paste("X_",xcol_sel))
    plots[[xcol_sel]] <- g4
    print(g4)
  }
  #multiplot(plotlist=plots, cols = 1, main = "Main title")
}
##########################################################################################

##########################################################################################
##partition storage and functionality.
analyze_partition<-function(rho){
  #Reframe the partition into non-gaping manner.
  result<-list();
  
  y.cluster<-rho[1,];
  x.cluster<-rho[2,];
  k<-unique(y.cluster);
  k_j<-matrix(NA,nrow=n,ncol=n);
  
  for(j in unique(y.cluster)){
    #jth_y.cluster<-cluster.y[cluster.y==k[j]]
    jth_cluster.x<-x.cluster[y.cluster==j]
    uni_cluster.x<-unique(jth_cluster.x)
    for(s in uni_cluster.x){
      #count how many element in each x cluster within jth cluster of y.
      k_j[j,s]<-sum(jth_cluster.x==s)
    }
  }
  #We work out the k and k_j's
  #The columns of k represents the unique levels of y clusters.
  #The rows of k_j represents the 
  #k_j[,1] represents how many unique x clusters within y clusters.
  #k_j[,2...] represents how many observations for each unique x.
  #
  result<-list(rho=rho,
               #The matrix representing the partition.
               k=k,
               #number of partitions in y variate.
               k_j=k_j);
  #number of partitions in x variate.
  
  return(result)
}
#
rho<-rbind(y.cluster,x.cluster);
result.rho<-analyze_partition(rho);
k_j<-result.rho$k_j;
#initialization of parameters.
theta<-matrix(1,nrow=dim(k_j)[1],ncol=1);
psi<-matrix(1,nrow=dim(k_j)[1],ncol=dim(result.rho$k_j)[2]);
#This is the partition matrix rho, each of its columns represents a partition.
#
#We need to define the hierarchical structure via partitioning the parameter space first.
#Here we can use two matrices, but in following sampling we need to put gamma hyperpriors on these precision alphas.
#These two matrices contains at most n columns since we cannot make more than n clusters out of n observations.
#
# Given such a partition matrix, we can compute its random partition according to [Wade et.al] Prop 1.
##########################################################################################

##########################################################################################
###Setup of kernels, replaceable.
K_Y_X<-function(yi,xi,theta){
  #I'd better add a dimension check for xi here...
  #if(length(xi)!=p){return(0);}
  beta.i=rep(theta,p);
  sigma.yi=1;
  return(dnorm(yi, mean=xi%*%beta.i,sd=sigma.yi) )
  #K(yi | xi, theta )=K(yi | xi, beta.i,sigma.y)
}
#only the prediction part requires MCMC since the joint density is not in closed form.
h_Y_X<-function(yi,xi,other_obs_X=matrix(),other_obs_Y=c()){
  require(cubature)
  if(is.null(dim(other_obs_X)) ){
    integrand <- function(x) {
      ret<-1;
      for(k in 1:dim(other_obs_X)[1] ){
        ret<-ret*K_Y_X(yi=other_obs_Y[k],xi=other_obs_X[k,],theta=x) 
      }
      return(ret)
    } 
  }else{
    integrand <- function(x) { K_Y_X(yi=yi,xi=xi,theta=x ) } 
  }
  #"x" is vector and can be fed with varying length
  int<-adaptIntegrate(integrand, lowerLimit = rep(-1e5,1), upperLimit = rep(1e5,1),absError=1e-6)
  return(int$integral)
}
#We can also compute the marginal
K_X<-function(xi,psi){
  mu.i=rep(psi,p);
  sigma.xi=diag(p);
  #require(mvtnorm)
  return(dmvnorm(xi,mean=mu.i,sigma=sigma.xi))
  #Use following code to check K_X is indeed a kernel.
}
h_X<-function(xi,other_obs=matrix()){
  require(cubature)
  if(is.null(dim(other_obs)) ){
    integrand <- function(x) {
      ret<-1;
      for(k in 1:dim(other_obs)[1] ){
        ret<-ret*K_X(xi=other_obs[k,],psi=x) 
      }
      return(ret)
    } 
  }else{
    integrand <- function(x) {K_X(xi=xi,psi=x)}
  }
  #"x" is vector and can be fed with varying length
  int<-adaptIntegrate(integrand, lowerLimit = -1e5, upperLimit = 1e5,absError=1e-6)
  return(int$integral)
}
#We can also compute the marginal
###
##########################################################################################

#Following is the MCMC algorithm that is implemented in [Wade et.al].
##########################################################################################
### A single step in MCMC simulation.
MCMC<-function(rho,k_j,x.data,y.data,
               #rho, k_j are two parameters of the hierachical model we want to keep track of in each MCMC step.
               #x.data and y.data is used for inputing dataset.
               theta=null,psi=null,method=EDP,iterparam=100){
  ##########
  #According to the analyzed result, we define two matrices to store parameters for each dimension.
  #theta is associted with Y, which has a coarser partition.
  #psi is associated with X, which is conditioned on the value of Y.
  x.data.new<-x.data;
  y.data.new<-y.data;
  rho.new<-rho;
  #Step 1: Update the partition structure for each leave-one-out observation.
  for(i in 1:n){
    x.loop<-x.data[i,]
    y.loop<-y.data[i]
    x.data.remove.i<-x.data[-i,];
    y.data.remove.i<-y.data[-i];
    record.i<-rho[,i];
    rho.remove.i<-rho[,-i];
    k_j.remove.i<-analyze_partition(rho.remove.i)$k_j;
    #####
    y.part<-record.i[1];
    x.part<-record.i[2];
    #print(record.i)
    #print(k_j.remove.i)
    if(is.na(k_j.remove.i[y.part,x.part])){next;}
    if(k_j.remove.i[y.part,x.part]==1){
      k_j.remove.i[y.part,x.part]<-NA;
      #if(sum(!is.na(k_j.remove.i[y.part,]))==0){
      #  k_j.remove.i<-k_j.remove.i[-y.part,];
      #}
      #This if structure is not necessarily needed since when the whole row is NA, then the probability given by omega is the same as new cluster.
    }else{
      k_j.remove.i[y.part,x.part]<-k_j.remove.i[y.part,x.part]-1;
    }
    
    
    #We analysed rho after removing the i-th column.
    #We sample s_i below using the CRP with the weight function omega.
    #
    
    #Again, according to what we sampled we update the nested partition structure on parameter.
    #
    alpha_psi<-function(psi.val=1){
      return(1/pi);
      #return(dexp(psi.val,rate=1)) This return if we want to put a hyperprior on theta as indicated in [Wade et.al]
    }
    alpha_theta<-function(theta.val,condition.psi){
      theta.val<-1/theta.val;
      if(is.na(condition.psi) | is.na(theta.val)){return(0)}
      return(dnorm(theta.val,mean=condition.psi,sd=1))
      #intensity rate for Dirichlet process to generate a new clutser on y|X.
    }
    omega<-function(j,l,yi,xi,theta1=theta,psi1=psi,k_j=k_j.remove.i,rho=rho.remove.i){
      require(mvtnorm)
      #message("omega_",j,'_',l)
      #This is the weight function such that the observation (yi,xi) will fall into the cluster
      #(j,l) for y clusters and x clusters respectively.
      theta.val<-theta[j];
      psi.val<-psi[j,l];
      n_j<-sum(k_j[j,!is.na(k_j[j,])]);
      if(!is.na(k_j[j,l])){
        #put it in existing x cluster.
        ret<-1;
        ret<-ret*n_j*sum(!is.na(k_j[j,]))/(n_j+alpha_theta(theta.val,condition.psi=psi.val));
        ret<-ret*K_Y_X(yi=yi,xi=xi,theta=theta.val);
        ret<-ret*K_X(xi=xi,psi=psi.val);
      }
      if(is.na(k_j[j,l]) & sum(!is.na(k_j[j,]))>0 ){
        #create new x cluster within existing y cluster.
        ret<-1;
        ret<-ret*n_j*alpha_theta(theta.val,condition.psi=psi.val)/(n_j+alpha_theta(theta.val,condition.psi=psi.val));
        ret<-ret*K_Y_X(yi=yi,xi=xi,theta=theta.val);
        Indx<-rho[,rho[1,]==j & rho[2,]==l];
        ret<-ret*h_X(xi=xi,other_obs = x.data.remove.i[Indx,]);
      }
      if(is.na(k_j[j,l]) & sum(!is.na(k_j[j,]))<=0 ){
        #create new x cluster within new y cluster.
        ret<-1;
        Indx<-rho[,rho[1,]==j & rho[2,]==l];
        ret<-alpha_theta(theta.val,condition.psi = psi.val)*h_Y_X(yi=yi,xi=xi,other_obs_X =x.data.remove.i[Indx,], other_obs_Y = y.data.remove.i[Indx,])*
          h_X(xi,other_obs = x.data.remove.i[Indx,]);
      }
      return(ret);
    }
    
    #This sampling is more like doing CRP(Chinese restaurant process) twice, with weights being omegas.
    possible.cluster<-matrix(NA,nrow=n*n+1,ncol=3);
    counter<-1;
    for(recluster.y in 1:(n-1) ){
      if(sum(!is.na(k_j.remove.i[recluster.y,]))>0 ){
        xseq<-which(!is.na(k_j.remove.i[recluster.y,])); 
        #New x.cluster is always the last one.
        #print(recluster.y)
        for(recluster.x in c(xseq,max(xseq)+1) ){
          possible.cluster[counter,1]<-recluster.y;
          possible.cluster[counter,2]<-recluster.x;
          possible.cluster[counter,3]<-omega(j=recluster.y,l=recluster.x,yi=y.loop,xi=x.loop,psi1 = psi,theta1 = theta,
                                             k_j=k_j.remove.i,rho=rho.remove.i);
          counter=counter+1;
        }
      }else{
        possible.cluster[counter,1]<-recluster.y;
        possible.cluster[counter,2]<-1;
        possible.cluster[counter,3]<-omega(j=recluster.y,l=1,yi=y.loop,xi=x.loop,psi1 = psi,theta1 = theta,
                                           k_j=k_j.remove.i,rho=rho.remove.i);
        #recluster.x=1 because we want to create a new partition for y.
        counter=counter+1;
      }
    }
    possible.cluster<-possible.cluster[complete.cases(possible.cluster),]
    norm.prob<-possible.cluster[,3]/sum(possible.cluster[,3])
    norm.prob[is.nan(norm.prob)]<-0;
    if(sum(norm.prob)==0){
      #message(i)
      # We do not want to change to another possible cluster.
    }else{
      #normalized.prob<-as.numeric(as.character(normalized.prob))
      #np=round(normalized.prob,6)
      new.obs<-sample(1:length(norm.prob),size=1,prob=norm.prob,replace=T);
      #Here we must add 'replcae=T' to avoid 'too few positive probabilities' error.
      #message(i)
      
      new.record<-possible.cluster[new.obs,]
      rho.new[,i]<-new.record[1:2];
    }
    #print(possible.cluster)
    k_j.new<-analyze_partition(rho.new)$k_j
    
  }
  if(method=="Image"){
    #####
    #Plot module 'EDP_visual.R' must be put under the same directory.
    source('EDP_visual.R')
    matdf<-as.data.frame(matrix2list(k_j))
    grp<-ggplot(data=matdf)+
      geom_text(aes(x=col,y=row),
                label=as.character(matdf$value) )+
      scale_x_continuous(breaks = 1:dim(k_j)[2])+
      scale_y_continuous(breaks = 1:dim(k_j)[1])
    #####
    print(grp)
  }
  ##########
  print("After Step1")
  print(k_j);
  ##########
  #Step 2: Metropolist-Hastings step.
  #MH step, not necessarily, but improves mixing greatly. Since our data is on a small scale, this is not necessary and is a repeatation of [Neal 2000] earlier work.
  #Move 1
  MHadd<-function(rho,k_j){
    rho.saved<-rho;
    k_j.saved<-k_j;
    count_partitions<-function(k_j){
      y_cluster_more_than_one<-c();
      y_cluster_only_one<-c();
      y_cluster_count<-rep(0,n)
      for(a in 1:dim(k_j)[1]){
        y_cluster_count[a]<-0;
        for(b in 1:dim(k_j)[2]){
          if(!is.na(k_j[a,b])){
            #message(a,";",b);
            y_cluster_count[a]<-y_cluster_count[a]+1;
          }
        }
        if(y_cluster_count[a]>1){
          y_cluster_more_than_one<-c(y_cluster_more_than_one,a)
        }else{
          if(y_cluster_count[a]==1){
            y_cluster_only_one<-c(y_cluster_only_one,a)
          }
        }
      }
      ret<-list(mto=y_cluster_more_than_one,ono=y_cluster_only_one)
      return(ret)
    }
    cps<-count_partitions(k_j)
    y_cluster_more_than_one<-cps$mto
    y_cluster_only_one<-cps$ono
    message("Move 1")
    if(length(y_cluster_more_than_one)>1){
      if(length(y_cluster_more_than_one)>1){
        c<-sample(y_cluster_more_than_one,size=1)
      }else{
        c<-y_cluster_more_than_one;
      }
      x_cluster_in_y<-which(!is.na(k_j[c,]))
      if(length(x_cluster_in_y)>1){d<-sample(x_cluster_in_y,size=1)}else{d<-x_cluster_in_y}
      candidate_y_cluster<-c(y_cluster_more_than_one,y_cluster_only_one);
      candidate_y_cluster<-candidate_y_cluster[candidate_y_cluster!=c];
      #Different y cluster.
      if(length(candidate_y_cluster)>1){ newc<-sample(candidate_y_cluster,size=1)}else{newc<-candidate_y_cluster}
      seqd<-1:n;
      seqd<-seqd[!seqd %in% (which(complete.cases(k_j[newc,])))];
      if(length(seqd)>1){
        newd<-sample( seqd,size=1 )
      }else{
        newd<-seqd
      }
      message("(",c,",",d,")","->","(",newc,",",newd,")")
      #print(k_j)
      #print(x_cluster_in_y)
      k_j[newc,newd]<-k_j[c,d]
      
      #update rho
      for(enum in which(rho[1,]==c & rho[2,]==d)){
        rho[1,enum]<-newc;
        rho[2,enum]<-newd;
      }
      k_j[c,d]<-NA
    
      #Proposal Move 1 accepted?
      pmove1<-1;
      n_j<-sum(k_j.saved[c,][complete.cases(k_j.saved[c,])])
      n_l_j<-k_j.saved[c,d]
      n_h<-sum(k_j.saved[newc,][complete.cases(k_j.saved[newc,])])
      message(n_j,"&",n_l_j,"&",n_h,"@",c)
      pmove1<-pmove1*gamma(n_j-n_l_j)*gamma(n_h+n_l_j)/(gamma(n_j)*gamma(n_h));
      pmove1<-pmove1*gamma(alpha_psi(theta[c])+n_j)*gamma(alpha_psi(theta[newc])+n_h)
      pmove1<-pmove1/gamma(alpha_psi(theta[c])+n_j-n_l_j)
      pmove1<-pmove1/gamma(alpha_psi(theta[newc])+n_h+n_l_j)
      pmove1<-pmove1*alpha_psi(theta[newc])/alpha_psi(theta[c])
      pmove1<-pmove1*length(y_cluster_more_than_one)/length(count_partitions(k_j)$mto)
      #The second k_j has already been updated. Move 1 is possible only when there are at least 2 more-than-one psi cluster.
      candidate_obs<-which(rho[1,]==c & rho[2,]==d);
      multiplier<-1;
      for(A in candidate_obs){
        multiplier<-multiplier*K_Y_X(xi=x.data[A],yi=y.data[A],theta = theta[newc])
        multiplier<-multiplier/K_Y_X(xi=x.data[A],yi=y.data[A],theta = theta[newc])
      }
      pmove1<-pmove1*multiplier
      
      if(runif(1)<min(pmove1,1)){
        message("Move 1 proposal accepted.=",pmove1)
        k_j.saved<-k_j;
        rho.saved<-rho;
      }else{
        message("Move 1 proposal rejected.=",pmove1)
        k_j<-k_j.saved
        rho<-rho.saved
      }
    }
    cps<-count_partitions(k_j)
    y_cluster_more_than_one<-cps$mto
    y_cluster_only_one<-cps$ono
    if(runif(1)<0.5){
      #Move 2
      message("Move 2")
      if(length(y_cluster_more_than_one)>1){
        cc<-sample(y_cluster_more_than_one,size=1);
      }else{
        cc<-y_cluster_more_than_one;
      }
      seqcc<-1:n;
      seqcc<-seqcc[! seqcc %in% y_cluster_more_than_one];
      seqcc<-seqcc[! seqcc %in% y_cluster_only_one];
      #New y cluster.
      if(length(seqcc)>1){
        newcc<-sample(seqcc,size=1);
      }else{
        newcc<-seqcc
      }
      choosedd<-which(complete.cases(k_j[cc,]))
      if(length(choosedd)>1){
        dd<-sample(choosedd,size=1)
      }else{
        dd<-choosedd;
      }
      
      
      message("(",cc,",",dd,")","->","(",newcc,",1)")
      #print(k_j)
      k_j[newcc,1]<-k_j[cc,dd]
      #update rho
      for(enum in which(rho[1,]==cc & rho[2,]==dd)){
        rho[1,enum]<-newcc;
        rho[2,enum]<-1;
      }
      k_j[cc,dd]<-NA;
      #Proposal Move 2 accepted?
      pmove2<-1;
      n_j<-sum(k_j.saved[cc,][complete.cases(k_j.saved[cc,])])
      n_l_j<-k_j.saved[cc,dd]
      n_h<-sum(k_j.saved[newcc,][complete.cases(k_j.saved[newcc,])])
      message(n_j,"&",n_l_j,"&",n_h,"@",cc,":",dd)
      #print(k_j)
      pmove2<-pmove2*gamma(n_j-n_l_j)*gamma(n_l_j)/gamma(n_j);
      
      pmove2<-pmove2*gamma(alpha_psi(theta[cc])+n_j)*gamma(alpha_psi(theta[newcc]) )
      pmove2<-pmove2/gamma(alpha_psi(theta[cc])+n_j-n_l_j)
      pmove2<-pmove2/gamma(alpha_psi(theta[newcc])+n_l_j)
      pmove2<-pmove2*alpha_psi(theta[newcc])/alpha_psi(theta[cc])*alpha_theta(theta.val =theta[cc] ,condition.psi =psi[cc,dd])
      pmove2<-pmove2*length(y_cluster_more_than_one)/length(count_partitions(k_j)$ono)
      pmove2<-pmove2/(length(count_partitions(k_j)$ono)+length(count_partitions(k_j)$mto)-1)
      #The second k_j has already been updated. Move 1 is possible only when there are at least 2 more-than-one psi cluster.
      candidate_obs<-which(rho[1,]==cc & rho[2,]==dd);
      multiplier<-1;
      
      for(A in candidate_obs){
        other_seq<-candidate_obs[candidate_obs!=A];
        multiplier<-multiplier*h_Y_X(xi=x.data[A],yi=y.data[A],other_obs_X=x.data[other_seq,],other_obs_Y=y.data[other_seq,] ,theta = theta[newcc])
        multiplier<-multiplier/K_Y_X(xi=x.data[A],yi=y.data[A],theta = theta[cc])
      }
      pmove2<-pmove2*multiplier
      
      if(runif(1)<min(pmove2,1)){
        message("Move 2 proposal accepted.=",pmove2)
        k_j.saved<-k_j;
        rho.saved<-rho;
      }else{
        message("Move 2 proposal rejected.=",pmove2)
        k_j<-k_j.saved
        rho<-rho.saved
      }
    }else{
      #Move 3
      message("Move 3")
      if(length(y_cluster_only_one)>0){
        if(length(y_cluster_only_one)>1){
          ccc<-sample(y_cluster_only_one,size=1);
        }else{
          ccc<-y_cluster_only_one;
        }
        ddd<-1;
        for(notNA in 1:dim(k_j)[2]){
          if(!is.na(k_j[ccc,notNA])){ddd<-notNA;}
        }
        candidate_y_cluster<-c(y_cluster_more_than_one,y_cluster_only_one);
        candidate_y_cluster<-candidate_y_cluster[candidate_y_cluster!=ccc];
        if(length(candidate_y_cluster)>1){
          newccc<-sample(candidate_y_cluster,size=1)
        }else{
          newccc<-candidate_y_cluster;
        }
        #Different y cluster.
        seqddd<-1:n
        seqddd<-seqddd[!seqddd %in% (which(complete.cases(k_j[newccc,])))];
        if(length(seqddd)>1){
          newddd<-sample(seqddd ,size=1 )
        }else{
          newddd<-seqddd
        }
        
        
        message("(",ccc,",",ddd,")","->","(",newccc,",",newddd,")")
        print(k_j)
        
        k_j[newccc,newddd]<-k_j[ccc,ddd]
        #update rho
        for(enum in which(rho[1,]==ccc & rho[2,]==ddd)){
          rho[1,enum]<-newccc;
          rho[2,enum]<-newddd;
        }
        k_j[ccc,]<-NA;
        #Proposal Move 3 accepted?
        pmove3<-1;
        n_j<-sum(k_j.saved[ccc,][complete.cases(k_j.saved[ccc,])])
        n_l_j<-k_j.saved[ccc,ddd]
        n_h<-sum(k_j.saved[newccc,][complete.cases(k_j.saved[newccc,])])
        message(n_j,"&",n_l_j,"&",n_h,"@",ccc,":",ddd)
        #print(k_j)
        pmove3<-pmove3*gamma(n_h+n_l_j)/( gamma(n_l_j)*gamma(n_h) );
        pmove3<-pmove3*gamma(alpha_psi(theta[ccc])+n_l_j)*gamma(alpha_psi(theta[newccc]+n_h) )
        pmove3<-pmove3/gamma(alpha_psi(theta[ccc]))
        pmove3<-pmove3/gamma(alpha_psi(theta[newccc])+n_h+n_l_j)
        pmove3<-pmove3*alpha_psi(theta[newccc])/alpha_psi(theta[ccc])*(1/alpha_theta(theta.val =theta[ccc] ,condition.psi =psi[ccc,ddd]))
        pmove3<-pmove3/length(count_partitions(k_j)$mto)
        pmove3<-pmove3*(length(y_cluster_only_one)*(length(y_cluster_only_one)+length(y_cluster_more_than_one))-1)
        #The second k_j has already been updated. Move 1 is possible only when there are at least 2 more-than-one psi cluster.
        candidate_obs<-which(rho[1,]==ccc & rho[2,]==ddd);
        multiplier<-1;
        
        for(A in candidate_obs){
          other_seq<-candidate_obs[candidate_obs!=A];
          multiplier<-multiplier/h_Y_X(xi=x.data[A],yi=y.data[A],other_obs_X=x.data[other_seq,],other_obs_Y=y.data[other_seq,] ,theta = theta[newccc])
          multiplier<-multiplier*K_Y_X(xi=x.data[A],yi=y.data[A],theta = theta[newccc])
        }
        pmove3<-pmove3*multiplier
        
        if(runif(1)<min(pmove3,1)){
          message("Move 3 proposal accepted.=",pmove3)
          k_j.saved<-k_j;
          rho.saved<-rho;
        }else{
          message("Move 3 proposal rejected.=",pmove3)
          k_j<-k_j.saved
          rho<-rho.saved
        }
      }
    }
    #print(k_j)
    ret<-list(rhoadd=rho,k_jadd=k_j)
    return(ret)
  }
  MHaddStep<-MHadd(rho.new,k_j.new)
  rho.new<-MHaddStep$rhoadd
  k_j.new<-analyze_partition(rho.new)$k_j
  #print("after MH")
  #print(rho.new)

  
  if(method=="ImageMH"){
    #####
    #Plot module 'EDP_visual.R' must be put under the same directory.
    source('EDP_visual.R')
    matdf<-as.data.frame(matrix2list(k_j.new))
    grp<-ggplot(data=matdf)+
      geom_text(aes(x=0,y=0),label="MHafter")+
      geom_text(aes(x=col,y=row),
                label=as.character(matdf$value) )+
      scale_x_continuous(breaks = 1:dim(k_j)[2])+
      scale_y_continuous(breaks = 1:dim(k_j)[1])
    #####
    print(grp)
  }
  ##########
  print("After Step2")
  print(k_j);
  ##########
  #Step 3: Update cluster parameters.
  if(TRUE){
    #This is a null stub if.
    #MCMC for theta
    #
    likelihood.theta <- function(y=y.data,x=x.data,theta_j,j){
      #likelihood of theta_{j} is only involving the observations in category j.
      sumll<-0;
      pre.data.x<-x[rho.new[1,]==j,];
      pre.data.y<-y[rho.new[1,]==j];
      if(all(is.na(pre.data.y))){return(sumll)}
      if(all(is.na(pre.data.x))){return(sumll)}
      
      if(is.null(dim(pre.data.x)[1])){
        iter.upper=-1;
      }else{
        iter.upper=dim(pre.data.x)[1];
      }
      if(iter.upper!=-1){
        for(i in 1:iter.upper) {
          
          singlelikelihood = K_Y_X(yi=pre.data.y[i], xi=pre.data.x[i,], theta=theta_j);
          sumll = sumll+singlelikelihood;
        }
      }
      return(sumll)   
    }
    prior.theta <- function(theta_j){
      sumpr<-0;
      for(i in 1:length(theta_j)){
        singleprior = dgamma(theta_j[i], shape=2.5, rate=1, log = T);
        sumpr<-sumpr+singleprior;
      }
      return(sumpr)
    }
    posterior.theta <- function(y=y.data,x=x.data,theta_j,j){
      return (likelihood.theta(y=y.data,x=x.data,theta_j,j) + prior.theta(theta_j))
    }
    proposalfunction.theta <- function(theta){
      return(rmvnorm(1,mean = theta, sigma=diag(length(theta))))
    }
    
    run_metropolis_MCMC_theta <- function(startvalue, iterations=iterparam){
      chain = array(dim = c(iterations+1,length(startvalue)))
      chain[1,] = startvalue
      for (iter in 1:iterations){
        proposal = proposalfunction.theta(chain[iter,])
        for(j in 1:length(startvalue)){
          logicalIndex<-rho.new[1,]==j;
          y_in_j<-y.data[logicalIndex];
          x_in_j<-x.data[logicalIndex,];
          probab = exp(posterior.theta(x=x_in_j,y=y_in_j,theta=proposal[j],j=j) - 
                         posterior.theta(x=x_in_j,y=y_in_j,theta=chain[iter,j],j=j))
          if (runif(1) < probab){
            chain[iter+1,j] = proposal[j]
          }else{
            chain[iter+1,j] = chain[iter,j]
          }
        }
      }
      return(chain)
    }
    #MCMC for psi
    #
    likelihood.psi <- function(y=y.data,x=x.data,psi_l_j,l,j){
      sumll<-0;
      pre.data<-x[rho.new[2,]==l & rho.new[1,]==j,];
      if(all(is.na(pre.data))){return(sumll)}
      if(is.null(dim(pre.data)[1])){
        iter.upper=-1;
      }else{
        iter.upper=dim(pre.data)[1];
      }
      if(iter.upper!=-1){
        for(i in 1:iter.upper){
          singlelikelihood = K_X(xi=x[i,], psi=psi_l_j);
          sumll = sumll+singlelikelihood;
        }
      }
      return(sumll)   
    }
    prior.psi <- function(psi_l_j){
      sumpr<-0;
      for(i in 1:length(psi_l_j)){
        singleprior = dgamma(psi_l_j[i], shape=2.5, rate=1, log = T);
        sumpr<-sumpr+singleprior;
      }
      return(sumpr)
    }
    posterior.psi <- function(y=y.data,x=x.data,psi_l_j,j,l){
      return (likelihood.psi(y=y.data,x=x.data,psi_l_j,j,l) + prior.psi(psi_l_j))
    }
    proposalfunction.psi <- function(psi){
      return(rmvnorm(1,mean = psi, sigma=diag(length(psi))))
    }
    
    run_metropolis_MCMC_psi <- function(startvalue, iterations=iterparam,j){
      #lseq<-which(!is.na(k_j[j,]))
      chain = array(dim = c(iterations+1,length(startvalue)))
      #if(length(startvalue)!=length(lseq)){message('dimension mismatch!');return(NA);}
      for(k in 1:length(startvalue)){chain[1,k] = startvalue[k]}
      for (iter in 1:iterations){
        proposal = proposalfunction.psi(chain[iter,])
        for(l in 1:length(startvalue)){
          #Loops observations in the same j,l cell.
          y_in_j_l<-y.data[rho.new[1,]==j & rho.new[2,]==l,];
          x_in_j_l<-x.data[rho.new[1,]==j & rho.new[2,]==l,];
          probab = exp(posterior.psi(x=x_in_j_l,y=y_in_j_l,psi=proposal[j],j=j,l=l) - 
                         posterior.psi(x=x_in_j_l,y=y_in_j_l,psi=chain[iter,j],j=j,l=l))
          if (runif(1) < probab){
            chain[iter+1,l] = proposal[l]
          }else{
            chain[iter+1,l] = chain[iter,l]
          }
        }
      }
      return(chain)
    }
  }
  
  #eliminate those parameter not in play.
  chain_theta<-run_metropolis_MCMC_theta(startvalue = theta);
  theta.draw<-tail(chain_theta,1);
  theta.draw<-theta.draw[!is.na(theta.draw)]
  theta[1:length(theta.draw)]<-theta.draw;
  #Update parameter theta.
  
  for(row in 1:dim(psi)[1]){
    psi.sv<-psi[row,]
    chain_psi<-run_metropolis_MCMC_psi(startvalue=psi[row,], iterations=iterparam,j=row)
    psi.draw<-tail(chain_psi,1);
    psi.draw<-psi.draw[!is.na(psi.draw)]
    psi[row,1:length(psi.draw)]<-psi.draw;
  }
  #Update parameter psi.
  #k.add=analyze_partition(rho.add)$k
  #k_j.add=analyze_partition(rho.add)$k_j
  #update k
  #print(rho.new);
  ##########
  returnList<-list(
    theta=theta,
    psi=psi,
    rho=rho.new,
    x.data=x.data,
    y.data=y.data,
    k_j=k_j.new)
  #Make a return list.
  
  return(returnList);
}
# The function 'MCMC' only realize one step in MCMC algorithm of EDP.
# The computational time is extremely long...
### Now we do a MCMC in-action.
##########################################################################################

##########################################################################################
### MCMC simulation with totalSteps length chain.
psi.tmp<-psi;
theta.tmp<-theta;
rho.tmp<-rho;
k_j.tmp<-k_j;
x.tmp<-x.data;
y.tmp<-y.data;

fileName<-paste("EDP_SIM_",Sys.Date(),"_",rnorm(1),".pdf",sep="");



pdf(fileName)
source('EDP_visual.R')
analyze_plot(rho,x.data,y.data)

pb <- txtProgressBar(min = 0, max = totalSteps, style = 3)
time0<-proc.time()
for(mcmc_steps in 1:totalSteps){
  #Let us pretend 101 steps leads to a perfect mixing..
  each_stepList<-MCMC(rho=rho.tmp,k_j=k_j.tmp,
                      x.data=x.tmp,y.data=y.tmp,
                      psi=psi.tmp,theta=theta.tmp,method="Image",iterparam = MHSteps);
  rho.tmp<-each_stepList$rho;
  psi.tmp<-each_stepList$psi;
  theta.tmp<-each_stepList$theta;
  x.tmp<-each_stepList$x.data;
  y.tmp<-each_stepList$y.data;
  k_j.tmp<-each_stepList$k_j;
  setTxtProgressBar(pb,mcmc_steps)
  
}
close(pb)
proc.time()-time0;
###Assume that the trace plots are perfect...
textplot(rho.tmp)
textplot(k_j.tmp)
psi.tmp[is.na(k_j.tmp)]<-NA
textplot(psi.tmp)
theta.tmp<-theta.tmp[rowSums(!is.na(k_j.tmp))>0]
textplot(theta.tmp)
analyze_plot(rho.tmp,x.data,y.data)
dev.off()
#In order to interpret psi and theta, you need to match them with the partition rho and throw those corresponding to NA.
##########################################################################################
