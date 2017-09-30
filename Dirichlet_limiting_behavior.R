N<-100;
theta<-rbeta(N,1,1)
P<-c(theta[1]);
for(j in 2:N){
  P[j]<-(1-P[j-1])*theta[j];
}
indicator<-function(a,b,V){
  #decide whether the value of random variable V is in (a,b]
  if(V>a && V<=b){
    return(1);
  }else{
    return(0);
  }
}
V<-rcauchy(N);
#We choose N(0,1)dx as the concentration measure on the real line.
DirichletP<-function(startpoint,endpoint){
  result<-0;
  for(i in 1:N){
    result<-result+P[i]*indicator(startpoint,endpoint,V[i])
  }
  return(result);
}
#Now we mesh the grid as very small to show how Dirichlet Process behaves.
DirichletP.value<-as.data.frame(matrix(NA,ncol=3));
colnames(DirichletP.value)<-c("start","end","value")
counter<-0;
for(mesh.startpoint in seq(-3,3,0.05)){
  for(mesh.endpoint in seq(mesh.startpoint,3,0.05)){
    DirichletP.value[counter,]<-c(mesh.startpoint,
                                 mesh.endpoint,
                                 DirichletP(mesh.startpoint,mesh.endpoint)
                                 )
    counter<-counter+1;
  }  
}
library(rgl)
plot3d(DirichletP.value$start, DirichletP.value$end, DirichletP.value$value,col=DirichletP.value$end-DirichletP.value$start+1) 
DirichletP.cdf<-function(x){
  return(DirichletP(1e-10,x))
}
DirichletP.density<-function(x,epsilon=1e-1){
  return(DirichletP.cdf(x+epsilon)-DirichletP.cdf(x-epsilon))/(2*epsilon)
}
DirichletP.plotcdf<-function(range=6){
  x<-seq(0,range,1e-3);
  y<-c();
  counter<-1;
  for(increment in x){
    y[counter]<-DirichletP.cdf(increment);
    counter<-counter+1;
  }
  plot(x,y,type="l")
}
DirichletP.plotpdf<-function(range=6){
  x<-seq(0,range,1e-3);
  y<-c();
  counter<-1;
  for(increment in x){
    y[counter]<-DirichletP.density(increment);
    counter<-counter+1;
  }
  plot(x,y,type="s")
}
par(mfrow=c(1,2));
DirichletP.plotpdf()
DirichletP.plotcdf()
Dirichlet.mixture.plotpdf<-function(howmanyDirichlets=2){
  x<-seq(0,range,1e-3);
  y<-c();
  counter<-1;
  for(increment in x){
    pointvalue<-0;
    breaks<-c(sort(runif(howmanyDirichlets-1,0,1)),1);
    parameters<-c(breaks[1],diff(breaks));
    for(convex.parameter in parameters){
      pointvalue<-DirichletP.density(increment)*convex.parameter
    }
    y[counter]<-pointvalue;
    counter<-counter+1;
  }
  plot(x,y,type="s")
}