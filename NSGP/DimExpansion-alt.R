#The data set I used is the same as
#Bornn, Luke, Gavin Shaddick, and James V. Zidek. "Modeling nonstationary processes through dimension expansion." Journal of the American Statistical Association 107.497 (2012): 281-289.
#Section 3.1(Solar radiation data set)
#Hay, John E. "An assessment of the mesoscale variability of solar radiation at the earth's surface." Solar Energy 32.3 (1984): 425-434.
#According to a personal email correspondence with Mr.Bornn, he directed me to p.115 of
#Sampson, P. D., & Guttorp, P. (1992). Nonparametric estimation of nonstationary spatial covariance structure. Journal of the American Statistical Association, 87(417), 108-119.
#and said they never used the full solar data set from Hay's paper.
#I reproduce this solar radiation data set from [Sampson&Guttorp]p.115
#dispersion matrix for solar stations
A<-NA;
solar_dispersion<-as.matrix(rbind(
                            c( A,34,36,50,41,50,63,38,42,47,49,49),
                            c(34, A, 5,14, 9,15,26,16,18,31,31,29),
                            c(36, 5, A,12, 6,12,25,14,17,31,31,29),
                            c(50,14,12, A, 7, 4,12,22,21,41,36,32),
                            c(41, 9, 6, 7, A, 6,17,14,15,33,30,27),
                            c(50,15,12, 4, 6, A,10,19,17,37,31,29),
                            c(63,26,25,12,17,10, A,26,24,40,34,31),
                            c(38,16,14,22,14,19,26, A, 6,14,14,14),
                            c(42,18,17,21,15,17,24, 6, A,15,11,10),
                            c(47,31,31,41,33,37,40,14,15, A, 6, 9),
                            c(49,31,31,36,30,31,34,14,11, 6, A, 3),
                            c(49,29,29,32,27,29,31,14,10, 9, 3, A)
                            )
                            )
colnames(solar_dispersion)<-c(1,2,3,12,11,4,5,6,10,7,8,9)
rownames(solar_dispersion)<-c(1,2,3,12,11,4,5,6,10,7,8,9)

v_star<-function(i,j){
  return(solar_dispersion[as.character(i),as.character(j)])
}
#L^2 norm function
#gamma function
gamma<-function(d,phi1,phi2,phi3){
  res<-0;
  res<-phi1*( 1-exp(-d/phi2) )+phi3
  return(res);
}

#objective function
OBJ<-function(phi1,phi2,phi3,Z,lambda1,X=c(0,0)){
  res<-0;
  #X is the known coordinates of locations in lower dimension from [Hay]1984.
  augmented_vec<-c(X,Z);
  augmented_dist<-as.matrix( dist(augmented_vec,p=2,diag=T) )
  for(j in 2:3){
    #2:3 because the augmented Z can only add one more dimension in this case.
    for(i in 1:(j-1) ){
      #message(i,j)
      temp_dist_ij<-augmented_dist[i,j]
      res<-res+( v_star(i,j)-
                 gamma(temp_dist_ij,phi1,phi2,phi3)
                )^2
    }
  }
  #Penalize by L1 norm of the expanded dimension.
  res<-res+lambda1*sum(abs(Z))
  return(res)
}


OBJ05<-function(vec,X){
  #This is a interface function for optim procedure.
  #We must write such an interface because optim procedure only optimize one vector parameter at one run.
  #vec[1]~vec[3] correspond to phi1~phi3
  #vec[4] correspond to Z
  return(OBJ(vec[1],vec[2],vec[3],vec[4],lambda1=.2,X))
  #lambda=0.2 is chosen from [Bornn et.al]2012
}
EXPAND<-function(location=c(0,0)){
  #Specify a location using observed data points in lower dimension.
  OPT<-optim(par=c(1,1,1,1),fn=OBJ05,X=location,method="Nelder-Mead")
  #According to [Bornn]2012, Nelder-Mead and BFGS works well for small amount of observations.
  #optimiser phi
  phi_hat<-OPT$par[1:3]
  #optimiser Z
  Z_hat<-OPT$par[4]
  #optimiser augmented [X,Z]
  expanded_coordinate<-c(location,Z_hat)
  return(expanded_coordinate)
}
#Now we consider the real dataset from [Hay]1984, Table 2(p.429)
locations2D<-as.data.frame(rbind(
              c(1,"UBC CLI. STN",49.16,123.15,93),
              c(2,"Langara",49.13,123.06,63),
              c(3,"Vancouver Airport",49.11,123.10,5),
              c(4,"B.C. Hydro",49.16,123.07,122),
              c(5,"Tsawwassen Ferry",49.00,123.08,3),
              c(6,"North Mount",49.19,123.04,114),
              c(7,"Grouse Mountain",49.23,123.05,1128),
              c(8,"Pitt Meadows Airport",49.13,122.42,5),
              c(9,"Mission - Habitat Apts",49.08,122.17,125),
              c(10,"Abbotsford (Library)",49.02,122.17,60),
              c(11,"Abbotsford Airport",49.01,122.22,61),
              c(12,"Langley",49.06,122.38,11) 
             )
             )
colnames(locations2D)<-c("num","station","lat","long","elev")


x_coord<-as.numeric(as.character(locations2D$lat))
y_coord<-as.numeric(as.character(locations2D$long))
elev_real<-as.numeric(as.character(locations2D$elev))
library(fields)
fit<- Tps(cbind(x_coord,y_coord), elev_real);
elev_1<-as.vector(predict(fit));

elev_2<-x_coord;
for(a in 1:length(x_coord)){
  elev_2[a]<-EXPAND(c(x_coord[a]/1000,y_coord[a]/1000))[3]*1000;
}


d<-sqrt(x_coord^2+y_coord^2)
par(mfrow=c(3,1))
plot(d,elev_real,ylim=c(0,1200))
plot(d,elev_1   ,ylim=c(0,1200))
plot(d,elev_2   ,ylim=c(0,1200))
plot.new()
par(mfrow=c(1,1))
par(new=T)
plot  (d,elev_real,ylim=c(0,1200),col=rgb(1,0,0,alpha=0.5),pch=16,cex=1.2,ylab="Estimate elevation",xlab="distance")
points(d,elev_1   ,ylim=c(0,1200),col=rgb(0,1,0,alpha=0.5),pch=16,cex=1.2)
points(d,elev_2   ,ylim=c(0,1200),col=rgb(0,0,1,alpha=0.5),pch=16,cex=1.2)
