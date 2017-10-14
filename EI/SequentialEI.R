# Sequential Design using Expected Improvement
library(DiceOptim)
#library(rgl)
#It is impossible to visualize 6-dim parameter space.
#get our initial starting design/from last run
csvnew<-read.csv('appendstudy1.csv',header=T)
csvnew<-csvnew[complete.cases(csvnew),]
design=csvnew[,2:7]
y.branin=csvnew[,8]
# Fit the GP model - built into the DiceOptim package

fit.gp=km(~1+X1+X2+X3+X4+X5+X6, design = design, 
          response = y.branin,covtype = "gauss")
#fit the model with the already observed data, each time we observe a new observation from RMIII,
#https://www.grc.nasa.gov/www/k-12/rocket/rktsim.html
#we append it to appendstudy1.csv for updating our GP model.
Ngrid=5
#Given the space-filling design is relatively sparse in 6-dim space, a finer grid might improve EI
#However due to limitation of computation power, we restrict ourselves to 5^6 grid, which seems good enough.

X=as.matrix(expand.grid(seq(0,1,length=Ngrid)*(25.147-0.127)+0.127,
                        seq(0,1,length=Ngrid)*(181.051-0.00)+0.000,
                        seq(0,1,length=Ngrid)*(12.575-0.249)+0.249,
                        seq(0,1,length=Ngrid)*(12.585-1.27)+1.27,
                        seq(0,1,length=Ngrid)*(6.01-1.27)+1.27,
                        seq(0,1,length=Ngrid)*(60.00-0.00)+0.00
)
)
colnames(X)<-c("X1","X2","X3","X4","X5","X6")
yhat=predict(fit.gp, newdata = X, 
             type = "UK")
# Predict on the new grid, should seem like the 'sausage plot' because we are using krigging here.
# Calculate EI
ego=apply(as.matrix(X),1,EI,fit.gp,type="UK",
          minimization=FALSE)
# Get next x that maximizes the EI within the range of parameter space.
LB<-c(0.127,0.00,0.249,1.27,1.27,0.00);
UB<-c(25.147,181.051,12.575,12.585,6.01,60.00);
INITB<-0.5*c((25.147-0.127),
             (181.051-0.00),
             (12.575-0.249),
             (12.585-1.27),
             (6.01-1.27),
             (60.00-0.00))
#The new observation shall be collected at the location that maximizes the EI(expected imrpovement)
x.new=max_EI(fit.gp,lower=LB,upper=UB, parinit = INITB, minimization=FALSE)$par
#Here we need to set up the lower and upper bound of optimization and its initial value, besides we set minimization=FALSE
#Because we are now caring about f_max instead of f_min in lectures.
yhat.xnew=predict(fit.gp,newdata=x.new,type="UK")$mean



# Update by evaluating our expensive function to yield a new response y.new, in this case is our RMIII!

#In this step we must use RocketSimulator to append n-m runs, which is expensive!

print(paste("X",1:6,"new=",x.new,sep="") )
print(paste("Ynew should be obtained from the RMIII",y.new))
print("++++++++++++++++++")