require("MASS")
require("nlme")

# covariance matrix of three correlated measurements
var1 <- 1 #variance of measurement 1
var2 <- 1
var3 <- 1
r12 <- 0.8 #correlation between measurement 1 and 2
r13 <- 0.5
r23 <- 0.8
cov12 <-var1*var2*r12
cov13 <-var1*var3*r13
cov23 <-var2*var3*r23
cov21 <-cov12
cov31 <-cov13
cov32 <-cov23
Sigma<-matrix( c(var1,cov12,cov13,
                 cov21,var2,cov23,
                 cov31,cov32,var3),
              3,3)
timetrend <-c(0,2,5) # time trend over measurement 1, 2, 3 in ctl arm
delta <- c(0,2,1)    # effect at measurement 1, 2, 3


#################function for #######################
# data generation for J subjects with K=3 measurements 
# having covariance matrix Sigma
# time trend timetrend
# effects over time delta
simulateData <- function (J,K,timetrend,Sigma,delta)
{
time <- rep(seq(1,K,length=K),J) # K measurements (1,..,K) per person
person <- rep(1:J, each=K)       # J person IDs
# make correlated residuals
a <- mvrnorm(n=J,mu=rep(0,K),Sigma=Sigma) # repeated measurements as rows in a matrix
resid <- c(t(a)) # into a stacked colum of repeated measurements
# fixed effects
trt <- as.numeric( (person %% 2 == 0)*(time > 1) )# every even person id has trt if beyond time 1
fu1 <- as.numeric(time==2) # dummy for change baseline to followup 1 (time=2)
fu2 <- as.numeric(time==3) # dummy for change baseline to followup 2 (time=3)
outcome <- trt*delta + timetrend + resid
return(data.frame(person,time,trt,fu1,fu2,outcome))
}

	# generate data
dat <- simulateData(J=100,K=3,timetrend, Sigma,delta)
	#analyze with compound symmetry correlation structure
glsFit <- gls(outcome ~ fu1 + fu2 + fu1:trt + fu2:trt , 
	data = dat, correlation = corCompSymm(form = ~ 1 | person))
summary(glsFit)
	# unstructured correlation structure, must use form= ~1 | person to use data order 
glsFit <- gls(outcome ~ fu1 + fu2 + fu1:trt + fu2:trt , 
	data = dat, correlation = corSymm(form = ~ 1 | person))
summary(glsFit)


# coefficient and se and p-value
# str(summary(glsFit))
summary(glsFit)$tTable["fu1:trt","Value"]
summary(glsFit)$tTable["fu1:trt","Std.Error"]
summary(glsFit)$tTable["fu1:trt","p-value"]
p1 <- summary(glsFit)$tTable["fu1:trt","p-value"]
signif <- as.numeric(p1<0.05)
signif

###############function to determine power ##############
power_marginalmodel <- function(J,K,timetrend,Sigma,delta,n.sims){
p <- rep(NA,n.sims)
est <- rep(NA, n.sims)
for (s in 1:n.sims){
	fake <- simulateData(J,K,timetrend, Sigma,delta)
	glsFit <- gls(outcome ~ fu1 + fu2 + fu1:trt + fu2:trt , 
			data = fake, correlation = corCompSymm(form = ~ 1 | person))
      est[s] <- summary(glsFit)$tTable["fu1:trt","Value"]
	p[s] <-  summary(glsFit)$tTable["fu1:trt","p-value"]
	}
signif <- as.numeric(p< 0.05)
return(data.frame(est,p,signif))
}

delta <- c(0,1.0,0.5)
res <- power_marginalmodel(J=15,K=3,timetrend,Sigma,delta,n.sims=100)
res
colMeans(res)





