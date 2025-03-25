
### Examples

#load codes
source("R/Erlang_mixture.R")
Rcpp::sourceCpp("src/rcpp_functions.cpp")

#Example 1
#sample
set.seed(1)
n = 500
pi = c(0.3,0.5,0.2)
true_shape = c(5,20, 30)
true_scale = 0.5

sample <- sample(true_shape,n,replace=TRUE,prob=pi)
x <- rgamma(n,shape = sample, scale = true_scale)

result = est_erlang(x=x,support=NULL,initial=TRUE,setseed=1,comp=length(true_shape),r_range=1:100) #initial=TRUE
print(result)
###result
#shape = (4, 17, 24)
#weight = (0.2740589, 0.5549624, 0.1709787)
#scale = 0.6073742
#log_likelihood = -1407.691




#Example 2
#sample
set.seed(1)
n=500
pi = c(0.15,0.35,0.35,0.15)
true_shape = c(5,20,35,50)
true_scale = 0.5

sample <- sample(true_shape,n,replace=TRUE,prob=pi)
x <- rgamma(n,shape = sample, scale = true_scale)

result = est_erlang(x=x,support=1:5,weight=rep(1/5,5),theta=0.5,initial=FALSE,r_range=1:100) #initial=FALSE
print(result)

###result
#shape = (5, 19, 27, 39, 54)
#weight = (0.1537457, 0.2974126, 0.1590853, 0.3271898, 0.0625667)
#scale = 0.5068856
#log_likelihood = -1653.493