#Bayesian Hierarchical Poisson Model Function:
library(lattice)
library(coda)
hlm <- function(nsamples, burnin, y, group.vector, a_alpha, b_alpha, a_beta, b_beta,
lambda.curr, alpha.curr, beta.curr)
{
# Gibbs Sampling:
"sample.lambda" <- function(nsamples, J, alpha, beta, sigma.y, n)
{
lambda_alpha = sapply(1:J, function(x) alpha + sigma.y[x])
lambda_beta = sapply(1:J, function(x) n[x] + beta)
return(sapply(1:J, function(x) rgamma(n = nsamples, lambda_alpha[x],
lambda_beta[x])))
}
"sample.alpha" <- function(nsamples, a_alpha, b_alpha)
{
u = runif(n=nsamples)
return(qgamma(u, a_alpha, b_alpha))
}
"sample.beta" <- function(nsamples, alpha, beta, J, a_beta, b_beta, sigma.lambda)
{
beta_alpha = J*alpha + a_beta
beta_beta = b_beta + sigma.lambda
return(rgamma(n=nsamples, beta_alpha, beta_beta))
}
# Extract sufficient statistics needed for sampling
N = length(y)
group.levels = unique(group.vector)
J = length(group.levels)
n = rep(NA,J)
sigma.y = rep(NA,J)
sigma.lambda = sum(lambda.curr)
for (i in 1:J){
group.i <- as.logical(group.vector == group.levels[i])
n[i] <- sum(group.i)
sigma.y[i] <- sum(y[group.i])
}
# Storage for the MCMC draws:
post.draws <- matrix(NA,nrow=nsamples,ncol=J+2)
colnames(post.draws) <- c("alpha","beta",paste("lambda_",1:J,sep=""))
# Gibbs sampler:
for (i in 1:(nsamples+burnin))
{
lambda.curr <- sample.lambda(nsamples = 1, J = J, alpha = alpha.curr, beta =
beta.curr, sigma.y = sigma.y, n = n)
alpha.curr <- sample.alpha(nsamples = 1,a_alpha =a_alpha, b_alpha = b_alpha)
beta.curr <- sample.beta(nsamples = 1, alpha = alpha.curr, beta = beta.curr, J =J,
a_beta = a_beta, b_beta = b_beta, sigma.lambda = sigma.lambda)
if (i > burnin){
post.draws[i-burnin,] <- c(lambda.curr, alpha.curr, beta.curr)
}
if (i %% 1000 == 0)
cat(paste("Finished iteration ",i,"...\n",sep=""))
}
# Reformat to MCMC draws:
post.draws <- mcmc(post.draws)
# Return:
return(post.draws)
}
#Bayesian Hierarchical Poisson Sampling:
#Reading in the Data
data = read.table("C:/Users/Ali/Desktop/stars.txt", header=TRUE, quote="\"")
y = data$photon.count
galaxy = data$galaxy
sapply(1:39, function(p) galaxy==p)
#Initial Plot of the Data
plot(y~galaxy, ylab = "Galaxy", main = "Photon Count by Galaxy", col = "blue")
# Prior for Alpha:
a_alpha = 1
b_alpha = 1
# Prior for Beta:
a_beta = 1
b_beta = 1
#Current State:
alpha.curr = 1
beta.curr = 1
lambda.curr = rep(mean(y)^2/var(y), length(unique(galaxy)))
#Number of Samples and Burnin Period
nsamples = 100000
burnin = 10000
lambdaMLE <- sapply(unique(galaxy), function(p) mean(data[c(which(p==galaxy)),2]))
alphabetaratio <- mean(lambdaMLE)
post.draws = hlm(nsamples = nsamples, burnin = burnin, y = y, group.vector = galaxy,
alphabetaratio, b_alpha, a_beta, b_beta, lambda.curr = lambda.curr, alpha.curr =
alphabetaratio, beta.curr = beta.curr)
post.drawsMLE = hlm(nsamples = nsamples, burnin = burnin, y = y, group.vector =
galaxy, 35, b_alpha, a_beta, b_beta, lambda.curr = lambdaMLE, alpha.curr = 35,
beta.curr = beta.curr)
#Sensitivity Analysis
summary(post.draws)
summary(post.drawsMLE)
# plot final 1000 samples:
plot(mcmc(post.draws[(nsamples-10000):nsamples,1:3]))
# Check diagnostics:
# Autocorrelation plots for Alpha, Beta, Lambda 1 and Lambda 2
lag.1.post.acf = apply(post.draws,2,acf,plot=FALSE)
par(mfrow = c(2,2))
lag.1.alpha = unlist(lapply(lag.1.post.acf,function(x){x[1][[1]][[1]][1]}))
lag.1.beta = unlist(lapply(lag.1.post.acf,function(x){x[2][[1]][[1]][1]}))
lag.1.lambda.1 = unlist(lapply(lag.1.post.acf,function(x){x[3][[1]][[1]][1]}))
lag.1.lambda.2 = unlist(lapply(lag.1.post.acf,function(x){x[4][[1]][[1]][1]}))
acf(lag.1.alpha, main = "ACF Plot for Alpha")
acf(lag.1.beta, main = "ACF Plot for Beta")
acf(lag.1.lambda.1, main = "ACF Plot for Lambda 1")
acf(lag.1.lambda.2, main = "ACF Plot for Lambda 2")
# Effective sample sizes:
ess = data.frame(effectiveSize(post.draws))