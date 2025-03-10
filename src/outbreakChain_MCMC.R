
##############################################################################
##############################################################################
## Preliminary Model code for chapter 2 project                             ##
##############################################################################
##############################################################################

####################
## Load necessary ##
##    packages    ##
####################

library(ggplot2)
library(viridis)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(colorspace)
library(igraph)
library(ggraph)
####################

## set working directory to project root folder. In this case './HostJump_model'
setwd("~/Desktop/Repos/Nipah_localPhi/") ## ensure this path matches the local path to the project root folder

#################
#################
## Load  model ##
##  functions  ##
#################
#################

## C++ simulation wrapper to generate transmission trees
#### or use data on spillover chain length to simulate
#### possible transmission trees of a given size
transmissionTree_sim <- function(n_reps, parameters){
  
  parvec = c( format(n_reps, scientific = F),
              parameters$model_type, # 0, 
              parameters$R0, # 1, 
              parameters$k, # 2, gamma shape parameter
              parameters$R0_SS, # 3, 
              parameters$prob_SS, # 4, 
              
              parameters$data_type, # 5, 
              parameters$reps, # 6, 
              format(parameters$N_threshold, scientific = F), # 7, 
              parameters$data_path # 9, 
  )
  
  strvec = format(parvec, digits = 5)
  
  setwd("~/Desktop/Repos/Nipah_localPhi/src") ## call has to be from location of .exe file or parameter read-in fails???
  
  ## Run the model
  ## The path to the bTB cpp binary file must be set correctly in the sys call below:
  nm = paste0("./transmissionTree_model.exe")
  r <- system2( nm, args = strvec, stdout = TRUE)
  
  ## read output from console
  out <- read.table(text = r, header = TRUE, sep = ';', check.names = FALSE) %>% mutate_all(as.numeric)
  
  return( out ) 
}


## MCMC framework for estimating branching tree model parameters
#### options for Poisson and Geometric offspring distributions (1 parameter),
#### negative binomial (2 parameters), or Poisson-mixture (3 parameters)
#### Calculates WAIC, DIC for assessing different numbers of parameters
mcmc.branching <- function(y, # observed/simulated data
                           model.parameters ## list of model parameters, including the type of branching model to fit
                           
){
  ## prep for MCMC
  accept = 0
  accept_tot = 0
  theta.save = matrix(NA, nrow = n.mcmc, ncol = length(pars.start)) %>% as.data.frame()
  colnames(theta.save) <- pars.start.names
  D_bar_i = rep(NA, n.mcmc) 
  w_it = matrix(NA, nrow = n.mcmc, ncol = length(y[,1])) # q_it to be converted into weights to preserve memory
  LL = rep(NA, n.mcmc+1) 
  adapt.t = 0 # time to next adaptive step
  
  # calculate likelihood of initial models parameters
  out <- do.call( model.function, list(model.parameters) ) 
  current.ll = do.call( disease_model_LL, list(observed_data = y, model_preds = out) )
  LL[1] = current.ll
  
  
  ## start MCMC
  for(iter in 1:n.mcmc){
  
    # print out every 100 iterations to track progress
    if(iter%%100 == 0){
      cat(iter," ")
    }
    
    
    
    # evaluate model likelihood with proposed values
    new.ll = do.call( disease_model_LL, list(observed_data = y, model_preds = out) )
    
    # print parameter values if new likelihood is undefined
    if(is.nan(new.ll)){new.ll = -Inf; print(paste0('error in new likelihood: ')); print(bio_transform(x = theta.star, transform = pars.transform)); }
    
    
    ## MH step
    mh1 = new.ll + do.call( ll.prior, list(pars.val = theta.star, hyper_pars = ll.prior.hyper) )
    mh2 = current.ll + do.call( ll.prior, list(pars.val = pars.start, hyper_pars = ll.prior.hyper))
    u = runif(1)
    
    # accept/reject accounting for proposal symmetry (or asymetry)
    if( u < exp(mh1 - mh2)*proposal.weight ){
      accept = accept + 1
      accept_tot = accept_tot + 1
      pars.start = bio_transform(x = theta.star, transform = pars.transform) # convert accepted values to biological scale
      current.ll <- new.ll
    }else{
      print(paste0(exp(mh1 - mh2)*proposal.weight, ' < ', u ))
      pars.start = bio_transform(x = pars.start, transform = pars.transform) # convert accepted values to biological scale
    }
    
    ## save values
    theta.save[iter,] = pars.start # pars.start should now be on the biologically realistic scale from if-else above
    
    LL[iter+1] = current.ll
    # calculate q_its for AIC
    w_it[iter,] = 1/current.ll
    
    # calculate DIC parameter component
    # log likelyhood of data given the ith iteration of MCMC parameter values
    D_bar_i[iter] = sum(current.ll)
    
    ## Adapt,
    if(n_adapt > 0){
      if( iter%%n_adapt == 0 ){
        ## move adapt counter up 1
        adapt.t = adapt.t+1
        ## new tuning parameters for beta
        adapt.vals = get.sigma(s2.tune = s2.tune, Sigma.tune = Sigma.tune, data = theta.save[(iter-n_adapt)+1:n_adapt,], accept = accept, t.adapt = adapt.t)
        Sigma.tune = adapt.vals$Sigma
        s2.tune = adapt.vals$s2
        ## resetting acceptances to 0
        accept = 0
      }
    }
  }# end MCMC
}

## code to implement log-adaptive tuning following Shaby and Wells
#### s2.tune=current variance scaling
#### Sigma.tune=current Covariance matrix of samples
#### data = matrix or vector of the most recent "k" MCMC iterations
#### accept = number of MH acceptances in the most recent "k" iterations
#### t.adapt = the number of times, including this one, that we have
####           applied the adaptive tuning procedure.
get.sigma = function(s2.tune, Sigma.tune, data, accept, t.adapt, c0 = 1, c1 = 0.8){
  
  
  if(is.matrix(Sigma.tune)){
    ## case for multiple params
    k = nrow(data) ## number of mcmc iterations
  }else{
    k = length(data)
  }
  r.hat = accept/k ## empirical acceptance prob
  S.hat = var(data)
  gamma1 = 1/t.adapt^c1
  gamma2 = c0*gamma1
  s2.new = exp( log(s2.tune) + gamma2*(r.hat-.234) )
  Sigma.new = Sigma.tune + gamma1*(S.hat-Sigma.tune)
  list(s2 = s2.new, Sigma = Sigma.new)
}


#################
#################



############################
############################
## Test block with simple ##
## hard-coded MCMC for    ##
## NB offspring model     ##
############################
############################

pars.true = data.frame(R0 = 1.1, k = 4)
model.pars <- data.frame(model_type = 0, R0 = pars.true$R0, R0_SS = 10, prob_SS = 0.05, k = pars.true$k,
                         data_type = 0, reps = 10, N_threshold = 50, 
                         data_path = "na"
)


library(dplyr)
library(coda)

N <- 20  # Number of spillover events

## generate outbreak data for each spillover
tst <- lapply(1:N, function(i) {
  # Generate transmission tree
  tmp <- transmissionTree_sim(n_reps = 1, parameters = model.pars)
  # print(sum(tmp$to == -1))
  tmp = tmp[!(tmp$to == -1),] ## remove terminal nodes with no offspring distribution due to simulation threshold
  
  # Calculate out-degree: Count how many times each node appears in "from"
  out_degree = tmp %>% count(from)
  out_degree$n[out_degree$from == 1] = out_degree$n[out_degree$from == 1] - 1 ## reduce focal node magnitude by one, since it contains a self reference for visualization purposes
  
  out_data = data.frame(node = 1:max(tmp$from), out_degree = 0)
  out_data$out_degree[unique(tmp$from)] = out_degree$n
  out_data$tree_ID <- i
  
  return(out_data)
  
}) %>% bind_rows()

dim(tst)[1]



n.mcmc = 50000
n_adapt = 500
pars.fit = data.frame(R0 = 1, k = 1) ## initial parameter values
model.pars[names(pars.fit)] = pars.fit
pars.save = matrix(NA, ncol = dim(pars.fit)[2], nrow = n.mcmc)
update.pars = 0

sample.size = 100
y = tst$out_degree
y = sample(x = tst$out_degree, size = min(length(tst$out_degree), sample.size), replace = F)
length(y)
# y = rnbinom(n=50, mu = pars.true$R0, size = pars.true$k)

get.p = function(R0, k){(1 + (R0/k))^(-1)} ## converts true model parameters into NB prob parameter


## R0~gamma; k~lim*Beta
## prior hyperparameters
k_shape = 1; k_scale = 0.5
R0_var = 0.1;

for(iter in 1:n.mcmc){
  
  # print out every 100 iterations to track progress
  if(iter%%100 == 0){
    cat(iter," ")
  }
  
  ## initial parameters
  k = as.numeric(pars.fit['k'])
  R0 = as.numeric(pars.fit['R0'])
  p = get.p(R0 = R0, k = k)
  
  adapt.t = 0
  
  ## improvise
  ## propose MH value for R0, k
  # k.star = rgamma(1, shape = k/k_var, scale = k_var) ## gamma proposal
  
  
  ## v = (k*(1-k)/k_var)-1; a = k*v; b = (1-k)*v ## compute beta parameters using mean-variance parameterization
  # k.star = rbeta(1, shape1 = a, shape2 = b) ## stretched beta proposal
  # k.star = k
  # v.star = (k.star*(1-k.star)/k_var)-1; a.star = (k.star)*v; b.star = (1-k.star)*v ## compute proposal beta parameters using mean-variance parameterization
  k.star = rgamma(1, shape = k_shape, scale = k_scale)
  
  R0.star = rgamma(1, shape = R0/R0_var, scale = R0_var)
  p.star = get.p(R0 = R0.star, k = k.star) # convert sampled parameters into NB parameter p
  
  
  ## gamma proposal MH ratio terms
  mh.num = sum( dnbinom(y, mu = R0.star, size = k.star, log = TRUE) )  +
    dgamma(k, shape = k_shape, scale = k_scale, log = TRUE) +
    dgamma(R0, shape = R0.star/R0_var, scale = R0_var, log = TRUE)

  mh.denom = sum( dnbinom(y, mu = R0, size = k, log = TRUE) ) +
    dgamma(k.star, shape = k_shape, scale = k_scale, log = TRUE) +
    dgamma(R0.star, shape = R0/R0_var, scale = R0_var, log = TRUE)
  
  ## model assuming homogeneous transmission
  mh.num = sum( dpois(y, lambda = R0.star, log = TRUE) )  +
    dgamma(R0, shape = R0.star/R0_var, scale = R0_var, log = TRUE)
  
  mh.denom = sum( dpois(y, lambda = R0, log = TRUE) ) +
    dgamma(R0.star, shape = R0/R0_var, scale = R0_var, log = TRUE)
  
  ## gamma proposal MH ratio terms
  # mh.num = sum( dnbinom(y, size = k.star, prob = p.star, log = TRUE) )  + 
  #   dgamma(k, shape = k.star/k_var, scale = k_var, log = TRUE) + dgamma(R0, shape = R0.star/R0_var, scale = R0_var, log = TRUE)
  # 
  # mh.denom = sum( dnbinom(y, size = k, prob = p, log = TRUE) ) +
  #               dgamma(k.star, shape = k/k_var, scale = k_var, log = TRUE) + dgamma(R0.star, shape = R0/R0_var, scale = R0_var, log = TRUE)
  # 
  
  
  
  
  if(is.na(exp(mh.num-mh.denom))){
    print(paste0("mh.num: ", mh.num))
    print(paste0("mh.denom: ", mh.denom))
    print(paste0("k: ", k.star))
    print(paste0("R0: ", R0.star))
    print(paste0("p: ", p.star))
  }
  
  if(runif(1) < exp(mh.num-mh.denom)){
    pars.fit$R0 = R0.star 
    pars.fit$k = k.star 
    update.pars = update.pars+1 ## count updated parameters
  }
  
  ## save out parameters
  pars.save[iter,1] = pars.fit$R0
  pars.save[iter,2] = pars.fit$k*lim
  
  
  ## adapt
  # if(iter%%n_adapt == 0){
  #   ## move adapt counter up 1
  #   adapt.t=adapt.t+1
  #   ## new tuning parameters for beta
  #   adapt.vals=get.sigma(s2.tune = s2.tune, Sigma.tune = S2, data = pars.save[(iter-n_adapt)+1:n_adapt,], accept = update.pars, t.adapt = adapt.t)
  #   S2=adapt.vals$Sigma
  #   s2=adapt.vals$s2
  #   ## resetting acceptances to 0
  #   update.pars=0
  # }
  
} ## overcome

update.pars/n.mcmc
  
effectiveSize(pars.save)


# matplot(pars.save,type="l")
burnin = 1000
pars.save = as.data.frame(pars.save)
pars.save = pars.save[-(1:burnin),]
pars.save_thin = pars.save[seq(0,dim(pars.save)[1], 500),]
names(pars.save) <- names(pars.true)
as.data.frame(lapply(pars.save, mean))
pars.true
abs(as.data.frame(lapply(pars.save, mean)) - pars.true)
dim(tst)[1]

############################
############################

# parameter validation plots
plot(pars.save$R0, type="l")
abline(h = mean(pars.save$R0), col = "red", lwd = 3)
abline(h = pars.true$R0, col = "gold", lwd = 3)
abline(h = quantile(pars.save$R0, c(.025, .975)), col = "aquamarine", lty = "dashed", lwd = 2)

plot(pars.save$k, type="l")
abline(h = mean(pars.save$k), col = "red", lwd = 3)
abline(h = pars.true$k, col = "gold", lwd = 3)
abline(h = quantile(pars.save$k, c(.025, .975)), col = "aquamarine", lty = "dashed", lwd = 2)
ecdf(pars.save$k)(pars.true$k)


abs( ecdf(pars.save$R0)(pars.true$R0) - 0.5 )
rate_est = 1 / (mean(pars.save$R0[pars.save$R0>1]) - 1)
df = as.data.frame(pars.save$R0[pars.save$R0>1])
names(df) = "R0"
aprox_dens = data.frame(R0 = seq(1,round(max(df$R0), digits = 1),0.001))
aprox_dens$density = dexp(x = aprox_dens$R0 - 1, rate = rate_est)
ggplot(data = df) +
  geom_histogram(aes( x = R0, y = after_stat(density)), colour = 1, fill = "white", bins = 20) +
  geom_line(data = aprox_dens, aes(x = R0, y = density))

df2 = as.data.frame(pars.save$R0)
names(df2) = "R0"
ggplot(data = df2) +
  geom_histogram(aes( x = R0, y = after_stat(density)), colour = 1, fill = "white")


df3 = data.frame( phi = 1 - (1 / (1 + rexp(10000, rate = rate_est))) )
ggplot(data = df3) +
  geom_histogram(aes( x = phi, y = after_stat(density)), colour = 1, fill = "white", bins = 20)
#################
#################


