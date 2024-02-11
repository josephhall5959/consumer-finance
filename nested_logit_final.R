
######PREAMBLE#########
setwd("/Users/josephhall/Dropbox/Projects/Household vs Corporate")
rm(list=ls())

#Load Packages
vec.pac = c("data.table", "dplyr", "ggplot2", "tabulator", "conflicted", "dfadjust",
            "readxl", "fixest", "readr", "doParallel", "stringr", "tidyr", "tictoc", 
            "kableExtra", "nleqslv", "etable", "numDeriv")

for (i in vec.pac){
  if(! i %in% installed.packages()){
    install.packages(i, dependencies = TRUE)
  }
}

lapply(vec.pac, require, character.only = TRUE)

source("Code/nlfunctions.R")

######CODE#######

alpha = -4

#Match the actual market structure 
oo = .92 #Outside option 
nf = 7 #7 total firms 
ts = .47 #Top firm has 47% conditional PQ share 
na = 2 #Two already-accepting firms
as = .17/na #Already accepting firms conditional share 
ns = 1 #Switching firms 
ss = .2/ns #Switching firms share -- reduce aggregate store credit by 4% from 20% = 20%
no = nf - 1 - na - ns #N others 

pq_0 = c(ts*(1-oo), #Large firm 
         rep((1-ts-ss*ns-as*na)*(1-oo)/no, no),  #Other firms
         rep(as*(1-oo), na), #Always accepting firms 
        rep(ss*(1-oo), ns)) #Switching firm 
#Check equal to 1-oo
pq_0
sum(pq_0)

#Consumer Usage 
u_0 = 0.37 
u_1 = u_0 + 0.055

#Pre-Shock Structure
a_0 = c(0,0,0,0,1,1,0)*u_0
g_0 = cbind(c(1,2,3,4, 5, 5,6), 
            c(1,1,1,1,1,1,1))
#Post-Shock Structure
a_1 = c(0,0,0,0,1,1,1)*u_1
g_1 = cbind(c(1,2,3,4, 5, 5,5), 
            c(1,1,1,1,1,1,1))


to_optimize = function(free_params, vector = F){
  
  #Free parameters are sigma and zeta 
  sigma = c(free_params[1],0)
  zeta = free_params[2]
  #oo =free_params[3]
  #Try with Outside Option as a free parameter 
  
  #Costs 
  c_0 = rep(1,7) + a_0* zeta
  c_1 = rep(1,7) + a_1* zeta

  #Solve model pre- and post- as function of params
  #Solve model Pre 
  #What prices and deltas (1) create these PQ shares 
  #and (2) set DPI/DP = 0? 
  #Recover the deltas -- these are structural 
  
  p_delta_0 = c(rep(2,7), c(1,rep(0,6)))
  
  #Function to set equal to zero 
  pq_moments = function(params){
     
    p = params[1:7]
    delta_0 = params[8:14]
    
    delta = delta_0 + alpha * p
    
    s = f.share(delta, sigma, g_0)$sj
    
    s_cond = f.share(delta, sigma, g_0)$sj_hg
    
    pq = s*p
    
    foc = s + (p - c_0)*(alpha/(1-sigma[1]))*s*(1 - sigma[1]*s_cond - (1-sigma[1])*s)
    
    return(c(pq - pq_0, foc))  
  }
  
  p_delta_1 = nleqslv(p_delta_0, pq_moments)
  
  #Check convergence
  p_delta_1
  
  #These are now the structural delta_0s forever  
  delta_0 = p_delta_1$x[8:14]
  p_0 = p_delta_1$x[1:7]
  s_0 = f.share(delta_0 + alpha*p_0, sigma, g_0)$sj
  
  #Recalculate endogenous prices 
  #Function to set equal to zero 
  p_moments = function(params){
    p = params[1:7]
    #delta_0 = params[8:14]
    
    delta = delta_0 + alpha * p
    
    s = f.share(delta, sigma, g_1)$sj
    
    s_cond = f.share(delta, sigma, g_1)$sj_hg
    
    #pq = s*p
    
    foc = s + (p - c_1)*(alpha/(1-sigma[1]))*s*(1 - sigma[1]*s_cond - (1-sigma[1])*s)
    
    return(foc)  
  }
  
  #New eqm prices and quantities 
  p_1 = nleqslv(p_0, p_moments)$x
  s_1 = f.share(delta_0 + alpha*p_1, sigma, g_1)$sj
  
  #Calculate the moments to compare to data
  #Two moments: avg change in sales 
  #And change in sales differential (accept v not)
  #delta_pi = mean((p_1 - c_1)/p_1 - (p_0 - c_0)/p_0)
  delta_t = mean(s_1*p_1/(s_0*p_0) - 1) 
  #Compounded (5 years) to mirror the data underlying the regression this will match 
  #delta_t = (1 + delta_t)^(1/5) - 1
  
  #Differential in Sales (Accept v Not)
  delta_a = mean(s_1[5:6]*p_1[5:6]/(s_0[5:6]*p_0[5:6])) - mean(s_1[-c(5:6)]*p_1[-c(5:6)]/(s_0[-c(5:6)]*p_0[-c(5:6)]))
  #Compounded (5 years) to mirror the data underlying the regression this will match 
  #delta_a = (1 + delta_a)^(1/5) - 1
  
  #Differential in total sales 
  #delta_sls = mean(s_1*p_1/s_0*p_0)
  
  if(vector){return(c(delta_t, delta_a))}
  
  errors = c(delta_t, delta_a) - target_moments
  
  return(errors %*% errors)
}

#Moment 1: delta pi is + 0.009
#Moment 2: delta a is -0.02
#Moment 3: delta sls is -0.07
target_moments = c(.0027, -.048)#, -0.007)

#Try re-running the profit regressions to see if 
#We can get something even smaller 

#Weighting matrix 
#w = matrix(c(1/.005, 0, 0, 1/.0027),2,2)

#Initialise
sigma_0 = 0
zeta_0 = 0 
#oo_0 = 0.1

to_optimize(c(sigma_0, zeta_0), vector = T)

#Optimise
P = optim(c(sigma_0, zeta_0), to_optimize)
#     method = "L-BFGS-B",
#      lower = c(0, -1), upper = c(1,1))$par
P
params = P$par
params
to_optimize(params, vector = T)

#The cost value to report, inclusive of usage f(g) = .37 
params[2]*.37

#Confirm intuitions
#Order is sigma, zeta 
#to_optim_vector = function(x){to_optimize(x,T)}
#to_optim_vector(c(0.01, -0.01))
#to_optim_vector(c(0.02, -0.01))
#Increasing sigma means worse profits and more heterogeneity

sigma = c(params[1],0)
zeta = params[2]
c_0 = rep(1,7) + a_0* zeta

p_delta_0 = c(rep(2,7), c(1,rep(0,6)))

#Function to set equal to zero 
pq_moments = function(params){
  p = params[1:7]
  delta_0 = params[8:14]
  
  delta = delta_0 + alpha * p
  
  s = f.share(delta, sigma, g_0)$sj
  
  s_cond = f.share(delta, sigma, g_0)$sj_hg
  
  pq = s*p
  
  foc = s + (p - c_0)*(alpha/(1-sigma[1]))*s*(1 - sigma[1]*s_cond - (1-sigma[1])*s)
  
  return(c(pq - pq_0, foc))  
}

p_delta_1 = nleqslv(p_delta_0, pq_moments)

delta_0 = p_delta_1$x[8:14]
p_0 = p_delta_1$x[1:7]
s_0 = f.share(delta_0 + alpha*p_0, sigma, g_0)$sj

#Counterfactuals 
#(1) Baseline 
b0 = mean(p_0)
b1 = mean(c_0)
b2 = s_0%*%(p_0 - c_0)
b2b = mean(s_0*(p_0 - c_0))
b3 = 1/exp(1-sum(s_0))
b4 = s_0%*%(p_0 - c_0) + 1/exp(1-sum(s_0))


#(2) Post-Marquette
#Change acceptance and costs
c_1 = rep(1,7) + a_1* zeta

#Calculate new prices 
p_moments = function(params){
  p = params[1:7]
  #delta_0 = params[8:14]
  
  delta = delta_0 + alpha * p
  
  s = f.share(delta, sigma, g_1)$sj
  
  s_cond = f.share(delta, sigma, g_1)$sj_hg
  
  #pq = s*p
  
  foc = s + (p - c_1)*(alpha/(1-sigma[1]))*s*(1 - sigma[1]*s_cond - (1-sigma[1])*s)
  
  return(foc)  
}

#New eqm prices and quantities 
p_1 = nleqslv(p_0, p_moments)$x
s_1 = f.share(delta_0 + alpha*p_1, sigma, g_1)$sj

#Results
mean(p_1)/b0 - 1
mean(c_1)/b1 - 1
#markup
mean(p_1)/mean(c_1) - mean(p_0)/mean(c_0)

s_1%*%(p_1 - c_1)/b2 - 1

mean(s_1*(p_1 - c_1))/b2b - 1

(1/exp(1-sum(s_1)))/b3-1
(s_1%*%(p_1 - c_1) + 1/exp(1-sum(s_1)))/b4-1

#(3) Change Nests But Not Costs
a_1 = c(0,0,0,0,1,1,0)*u_1
g_1 = cbind(c(1,2,3,4, 5,5 ,5), 
            c(1,1,1,1,1,1,1))
c_1 = rep(1,7) + a_0* zeta

#Calculate new prices 
p_moments = function(params){
  p = params[1:7]
  #delta_0 = params[8:14]
  
  delta = delta_0 + alpha * p
  
  s = f.share(delta, sigma, g_1)$sj
  
  s_cond = f.share(delta, sigma, g_1)$sj_hg
  
  #pq = s*p
  
  foc = s + (p - c_1)*(alpha/(1-sigma[1]))*s*(1 - sigma[1]*s_cond - (1-sigma[1])*s)
  
  return(foc)  
}

#New eqm prices and quantities 
p_1 = nleqslv(p_0, p_moments)$x
s_1 = f.share(delta_0 + alpha*p_1, sigma, g_1)$sj

#Results
mean(p_1)/b0 - 1
mean(c_1)/b1 - 1
s_1%*%(p_1 - c_1)/b2 - 1
(1/exp(1-sum(s_1)))/b3-1
(s_1%*%(p_1 - c_1) + 1/exp(1-sum(s_1)))/b4-1

#(4) Change costs but not nests 
a_1 = c(0,0,0,0,1,1,1)*u_1
g_1 = cbind(c(1,2,3,4, 5, 5,6), 
            c(1,1,1,1,1,1,1))
c_1 = rep(1,7) + a_1* zeta

#Calculate new prices 
p_moments = function(params){
  p = params[1:7]
  #delta_0 = params[8:14]
  
  delta = delta_0 + alpha * p
  
  s = f.share(delta, sigma, g_1)$sj
  
  s_cond = f.share(delta, sigma, g_1)$sj_hg
  
  #pq = s*p
  
  foc = s + (p - c_1)*(alpha/(1-sigma[1]))*s*(1 - sigma[1]*s_cond - (1-sigma[1])*s)
  
  return(foc)  
}

#New eqm prices and quantities 
p_1 = nleqslv(p_0, p_moments)$x
s_1 = f.share(delta_0 + alpha*p_1, sigma, g_1)$sj

#Results
mean(p_1)/b0 - 1
mean(c_1)/b1 - 1
s_1%*%(p_1 - c_1)/b2 - 1
(1/exp(1-sum(s_1)))/b3-1
(s_1%*%(p_1 - c_1) + 1/exp(1-sum(s_1)))/b4-1


#(5) Set Acceptance to Zero
a_1 = c(0,0,0,0,0,0,0)
g_1 = cbind(c(1,2,3,4, 5, 6,7), 
            c(1,1,1,1,1,1,1))
c_1 = rep(1,7) + a_1* zeta

#Calculate new prices 
p_moments = function(params){
  p = params[1:7]
  #delta_0 = params[8:14]
  
  delta = delta_0 + alpha * p
  
  s = f.share(delta, sigma, g_1)$sj
  
  s_cond = f.share(delta, sigma, g_1)$sj_hg
  
  #pq = s*p
  
  foc = s + (p - c_1)*(alpha/(1-sigma[1]))*s*(1 - sigma[1]*s_cond - (1-sigma[1])*s)
  
  return(foc)  
}

#New eqm prices and quantities 
p_1 = nleqslv(p_0, p_moments)$x
s_1 = f.share(delta_0 + alpha*p_1, sigma, g_1)$sj

#Results
mean(p_1)/b0 - 1
mean(c_1)/b1 - 1
s_1%*%(p_1 - c_1)/b2 - 1
(1/exp(1-sum(s_1)))/b3-1
(s_1%*%(p_1 - c_1) + 1/exp(1-sum(s_1)))/b4-1

#(6) Set Acceptance and Usage to 1 
#"Homogeneous payment network"
a_1 = c(1,1,1,1,1,1,1)
g_1 = cbind(c(1,1,1,1,1,1,1), 
            c(1,1,1,1,1,1,1))
c_1 = rep(1,7) + a_1* zeta

#Calculate new prices 
p_moments = function(params){
  p = params[1:7]
  #delta_0 = params[8:14]
  
  delta = delta_0 + alpha * p
  
  s = f.share(delta, sigma, g_1)$sj
  
  s_cond = f.share(delta, sigma, g_1)$sj_hg
  
  #pq = s*p
  
  foc = s + (p - c_1)*(alpha/(1-sigma[1]))*s*(1 - sigma[1]*s_cond - (1-sigma[1])*s)
  
  return(foc)  
}

#New eqm prices and quantities 
p_1 = nleqslv(p_0, p_moments)$x
s_1 = f.share(delta_0 + alpha*p_1, sigma, g_1)$sj

#Results
mean(p_1)/b0 - 1
mean(c_1)/b1 - 1
s_1%*%(p_1 - c_1)/b2 - 1
(1/exp(1-sum(s_1)))/b3-1
(s_1%*%(p_1 - c_1) + 1/exp(1-sum(s_1)))/b4-1



