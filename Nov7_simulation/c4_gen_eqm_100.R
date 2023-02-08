## Packages ----
#install.packages("mvtnorm")
#install.packages("reshape2")
library(MASS)
library(dplyr)
library(mvtnorm)
library(reshape2)
library(MCMCpack)

## Functions ----
#' Simulate Longitudinal Data with Monotone Missing from Multivariate Normal Distribution
#'
#' @param n Number of Subjects
#' @param mu Mean of the profile
#' @param sigma Variance-covariance matrix of the profile
#' @param pi Dropout probability for MCAR simulation
#' @param pi_par Dropout probability parameters of logistics model for MAR or MNAR simulation
#' @param trt Treatment group indicator
#' @param missing_type Missing data type, MCAR, MAR or MNAR
#'
#' @export
sim_one_group <- function(n, mu, sigma, pi = NULL, pi_par = NULL, trt = 1, missing_type){
  p <- length(mu)
  n_pattern <- length(pi)
  db <- data.frame(rmvnorm(n = n, mean = mu, sigma = sigma))
  colnames(db) <- paste0("y",1:p)
  db$num <- c(1:n)
  db$id <- paste0(trt, "-", 1:n)
  
  if(missing_type == "MCAR"){
    db$pattern <- sample(1:p, size = n, replace = TRUE, prob = pi)
  }
  
  if(missing_type %in% c("MAR","MNAR")){
    if(pi_par[3] != 0 & missing_type == "MAR"){
      stop("This is not a MAR setting")
    }
    logit_inv <- function(x){
      exp(x) / (1 + exp(x))
    }
    
    
    .pattern <- rep(1, n)
    for(i_missing in p:2){
      .score <- as.matrix(data.frame(1, db[,c(i_missing - 1,i_missing)])) %*% pi_par
      .pi <- logit_inv(as.numeric(.score))
      .pattern <- ifelse( rbinom(n = n, size = 1, prob = .pi) == 1, i_missing, .pattern)
    }
    db$pattern <- .pattern
    
  }
  
  db_comb <- db
  db_comb$trt <- trt
  
  for(i in 1:nrow(db_comb)){
    for(j in 2:p){
      if(db_comb$pattern[i] <= j & db_comb$pattern[i] > 1){
        db_comb[i, j] <- NA
      }
    }
  }
  
  db_long <- melt(db, id.vars = c("id","pattern", "num"),
                  variable.name = c("time") , value.name = "aval")
  db_long <- db_long %>% group_by(id) %>%
    mutate(
      time = as.numeric(time),
      trt = trt) %>%
    ungroup()
  
  return(list(db_comb = db_comb, db_long = db_long))
}

# posterior of mu
#Y=(Yi1....Yit)'; mu is the prior mean on theta, Omega is the prior variance of theta, 
#Sigma is the covariance of Yij, nj is the number of complete Yij in the group 
#each row is a subject
#input=matrix with all subjects' Y; output:matrix including only complete Yi
get_complete<-function(Y_matrix,Y_long){
  T <- length(unique(Y_long$time))
  Y=Y_matrix[,1:T]
  out<-list()
  Ycomp=Y[is.na(Y[,T])==FALSE,]
  Yincomp=Y[is.na(Y[,T])==TRUE,]
  out$complete<-Ycomp
  out$miss<-Yincomp
  return(out)
}


#input:matrix with all subjects' Y; 
#mu:prior mean and Omega: prior covariance of theta;
#Sigma:covariance of each Yi
#output:random sample of theta from posterior distribution
get_theta<-function(Y_matrix,Y_long,mu,Sigma,Omega){
  Ycomp=t(get_complete(Y_matrix,Y_long)$complete)
  nj=ncol(Ycomp)
  Ybar=rowMeans(Ycomp)
  cov=solve(solve(Omega)+nj*solve(Sigma))
  mean=cov%*%(solve(Omega)%*%mu+nj*solve(Sigma)%*%Ybar)
  return(mvrnorm(n = 1,mean,cov))
}


# posterior of Sigma
#return the posterior of Sigma
Sigma_posterior<-function(Y_matrix,Y_long,T_length, theta, v0, S0){
  # T_length is the period length 
  # theta is the mean of the MVN that generates Y   
  
  # nj is the number of the generaed data Y
  # v0 is the first para of Inv-Wishart that generates Y
  # S0 it the INVERSE of the second para of Inv-Wishart that generates Y
  
  #Y_data <- Y_generate_data(mu=mu, Sigma=Sigma, nj=nj, v0=v0, S0=S0)    # Y_data is the complete data,  MATRIX n*p here!
  Y_data=get_complete(Y_matrix,Y_long)$complete
  nj=nrow(Y_data)
  T=ncol(Y_data)
  first_para <- v0 + nj
  S_theta<-matrix(c(rep(0,T*T)), ncol=T)
  tmp_mat <- as.matrix(Y_data)
  for (i in 1: nj){
    S_theta <- S_theta + (tmp_mat[i,]-theta) %*% (t(tmp_mat[i,]-theta))
  }
  second_para <- solve(S0)+S_theta
  inv_sigma_post <- rWishart(1, first_para, solve(second_para))[,,1]
  # Sigma_post_para <- riwish(first_para, second_para)
  Sigma_post_para <- solve(inv_sigma_post)
  
  return(Sigma_post_para)
}

# Ymin given Yobs
# R control covariance
input_mis<- function(x,mu_control,mu_treatment,R){
  library(MASS)
  t<- length(x)
  obs<- x[is.na(x)==FALSE]
  n_missing<- sum(is.na(x)==TRUE)
  mu1<- mu_treatment[is.na(x)==FALSE]
  mu2<- mu_control[is.na(x)==TRUE]
  R11<- R[is.na(x)==FALSE,is.na(x)==FALSE]
  R22<- R[is.na(x)==TRUE,is.na(x)==TRUE]
  R12<- R[is.na(x)==FALSE,is.na(x)==TRUE]
  R21<- R[is.na(x)==TRUE,is.na(x)==FALSE]
  mu<- mu2+R21%*%solve(R11)%*%(obs-mu1)
  sigma<- R22-R21%*%solve(R11)%*%R12
  miss<- mvrnorm(1,mu=mu,Sigma=sigma)
  return(c(obs,miss))
}

mcmc_sim <- function(db_comb, db_long, burn_in_time, rep_time, B){
  db_comb_crtl <- db_comb[which(db_comb$trt == 1),]
  db_comb_trt <- db_comb[which(db_comb$trt == 2),]
  db_long_crtl <- db_long[which(db_long$trt == 1),]
  db_long_trt <- db_long[which(db_long$trt == 2),]
  
  db_imp_ctl <- db_comb_crtl
  db_imp_trt <- db_comb_trt
  
  #get the number of days :T
  T_length=length(unique(db_long_trt$time))
  
  R_crtl <- ifelse(is.na(db_comb_crtl[,T_length]), 0, 1)
  R_trt <- ifelse(is.na(db_comb_trt[,T_length]), 0, 1)
  
  #start imputing
  #Initial value
  # mu = rep(1, T_length)
  mu = rep(0, T_length)
  # Omega=rep(1,T_length)%*%t(rep(1,T_length))
  # Omega = diag(T_length)
  Omega = 100^2*diag(T_length)
  theta=rep(1, T_length)
  Sigma=diag(T_length)
  v0=T_length+1
  S0=v0*diag(T_length)
  
  #initialize
  theta_trt <- get_theta(db_comb_trt, db_long_trt,mu,Sigma,Omega)
  theta_crtl <- get_theta(db_comb_crtl, db_long_crtl,mu,Sigma,Omega)
  
  Sigma_trt <- Sigma_posterior(db_comb_trt, db_long_trt,T_length, theta_trt, v0, S0)
  Sigma_crtl <- Sigma_posterior(db_comb_crtl, db_long_crtl,T_length, theta_crtl, v0, S0)
  
  #estimate pi
  Ytrt_comp=get_complete(db_comb_trt, db_long_trt)$complete
  Ycrtl_comp=get_complete(db_comb_crtl, db_long_crtl)$complete
  pi_hat_trt=nrow(Ytrt_comp)/nrow(db_comb_trt)
  pi_hat_crtl=nrow(Ycrtl_comp)/nrow(db_comb_crtl)
  
  #get the average of Yobs
  Y_trt_avg_T=mean(Ytrt_comp[,T_length])
  Y_crtl_avg_T=mean(Ycrtl_comp[,T_length])
  
  U <- matrix(0, 3, rep_time)
  # print("burn-in")
  # pb <- txtProgressBar(min = 0, max = burn_in_time, style = 3)
  for(i in 1:burn_in_time){
    #get the estimate of theta and Sigma
    theta_trt<-get_theta(db_comb_trt, db_long_trt,mu,Sigma_trt,Omega)
    theta_crtl<-get_theta(db_comb_crtl, db_long_crtl,mu,Sigma_crtl,Omega)
    
    Sigma_trt<-Sigma_posterior(db_comb_trt, db_long_trt,T_length,theta_trt, v0, S0)
    Sigma_crtl<-Sigma_posterior(db_comb_crtl, db_long_crtl,T_length,theta_crtl, v0, S0)
    # setTxtProgressBar(pb, i)
  }
  # close(pb)
  
  vtrt <- c(0)
  vcon <- c(0)
  vdiff <- c(0)
  mis_imp_ctl <- matrix(0, nrow(db_comb_crtl) - nrow(Ycrtl_comp), rep_time)
  mis_imp_trt <- matrix(0, nrow(db_comb_trt) - nrow(Ytrt_comp), rep_time)
  for (k in 1:rep_time){
    #get the estimate of theta and Sigma
    theta_trt<-get_theta(db_comb_trt, db_long_trt,mu,Sigma_trt,Omega)
    theta_crtl<-get_theta(db_comb_crtl, db_long_crtl,mu,Sigma_crtl,Omega)
    
    Sigma_trt<-Sigma_posterior(db_comb_trt, db_long_trt,T_length,theta_trt, v0, S0)
    Sigma_crtl<-Sigma_posterior(db_comb_crtl, db_long_crtl,T_length,theta_crtl, v0, S0)
    
    #get imputed trt group Y
    Ymis=get_complete(db_comb_trt, db_long_trt)$miss
    a<- nrow(Ymis)
    b<- ncol(Ymis)
    t<- matrix(rep(0,a*b),a,b)
    for(l in 1:a){
      t[l,]<- input_mis(Ymis[l,],theta_crtl,theta_trt,Sigma_crtl)
    }
    colnames(t)<- c('y1','y2','y3','y4','y5')
    Y_mistrt_imputed<- t
    comp<- get_complete(db_comb_trt, db_long_trt)$complete
    vtrt[k]<- var(t[,5])/nrow(t)* (1-pi_hat_trt)^2+var(comp[,5])/nrow(comp)*pi_hat_trt^2
    Ytrt_complete=rbind(get_complete(db_comb_trt, db_long_trt)$complete,Y_mistrt_imputed)
    
    #get imputed control group data
    Ymis=get_complete(db_comb_crtl, db_long_crtl)$miss
    a<- nrow(Ymis)
    b<- ncol(Ymis)
    t<- matrix(rep(0,a*b),a,b)
    for(l in 1:a){
      t[l,]<- input_mis(Ymis[l,],theta_crtl,theta_crtl,Sigma_crtl)
    }
    colnames(t)<- c('y1','y2','y3','y4','y5')
    Y_miscrtl_imputed<- t
    comp<- get_complete(db_comb_crtl, db_long_crtl)$complete
    vcon[k]<- var(t[,5])/nrow(t)* (1-pi_hat_crtl)^2+var(comp[,5])/nrow(comp)*pi_hat_crtl^2
    Ycrtl_complete=rbind(get_complete(db_comb_crtl, db_long_crtl)$complete,Y_miscrtl_imputed)
    
    vdiff[k] <- vcon[k] + vtrt[k]
    
    mis_imp_ctl[,k] <- Y_miscrtl_imputed[,T_length]
    mis_imp_trt[,k] <- Y_mistrt_imputed[,T_length]
    
    #estimate u
    Y_mistrt_avg_T=mean(Y_mistrt_imputed[,T_length])
    Y_miscrtl_avg_T=mean(Y_miscrtl_imputed[,T_length])
    u_trt=pi_hat_trt*Y_trt_avg_T + (1-pi_hat_trt)*Y_mistrt_avg_T
    u_crtl=pi_hat_crtl*Y_crtl_avg_T + (1-pi_hat_crtl)*Y_miscrtl_avg_T
    u_diff = u_trt-u_crtl
    U[,k] = c(u_trt,u_crtl,u_diff)
    
    db_imp_ctl[is.na(db_comb_crtl[,T_length]), T_length] <- Y_miscrtl_imputed[,T_length]
    db_imp_trt[is.na(db_comb_trt[,T_length]), T_length] <- Y_mistrt_imputed[,T_length]
    db_imp_total <- c(db_imp_ctl[,T_length], db_imp_trt[,T_length])
    
    
    # save the data
    if(k == 1){
      total_yt_mat <- db_imp_total
    }
    else{
      total_yt_mat <- cbind(total_yt_mat, db_imp_total)
    }
  }
  rownames(total_yt_mat) <- db_comb$id
  
  # (1) Rubin's variance estimate
  W_crtl <- mean(vcon)
  W_trt <- mean(vtrt)
  W_diff <- mean(vdiff)
  var_rubin_est <- c(W_trt, W_crtl, W_diff) + (1 + 1/rep_time)*apply(U, 1, var)
  
  # (2) Using weighted bootstrap to estimate variance
  mu_ctl <- c(0)
  mu_trt <- c(0)
  mu_diff <- c(0)
  for(b in 1:B){
    w_boot_ctl <- rexp(nrow(db_comb_crtl))
    w_boot_trt <- rexp(nrow(db_comb_trt))
    
    pi_ctl <- sum(w_boot_ctl*R_crtl)/sum(w_boot_ctl)
    pi_trt <- sum(w_boot_trt*R_trt)/sum(w_boot_trt)
    
    w_boot_obs_ctl <- w_boot_ctl[which(R_crtl == 1)]
    w_boot_mis_ctl <- w_boot_ctl[which(R_crtl == 0)]
    w_boot_obs_trt <- w_boot_trt[which(R_trt == 1)]
    w_boot_mis_trt <- w_boot_trt[which(R_trt == 0)]
    
    mis_mean_ctl <- sum(w_boot_mis_ctl*rowMeans(mis_imp_ctl))/sum(w_boot_mis_ctl)
    mis_mean_trt <- sum(w_boot_mis_trt*rowMeans(mis_imp_trt))/sum(w_boot_mis_trt)
    
    obs_mean_ctl <- sum(w_boot_obs_ctl*Ycrtl_comp[,T_length])/sum(w_boot_obs_ctl)
    obs_mean_trt <- sum(w_boot_obs_trt*Ytrt_comp[,T_length])/sum(w_boot_obs_trt)
    
    mu_ctl[b] <- obs_mean_ctl*pi_ctl + mis_mean_ctl*(1 - pi_ctl)
    mu_trt[b] <- obs_mean_trt*pi_trt + mis_mean_trt*(1 - pi_trt)
    mu_diff[b] <- mu_trt[b] - mu_ctl[b]
  }
  var_wb_est <- c(var(mu_trt), var(mu_ctl), var(mu_diff))
  
  ## nonparametric bootstrap
  get_mean_fn <- function(yt, db_comb_boot){
    db_comb_ctl <- db_comb_boot[which(db_comb_boot$trt == 1), ]
    db_comb_trt <- db_comb_boot[which(db_comb_boot$trt == 2), ]
    
    mu_ctl <- mean(yt[which(db_comb_boot$trt == 1)])
    mu_trt <- mean(yt[which(db_comb_boot$trt == 2)])
    
    return(c(mu_ctl, mu_trt, mu_trt - mu_ctl))
  }
  
  mean_boot <- matrix(0, 3, B)
  for(b in 1:B){
    id_boot <- sample(unique(db_comb$id), replace = TRUE)
    data_boot <- data.frame(id_new = 1:length(id_boot), id = id_boot)
    db_comb_boot <- merge(data_boot, db_comb,  all.x = TRUE)[,-c(1:2)]
    db_long_boot <- merge(data_boot, db_long,  all.x = TRUE)
    mean_boot[,b] <- rowMeans(apply(total_yt_mat[id_boot,], 2, function(x) get_mean_fn(x, db_comb_boot)))
  }
  var_nonpara_est <- apply(mean_boot, 1, var)
  
  return(list(U = U, 
              var_rubin_est = var_rubin_est, 
              var_wb_est = var_wb_est,
              var_nonpara_est = var_nonpara_est))
}

## MLE
mle_finder <- function(db_comb, db_long, method = "naive"){
  n1 <- length(which(db_comb$trt == 1))
  n2 <- length(which(db_comb$trt == 2))
  p <- length(unique(db_long$time))
  obs_pi <- 1 - c(sum(is.na(db_comb[,p] & db_comb$trt == 1))/n1, 
                  sum(is.na(db_comb[,p] & db_comb$trt == 2))/n2)
  db_comb_ctl <- na.omit(db_comb[which(db_comb$trt == 1),])
  db_comb_trt <- na.omit(db_comb[which(db_comb$trt == 2),])
  
  if(method == "naive"){
    mu_ctl <- mean(db_comb_ctl[,p])
    mu_trt <- mean(db_comb_trt[,p])*obs_pi[2] + mu_ctl*(1 - obs_pi[2])
    mu_diff <- mu_trt - mu_ctl
  }
  est_value <- c(mu_ctl, mu_trt, mu_diff)
  
  return(est_value)
}

## Bootstrap to compute the standard error
sd_bootstrap <- function(db_comb, db_long, B){
  mean_boot_mle <- matrix(0, 3, B)
  for(b in 1:B){
    id_boot <- sample(unique(db_comb$id), replace = TRUE)
    data_boot <- data.frame(id_new = 1:length(id_boot), id = id_boot)
    db_comb_boot <- merge(data_boot, db_comb,  all.x = TRUE)[,-c(1:2)]
    db_long_boot <- merge(data_boot, db_long,  all.x = TRUE)
    mean_boot_mle[,b] <- mle_finder(db_comb_boot, db_long_boot, method = "naive")
  }
  sd_mle <- apply(mean_boot_mle, 1, sd)
  return(sd_mle)
}

## Data ----
# N <- 100
P <- 5
Mu1  = c(0, 1.0, 1.8, 2.5, 3)
Mu2  = c(0, 1.0, 1.8, 2.5, 3)
# Mu2  = c(0, 1.3, 2.3, 3.2, 4)
sd  <- c(2.0, 1.8, 2.0, 2.1, 2.2)
corr   <- matrix(
  c(1, 0.6, 0.3, 0.2, 0.1,
    0.6, 1, 0.7, 0.5, 0.2,
    0.3, 0.7, 1, 0.6, 0.4,
    0.2, 0.5, 0.6, 1, 0.5,
    0.1, 0.2, 0.4, 0.5, 1), 5, 5)

Sigma <- diag(sd) %*% corr %*% diag(sd)


# Missing types
# (1) MCAR
# Pi <- c(60, rep(10,4)) / 100
Pi <- c(80, rep(5, 4))/100

# (2) MAR 
# pi_par <- list(
#   pi_par1 = c(-2.3, 0.2, 0),
#   pi_par2 = c(-2.8, 0.2, 0)
# )

# (3) MNAR
# pi_par <- list(pi_par1 = c(-3.2, 0.2, 0.1), pi_par2 = c(-3.2, 0.2, 0.1))


# ## Convergence check ----
# #trace plot and autocorrelation plot
# par(mfrow=c(3,2))
# x<- U[1,,(burn_in_time+1):rep_time][1,]#sample=50
# acf(x,ylab='sample=50,trt')
# plot(x,type="l",ylab='sample=50,trt')
# x<- U[1,,][2,]#sample=50
# acf(x,ylab='sample=50,control')
# plot(x,type="l",ylab='sample=50,control')
# x<- U[1,,][3,]#sample=50
# acf(x,ylab='sample=50,diff')
# plot(x,type="l",ylab='sample=50,diff')
# 
# x<- U[2,,][1,]#sample=100
# acf(x,ylab='sample=100,trt')
# plot(x,type="l",ylab='sample=100,trt')
# x<- U[2,,][2,]#sample=100
# acf(x,ylab='sample=100,control')
# plot(x,type="l",ylab='sample=100,control')
# x<- U[2,,][3,]#sample=100
# acf(x,ylab='sample=100,diff')
# plot(x,type="l",ylab='sample=100,diff')
# 
# x<- U[3,,][1,]#sample=500
# acf(x,ylab='sample=500,trt')
# plot(x,type="l",ylab='sample=500,trt')
# x<- U[3,,][2,]#sample=500
# acf(x,ylab='sample=500,control')
# plot(x,type="l",ylab='sample=500,control')
# x<- U[3,,][3,]#sample=500
# acf(x,ylab='sample=500,diff')
# plot(x,type="l",ylab='sample=500,diff')

#' Find the mle using only observed data
#' @param n sample size
#' @param nsim Simulation times
#' @param method = "naive": simply taking the average of the observed data at last time point
#' @param missing_type "MCAR", "MAR", "MNAR"
#' @export
mc_sim <- function(N, nsim, burn_in_time, rep_time, B, method, missing_type){
  mle_mean_res <- matrix(0, 3, nsim)
  mle_sd_res <- matrix(0, 3, nsim)
  mean_ctl_mcmc <- c()
  mean_trt_mcmc <- c()
  mean_diff_mcmc <- c()
  var_rubin_est <- matrix(0, 3, nsim)
  var_wb_est <- matrix(0, 3, nsim)
  var_nonpara_est <- matrix(0, 3, nsim)
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for(s in 1:nsim){
    set.seed(s)
    if(missing_type == "MCAR"){
      tmp1 <- sim_one_group(n = N, mu = Mu1, sigma = Sigma, pi = Pi, trt = 1, missing_type = missing_type)
      tmp2 <- sim_one_group(n = N, mu = Mu2, sigma = Sigma, pi = Pi, trt = 2, missing_type = missing_type)
      
      db_comb <- rbind(tmp1[[1]], tmp2[[1]]) # to fit lm
      db_long <- rbind(tmp1[[2]], tmp2[[2]]) # to fit mmrm
    }
    if(missing_type %in% c("MAR", "MNAR")){
      tmp1 <- sim_one_group(n = N, mu = Mu1, sigma = Sigma, pi_par = pi_par$pi_par1, trt = 1, missing_type = missing_type)
      tmp2 <- sim_one_group(n = N, mu = Mu2, sigma = Sigma, pi_par = pi_par$pi_par1, trt = 2, missing_type = missing_type)
      
      db_comb <- rbind(tmp1[[1]], tmp2[[1]]) # to fit lm
      db_long <- rbind(tmp1[[2]], tmp2[[2]]) # to fit mmrm
    }
    
    res_mcmc <- mcmc_sim(db_comb, db_long, burn_in_time, rep_time, B)
    mean_mcmc <- res_mcmc$U
    mean_ctl_mcmc <- c(mean_ctl_mcmc, mean_mcmc[2,])
    mean_trt_mcmc <- c(mean_trt_mcmc, mean_mcmc[1,])
    mean_diff_mcmc <- c(mean_diff_mcmc, mean_mcmc[3,])
    var_rubin_est[,s] <- res_mcmc$var_rubin_est
    var_wb_est[,s] <- res_mcmc$var_wb_est
    var_nonpara_est[,s] <- res_mcmc$var_nonpara_est
    mle_mean_res[,s] <- mle_finder(db_comb, db_long, method = "naive")
    mle_sd_res[,s] <- sd_bootstrap(db_comb, db_long, B)
    
    setTxtProgressBar(pb, s)
  }
  close(pb)
  
  mc_mean_mle <- apply(mle_mean_res, 1, mean)
  mc_sd_mle <- apply(mle_sd_res, 1, mean)
  return(list(mean_ctl_mcmc = mean_ctl_mcmc,
              mean_trt_mcmc = mean_trt_mcmc,
              mean_diff_mcmc = mean_diff_mcmc,
              var_rubin_est = var_rubin_est,
              var_wb_est = var_wb_est,
              var_nonpara_est = var_nonpara_est,
              mean_mle = mc_mean_mle, 
              sd_mle = mc_sd_mle, 
              mle_mean_res = mle_mean_res,
              mle_sd_res = mle_sd_res))
}

## Storing simulation results ----

# different sample size
#nn = c(50,100,1000)
nn <- 100
nsim <- 500
# small nsim <- 100
burn_in_time <- 100
rep_time <- 500
B <- 100

tmp1 <- mc_sim(nn, nsim, burn_in_time, rep_time, B, "naive", "MCAR")

#tmp2 <- mc_sim(nn[2], nsim, burn_in_time, rep_time, B, "naive", "MCAR")

#tmp3 <- mc_sim(nn[3], nsim, burn_in_time, rep_time, B, "naive", "MCAR")

c4_gen_eqm_100 <- tmp1

save(c4_gen_eqm_100, file = "c4_gen_eqm_100.RData")

## Assess the estimator (posterior) ----
# par(mfrow=c(2,2))
# x1 <- tmp1$mean_trt_mcmc
# h<-hist(x1, breaks=20, xlab="value", main="n = 50, trt", ylim = c(0, max(density(tmp1$mle_mean_res[2,])$y)),  freq = FALSE)
# xfit<-seq(min(x1),max(x1),length=200)
# yfit<-dnorm(xfit,mean=tmp1$mean_mle[2],sd=tmp1$sd_mle[2])
# lines(density(x1), col = "red")
# lines(density(tmp1$mle_mean_res[2,]), col = "purple", lwd = 2)
# lines(xfit, yfit, col="blue", lwd=2)  #naive simulation
# # quantile(x1)
# 
# x4 <- tmp2$mean_trt_mcmc
# h<-hist(x4, breaks=20, xlab="value", main="n = 100, trt", ylim = c(0, max(density(tmp2$mle_mean_res[2,])$y)), freq = FALSE)
# xfit<-seq(min(x4),max(x4),length=200)
# yfit<-dnorm(xfit,mean=tmp2$mean_mle[2],sd=tmp2$sd_mle[2])
# lines(density(x4), col = "red")
# lines(density(tmp2$mle_mean_res[2,]), col = "purple")
# lines(xfit, yfit, col="blue", lwd=2)  #naive simulation
# 
# x7 <- tmp3$mean_trt_mcmc
# h<-hist(x7, breaks=20, xlab="value", main="n = 500, trt", ylim = c(0, max(density(tmp3$mle_mean_res[2,])$y)), freq = FALSE)
# xfit<-seq(min(x7),max(x7),length=200)
# yfit<-dnorm(xfit,mean=tmp3$mean_mle[2],sd=tmp3$sd_mle[2])
# lines(density(x7), col = "red")
# lines(density(tmp3$mle_mean_res[2,]), col = "purple")
# lines(xfit, yfit, col="blue", lwd=2)  #naive simulation
# 
# par(mfrow=c(2,2))
# x2 <- tmp1$mean_ctl_mcmc
# xfit<-seq(min(x2),max(x2),length=200)
# yfit<-dnorm(xfit,mean=tmp1$mean_mle[1],sd=tmp1$sd_mle[1])
# h<-hist(x2, breaks=20, xlab="value", main="n = 50, control", ylim = c(0, max(density(tmp1$mle_mean_res[1,])$y)),freq = FALSE)
# lines(density(x2), col = "red")
# lines(density(tmp1$mle_mean_res[1,]), col = "purple", lwd = 2)
# lines(xfit, yfit, col="blue", lwd=2)  #naive simulation
# 
# x5 <- tmp2$mean_ctl_mcmc
# xfit<-seq(min(x5),max(x5),length=200)
# yfit<-dnorm(xfit,mean=tmp2$mean_mle[1],sd=tmp2$sd_mle[1])
# h<-hist(x5, breaks=20, xlab="value", main="n = 100, control",ylim = c(0, max(density(tmp2$mle_mean_res[1,])$y)), freq = FALSE)
# lines(density(x5), col = "red")
# lines(density(tmp2$mle_mean_res[1,]), col = "purple", lwd = 2)
# lines(xfit, yfit, col="blue", lwd=2)  #naive simulation
# 
# x8 <- tmp3$mean_ctl_mcmc
# xfit<-seq(min(x8),max(x8),length=200)
# yfit<-dnorm(xfit,mean=tmp3$mean_mle[1],sd=tmp3$sd_mle[1])
# h<-hist(x8, breaks=20, xlab="value", main="n = 500, control", ylim = c(0, max(density(tmp3$mle_mean_res[1,])$y)), freq = FALSE)
# lines(density(x8), col = "red")
# lines(density(tmp3$mle_mean_res[1,]), col = "purple", lwd = 2)
# lines(xfit, yfit, col="blue", lwd=2)  #naive simulation
# 
# par(mfrow=c(2,2))
# x3 <- tmp1$mean_diff_mcmc
# h<-hist(x3, breaks=20, xlab="value", main="n = 50, ATE", ylim = c(0, max(density(tmp1$mle_mean_res[3,])$y)),freq = FALSE)
# xfit<-seq(min(x3),max(x3),length=200)
# yfit<-dnorm(xfit,mean=tmp1$mean_mle[3],sd=tmp1$sd_mle[3])
# lines(density(x3), col = "red")
# lines(density(tmp1$mle_mean_res[3,]), col = "purple", lwd = 2)
# lines(xfit, yfit, col="blue", lwd=2)  #naive simulation
# 
# x6 <- tmp2$mean_diff_mcmc
# h<-hist(x6, breaks=20, xlab="value", main="n = 100, ATE", ylim = c(0, max(density(tmp2$mle_mean_res[3,])$y)),freq = FALSE)
# xfit<-seq(min(x6),max(x6),length=200)
# yfit <- dnorm(xfit,mean=tmp2$mean_mle[3],sd=tmp2$sd_mle[3])
# lines(density(x6), col = "red")
# lines(density(tmp2$mle_mean_res[3,]), col = "purple", lwd = 2)
# lines(xfit, yfit, col="blue", lwd=2)  #naive simulation
# 
# x9 <- tmp3$mean_diff_mcmc
# h<-hist(x9, breaks=20, xlab="value", main="n = 500, ATE", ylim = c(0, max(density(tmp3$mle_mean_res[3,])$y)),freq = FALSE)
# xfit<-seq(min(x9),max(x9),length=200)
# yfit <- dnorm(xfit,mean=tmp3$mean_mle[3],sd=tmp3$sd_mle[3])
# lines(density(x9), col = "red")
# lines(density(tmp3$mle_mean_res[3,]), col = "purple", lwd = 2)
# lines(xfit, yfit, col="blue", lwd=2)  #naive simulation
# 

## Data summary ----
library(e1071)
data_summary <- function(n, mcmc_data, var_rubin, var_wb, var_nonpara, true_mean){
  mean_mcmc <- c(0)
  var_mcmc_fs <- c(0)
  skew_mcmc <- c(0)
  # CI_quantile <- matrix(0, 2, nsim)
  CI_rubin_wald <- matrix(0, 2, nsim)
  CI_wb_wald <- matrix(0, 2, nsim)
  CI_nonpara_wald <- matrix(0, 2, nsim)
  # count_quantile <- c(0)
  count_CI_rubin_wald <- c(0)
  count_CI_wb_wald <- c(0)
  count_CI_nonpara_wald <- c(0)  
  count_power_rubin <- c(0)
  count_power_wb <- c(0)
  count_power_nonpara <- c(0)
  for(s in 1:nsim){
    mean_mcmc[s] <- mean(mcmc_data[((s-1)*rep_time + 1): (s*rep_time)])
    var_mcmc_fs[s] <- var(mcmc_data[((s-1)*rep_time + 1): (s*rep_time)])
    skew_mcmc[s] <- skewness(mcmc_data[((s-1)*rep_time + 1): (s*rep_time)])
    # CI_quantile[,s] = c(quantile(mcmc_data[((s-1)*rep_time + 1): (s*rep_time)], 0.025),quantile(mcmc_data[((s-1)*rep_time + 1): (s*rep_time)], 0.975))
    # count_quantile[s] <- ifelse(CI_quantile[1,s] < true_mean & CI_quantile[2,s] > true_mean, 1, 0)
    CLB_rubin = mean_mcmc[s]-sqrt(var_rubin[s])*qnorm(.025,lower.tail=FALSE)
    CUB_rubin = mean_mcmc[s]+sqrt(var_rubin[s])*qnorm(.025,lower.tail=FALSE)
    CI_rubin_wald[,s]=c(CLB_rubin,CUB_rubin)
    count_CI_rubin_wald[s] <- ifelse(CI_rubin_wald[1,s] < true_mean & CI_rubin_wald[2,s] > true_mean, 1, 0)
    CLB_wb = mean_mcmc[s]-sqrt(var_wb[s])*qnorm(.025,lower.tail=FALSE)
    CUB_wb = mean_mcmc[s]+sqrt(var_wb[s])*qnorm(.025,lower.tail=FALSE)
    CI_wb_wald[,s]=c(CLB_wb,CUB_wb)
    count_CI_wb_wald[s] <- ifelse(CI_wb_wald[1,s] < true_mean & CI_wb_wald[2,s] > true_mean, 1, 0)
    CLB_nonpara = mean_mcmc[s]-sqrt(var_nonpara[s])*qnorm(.025,lower.tail=FALSE)
    CUB_nonpara = mean_mcmc[s]+sqrt(var_nonpara[s])*qnorm(.025,lower.tail=FALSE)
    CI_nonpara_wald[,s]=c(CLB_nonpara,CUB_nonpara)
    count_CI_nonpara_wald[s] <- ifelse(CI_nonpara_wald[1,s] < true_mean & CI_nonpara_wald[2,s] > true_mean, 1, 0)
    count_power_rubin[s] <- ifelse(abs(mean_mcmc[s] - true_mean) > sqrt(var_rubin[s])*qt(0.975, n-1), 1, 0)
    count_power_wb[s] <- ifelse(abs(mean_mcmc[s] - true_mean) > sqrt(var_wb[s])*qt(0.975, n-1), 1, 0)
    count_power_nonpara[s] <- ifelse(abs(mean_mcmc[s] - true_mean) > sqrt(var_nonpara[s])*qt(0.975, n-1), 1, 0)
  }
  # coverage_prob_quantile <- mean(count_quantile)
  coverage_prob_rubin_wald <- mean(count_CI_rubin_wald)
  coverage_prob_wb_wald <- mean(count_CI_wb_wald)
  coverage_prob_nonpara_wald <- mean(count_CI_nonpara_wald)
  
  # expected_length_quantile <- apply(CI_quantile, 2, function(x) max(x) - min(x))
  expected_length_rubin_wald <- apply(CI_rubin_wald, 2, function(x) max(x) - min(x))
  expected_length_wb_wald <- apply(CI_wb_wald, 2, function(x) max(x) - min(x))
  expected_length_nonpara_wald <- apply(CI_nonpara_wald, 2, function(x) max(x) - min(x))
  
  # power
  power_rubin <- mean(count_power_rubin)
  power_wb <- mean(count_power_wb)
  power_nonpara <- mean(count_power_nonpara)
  
  mean_mean_mcmc <- mean(mean_mcmc)
  mean_var_rubin_mcmc <- mean(var_rubin)
  mean_var_wb_mcmc <- mean(var_wb)
  mean_var_nonpara_mcmc <- mean(var_nonpara)
  # mean_skew_mcmc <- mean(skew_mcmc)
  # mean_expected_length_quantile <- mean(expected_length_quantile)
  mean_expected_length_rubin_wald <- mean(expected_length_rubin_wald)
  mean_expected_length_wb_wald <- mean(expected_length_wb_wald)
  mean_expected_length_nonpara_wald <- mean(expected_length_nonpara_wald)
  
  rela_bias_var_rubin <- (var(mean_mcmc) - mean_var_rubin_mcmc)/var(mean_mcmc)
  rela_bias_var_wb <- (var(mean_mcmc) - mean_var_wb_mcmc)/var(mean_mcmc)
  rela_bias_var_nonpara <- (var(mean_mcmc) - mean_var_nonpara_mcmc)/var(mean_mcmc)
  
 # rela_bias_skew <- (skewness(mean_mcmc) - mean_skew_mcmc)/skewness(mean_mcmc)
  
  
  return(list(mean_mcmc = mean_mcmc,
              var_mcmc_fs = var_mcmc_fs,
              skew_mcmc = skew_mcmc,
              mean_mean_mcmc = mean_mean_mcmc,
              mean_var_rubin_mcmc = mean_var_rubin_mcmc,
              mean_var_wb_mcmc = mean_var_wb_mcmc,
              mean_var_nonpara_mcmc = mean_var_nonpara_mcmc,
              coverage_prob_rubin_wald = coverage_prob_rubin_wald,
              coverage_prob_wb_wald = coverage_prob_wb_wald,
              coverage_prob_nonpara_wald = coverage_prob_nonpara_wald,
              power_rubin = power_rubin,
              power_wb = power_wb,
              power_nonpara = power_nonpara,
              mean_expected_length_rubin_wald = mean_expected_length_rubin_wald,
              mean_expected_length_wb_wald = mean_expected_length_wb_wald,
              mean_expected_length_nonpara_wald = mean_expected_length_nonpara_wald,
              rela_bias_var_rubin = rela_bias_var_rubin,
              rela_bias_var_wb = rela_bias_var_wb,
              rela_bias_var_nonpara = rela_bias_var_nonpara
              ))
}

res <- data_summary(nn[1], tmp1$mean_ctl_mcmc,tmp1$var_rubin_est[1,], tmp1$var_wb_est[1,], tmp1$var_nonpara_est[1,] , 3)
res$coverage_prob_rubin_wald
res$coverage_prob_wb_wald
res$coverage_prob_nonpara_wald
res$rela_bias_var_rubin
res$rela_bias_var_wb
res$rela_bias_var_nonpara
res$power_rubin
res$power_wb
res$power_nonpara

