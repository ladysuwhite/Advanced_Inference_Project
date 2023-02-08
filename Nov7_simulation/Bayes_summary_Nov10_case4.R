##bayes estimator summary

## Data summary ----
library(e1071)
data_summary <- function(n, mcmc_data, var_rubin, var_bayes, true_mean, null_value){
  mean_mcmc <- c(0)
  var_mcmc_fs <- c(0)
  skew_mcmc <- c(0)
  CI_rubin_wald <- matrix(0, 2, nsim)
  CI_bayes_wald <- matrix(0, 2, nsim)
  count_CI_rubin_wald <- c(0)
  count_CI_bayes_wald <- c(0)  
  count_power_rubin <- c(0)
  count_power_bayes <- c(0)
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
    CLB_bayes = mean_mcmc[s]-sqrt(var_bayes[s])*qnorm(.025,lower.tail=FALSE)
    CUB_bayes = mean_mcmc[s]+sqrt(var_bayes[s])*qnorm(.025,lower.tail=FALSE)
    CI_bayes_wald[,s]=c(CLB_bayes,CUB_bayes)
    count_CI_bayes_wald[s] <- ifelse(CI_bayes_wald[1,s] < true_mean & CI_bayes_wald[2,s] > true_mean, 1, 0)
    count_power_rubin[s] <- ifelse(abs(mean_mcmc[s] - null_value) > sqrt(var_rubin[s])*qt(0.975, n-1), 1, 0)
    count_power_bayes[s] <- ifelse(abs(mean_mcmc[s] - null_value) > sqrt(var_bayes[s])*qt(0.975, n-1), 1, 0)
  }
  # coverage_prob_quantile <- mean(count_quantile)
  coverage_prob_rubin_wald <- mean(count_CI_rubin_wald)
  coverage_prob_bayes_wald <- mean(count_CI_bayes_wald)
  
  # expected_length_quantile <- apply(CI_quantile, 2, function(x) max(x) - min(x))
  expected_length_rubin_wald <- apply(CI_rubin_wald, 2, function(x) max(x) - min(x))
  expected_length_bayes_wald <- apply(CI_bayes_wald, 2, function(x) max(x) - min(x))
  
  # power
  power_rubin <- mean(count_power_rubin)
  power_bayes <- mean(count_power_bayes)
  
  mean_mean_mcmc <- mean(mean_mcmc)
  mean_var_rubin_mcmc <- mean(var_rubin)
  mean_var_bayes_mcmc <- mean(var_bayes)
  mean_expected_length_rubin_wald <- mean(expected_length_rubin_wald)
  mean_expected_length_bayes_wald <- mean(expected_length_bayes_wald)
  
  rela_bias_var_rubin <- (var(mean_mcmc) - mean_var_rubin_mcmc)/var(mean_mcmc)
  rela_bias_var_bayes <- (var(mean_mcmc) - mean_var_bayes_mcmc)/var(mean_mcmc)
  
  # true variance
  true_var_mcmc <- var(mean_mcmc)
  
  return(list(mean_mcmc = mean_mcmc,
              var_mcmc_fs = var_mcmc_fs,
              skew_mcmc = skew_mcmc,
              mean_mean_mcmc = mean_mean_mcmc,
              mean_var_rubin_mcmc = mean_var_rubin_mcmc,
              mean_var_bayes_mcmc = mean_var_bayes_mcmc,
              coverage_prob_rubin_wald = coverage_prob_rubin_wald,
              coverage_prob_bayes_wald = coverage_prob_bayes_wald,
              power_rubin = power_rubin,
              power_bayes = power_bayes,
              mean_expected_length_rubin_wald = mean_expected_length_rubin_wald,
              mean_expected_length_bayes_wald = mean_expected_length_bayes_wald,
              rela_bias_var_rubin = rela_bias_var_rubin,
              rela_bias_var_bayes = rela_bias_var_bayes, 
              true_var_mcmc = true_var_mcmc
  ))
}

get_all<-function(res){
  out=as.matrix(c(res$mean_mean_mcmc,
                  res$true_var_mcmc,
                  res$mean_var_bayes_mcmc,
                  res$rela_bias_var_bayes,
                  res$coverage_prob_bayes_wald,
                  res$power_bayes))
  rownames(out)=c("mean","var_true","var_bayes","RB_Bayes","Coverage_Bayes","Power_Bayes")
  return(out)
}

nsim=500
rep_time=200


T3<-tmp3
res3 <- data_summary(100, T3$mean_diff_mcmc,T3$var_rubin_est[3,], T3$var_bayes_est[3,],0,0)
get_all(res3)

tmp1<-tmp2
tmp2<-bc4_dia_eqm_100
tmp3<-bc4_comp_eqm_100


par(mfrow=c(1,3))
x1 <- tmp1$mean_trt_mcmc
h<-hist(x1, breaks=20, xlab="value", main="n = 100, Bayesian ATE, General", ylim = c(0, max(density(tmp1$mle_mean_res[2,])$y)),  freq = FALSE)
xfit<-seq(min(x1),max(x1),length=200)
yfit<-dnorm(xfit,mean=tmp1$mean_mle[2],sd=tmp1$sd_mle[2])
lines(density(x1), col = "red")
lines(density(tmp1$mle_mean_res[2,]), col = "purple", lwd = 2)
lines(xfit, yfit, col="blue", lwd=2)  #naive simulation
legend(x="right",c("mcmc","mle", "mle_normal"),
       lty=c(1,1,1),
       lwd=1,bty="n", 
       col=c("red","purple","blue"))


x4 <- tmp2$mean_trt_mcmc
h<-hist(x4, breaks=20, xlab="value", main="n = 100, Bayesian ATE, Diagonal", ylim = c(0, 1.52), freq = FALSE)
xfit<-seq(min(x4),max(x4),length=200)
yfit<-dnorm(xfit,mean=tmp2$mean_mle[2],sd=tmp2$sd_mle[2])
lines(density(x4), col = "red")
lines(density(tmp2$mle_mean_res[2,]), col = "purple")
lines(xfit, yfit, col="blue", lwd=2)  #naive simulation
legend(x="right",c("mcmc","mle", "mle_normal"),
       lty=c(1,1,1),
       lwd=1,bty="n", 
       col=c("red","purple","blue"))



x7 <- tmp3$mean_trt_mcmc
h<-hist(x7, breaks=20, xlab="value", main="n = 100, Bayesian ATE, Composite", ylim = c(0,5), freq = FALSE)
xfit<-seq(min(x7),max(x7),length=200)
yfit<-dnorm(xfit,mean=tmp3$mean_mle[2],sd=tmp3$sd_mle[2])
lines(density(x7), col = "red")
lines(density(tmp3$mle_mean_res[2,]), col = "purple")
lines(xfit, yfit, col="blue", lwd=2)  #naive simulation
legend(x="right",c("mcmc","mle", "mle_normal"),
       lty=c(1,1,1),
       lwd=1,bty="n", 
       col=c("red","purple","blue"))



