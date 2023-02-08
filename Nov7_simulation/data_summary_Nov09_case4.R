## Data summary ----
library(e1071)
data_summary <- function(n, mcmc_data, var_rubin, var_wb, var_nonpara, true_mean, null_value){
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
    count_power_rubin[s] <- ifelse(abs(mean_mcmc[s] - null_value) > sqrt(var_rubin[s])*qt(0.975, n-1), 1, 0)
    count_power_wb[s] <- ifelse(abs(mean_mcmc[s] - null_value) > sqrt(var_wb[s])*qt(0.975, n-1), 1, 0)
    count_power_nonpara[s] <- ifelse(abs(mean_mcmc[s] - null_value) > sqrt(var_nonpara[s])*qt(0.975, n-1), 1, 0)
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
  var_true<-var(mean_mcmc)
  
  
  return(list(mean_mcmc = mean_mcmc,
              var_true = var_true,
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

get_all<-function(res){
  out=as.matrix(c(res$mean_mean_mcmc,
                  res$var_true,
                  res$mean_var_rubin_mcmc,
                  res$mean_var_wb_mcmc,
                  res$mean_var_nonpara_mcmc,
                  res$coverage_prob_rubin_wald,
                  res$coverage_prob_wb_wald,
                  res$coverage_prob_nonpara_wald,
                  res$rela_bias_var_rubin,
                  res$rela_bias_var_wb,
                  res$rela_bias_var_nonpara,
                  res$power_rubin,
                  res$power_wb,
                  res$power_nonpara))
  rownames(out)=c("mean","var_true","var_rubin","ar_wb","var_nonpara","coverage_rubin","coverage_wb","coverage_nonpara","RB_rubin",
                  "RB_wb","RB_nonpara","power_rubin","power_wb","power_nonpara")
  return(out)
}
nsim=500
rep_time=500

tmp1<-c4_gen_eqm_50
tmp1_alt <- c4_gen_ineqm_50

tmp1<-c4_gen_eqm_100
tmp1_alt <- c4_dia_ineqm_100

tmp1<-c4_gen_eqm_1000
tmp1_alt <- c4_gen_ineqm_1000

res <- data_summary(100, tmp1$mean_diff_mcmc,tmp1$var_rubin_est[3,], tmp1$var_wb_est[3,], tmp1$var_nonpara_est[3,] , 0, 0)
get_all(res)

res_alt <- data_summary(1000, tmp1_alt$mean_diff_mcmc,tmp1_alt$var_rubin_est[3,], tmp1_alt$var_wb_est[3,], tmp1_alt$var_nonpara_est[3,] , 1, 0)
get_all(res_alt)





tmp1 <- c4_gen_eqm_100
tmp2 <- c4_dia_eqm_100
tmp3 <- c4_comp_eqm_100


par(mfrow=c(1,3))
x1 <- tmp1$mean_trt_mcmc
h<-hist(x1, breaks=20, xlab="value", main="n = 100, Freq ATE, General", ylim = c(0, max(density(tmp1$mle_mean_res[2,])$y)),  freq = FALSE)
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
h<-hist(x4, breaks=20, xlab="value", main="n = 100, Freq ATE, Diagonal", ylim = c(0, 2.1), freq = FALSE)
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
h<-hist(x7, breaks=20, xlab="value", main="n = 100, Freq ATE, Composite", ylim = c(0,4.5), freq = FALSE)
xfit<-seq(min(x7),max(x7),length=200)
yfit<-dnorm(xfit,mean=tmp3$mean_mle[2],sd=tmp3$sd_mle[2])
lines(density(x7), col = "red")
lines(density(tmp3$mle_mean_res[2,]), col = "purple")
lines(xfit, yfit, col="blue", lwd=2)  #naive simulation
legend(x="right",c("mcmc","mle", "mle_normal"),
       lty=c(1,1,1),
       lwd=1,bty="n", 
       col=c("red","purple","blue"))

