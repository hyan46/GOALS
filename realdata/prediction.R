rm(list=ls())
library(doParallel)
library(foreach)
library(quantreg)
library(xtable)
library(latex2exp)
registerDoParallel(6)
source("../Util_pool.R")
source("../Util_SN.R")
source("../Util_BIC.R")
source("./Util_noncrossing.R")
# source("./Util_pred.R")
Sys.setlocale("LC_TIME", "English")

######## Basic setting
code <- 'US'

######## multiple grid-trim
parameter_list <- cbind(c(0.06,rep(c(0.08,0.1,0.12,0.15),each=2)), c(0.01,rep(c(0.01,0.02),4)))
critical_table <- read.csv("../SN_MultipleSlope_CriticalValue.csv")
tau_sets <- c(0.1,0.5,0.9)
num_tau <- length(tau_sets)

######## Data
full_data <- read.csv("./owid-covid-data.csv")
country_name <- "United States"
country_data <- subset(full_data, (location==country_name)&(total_cases>1000)&(new_cases>=1))
date <- format(as.Date(country_data$date),format="%Y-%b-%d")
y <- log(country_data$new_cases)
y[is.na(y)] <- 0
n <- length(y)
last_forecast_point <- '2020-Nov-09'
last_observe <- '2020-Nov-21'
last_point <- which(date==last_observe)
print(last_point)
y <- y[1:last_point]
reference_point <- which(date==last_forecast_point)

######## Prediction in-sample fit
num_weeks_back <- (reference_point-which(date=='2020-Aug-03'))/7 # backtest how many weeks (Aug-03 CDC data available)
prediction_result <- foreach(pred_grid=0:num_weeks_back, .packages='quantreg') %dopar% {
  start_point <- 1
  end_point <- reference_point-pred_grid*7
  datasub <- y[start_point:end_point]
  ### Sensitivity analysis
  result <- list()
  for (par in 1:nrow(parameter_list)){
    grid <- parameter_list[par,1]
    trim <- parameter_list[par,2]
    critical_sets <- subset(critical_table, epsilon==grid&delta==trim)
    critical_sets <- unlist(critical_sets)[3:5]
    grid_size <- floor(length(datasub)*grid) # round grid size
    trim_size <- ceiling(length(datasub)*trim) # ceiling trim size
    SN_sweep_result <- tryCatch({SN_sweep_quantile_pool(data=datasub, grid_size=grid_size, trim_size=trim_size, tau_sets=tau_sets)},
                                error=function(cond){return(NULL)})
    if(is.null(SN_sweep_result)){
      SN_result90 <- SN_result95 <- NULL
    }else{
      max_SN <- max_SNsweep(SN_sweep_result)
      SN_result90 <- SN_local_peak(grid_size, max_SN, critical_value=critical_sets[1])
    }
    result[[par]] <- list(grid=grid, trim=trim, SN_result90=SN_result90, SN_sweep_result=SN_sweep_result)
  }
  list(result=result)
}

######## SNCP selected by BIC
prediction_result_BIC <- foreach(pred_grid=0:num_weeks_back, .packages='quantreg') %dopar% {
  start_point <- 1
  end_point <- reference_point-pred_grid*7
  datasub <- y[start_point:end_point]
  n <- length(datasub)
  # BIC selection
  BIC_factor <- 0.1
  sum_90 <-  list()
  BIC90_individual_model_mtau <-  c()

  for (par in 1:nrow(parameter_list)){
    SN_result90 <- prediction_result[[pred_grid+1]]$result[[par]]$SN_result90
    # BIC for individual models at tau=(0.1,0.5,0.9)
    BIC90_individual_model_mtau[par] <- BIC_cal_LeeNohPark(datasub, c(0,SN_result90,n), tau_sets=tau_sets, BIC_factor=BIC_factor)
    sum_90[[par]] <- SN_result90
  }
  estcp_BIC90_mtau <- sum_90[[which.min(BIC90_individual_model_mtau)]]
  
  min_length <- floor(min(parameter_list[,1])*n)
  # Multi-bandwidth BIC at tau=(0.1,0.5,0.9)
  estcp_multiband_BIC90_mtau <- multiband_BIC(candidates=sum_90, y=datasub, min_length=min_length, tau_sets=tau_sets, BIC_factor=BIC_factor)
  # Summary
  # Summarize result
  list(BIC_factor=BIC_factor, tau_sets=tau_sets,
       # multiple tau at tau=(0.1,0.5,0.9)
       BIC90_individual_model_mtau=BIC90_individual_model_mtau, 
       estcp_BIC90_mtau=estcp_BIC90_mtau, 
       estcp_multiband_BIC90_mtau=estcp_multiband_BIC90_mtau)
}

######## Out-of-sampe prediction (without non-crossing constraint)
contcolor <- c("black","red","#3E9F05")
forecast_horizon <- 12 # 12 is for one and two week ahead prediction

prediction_result_oneweek <- prediction_result_twoweek <- c()
prediction_result_oneweek_plot <- prediction_result_twoweek_plot <- c()
for(pred_grid in 0:num_weeks_back){
  start_point <- 1
  end_point <- reference_point-pred_grid*7
  datasub <- y[start_point:end_point]
  # print(pred_grid)
  # print(date[reference_point-pred_grid*7])
  len_datasub <- length(datasub)
  # print(paste0('Data length=', len_datasub))
  sweep_res <- prediction_result[[pred_grid+1]]$result[[5]]$SN_sweep_result # tuning parameter (0.1,0.02)
  date <- format(as.Date(country_data$date),format="%b-%d")
  plot(y, xaxt="n", xlab="", ylab="", main=paste("Forecast Made on", date[end_point]), cex.main=2, cex.lab=1.5)
  mtext(TeX("Y_t"), side = 2, line = 2,cex=1.5)
  axis(1,c(1,floor(length(y)/5)*(1:4),length(y)), date[c(1,floor(length(y)/5)*(1:4),length(y))])
  abline(v=(reference_point-pred_grid*7))
  abline(v=(reference_point-pred_grid*7+5),lty=3)
  abline(v=(reference_point-pred_grid*7+12),lty=3)
  # detected cp
  critical_sets <- c(49.88713,56.58763,68.21443)  ### tau=0.1,0.5,0.9 (grid=0.1, trim=0.02)
  x1 <- SN_local_peak(grid_size=floor(0.1*len_datasub), max_sweep_result=max_SNsweep(sweep_res), critical_value=critical_sets[1])
  mylims <- par("usr")
  points(x1, rep(mylims[3],length(x1)), col='blue', pch=17,  cex=2.5, bg='blue')
  date <- format(as.Date(country_data$date),format="%Y-%b-%d")
  x1_bic <- prediction_result_BIC[[pred_grid+1]]$estcp_multiband_BIC90_mtau$est_cp
  points(x1_bic, rep(mylims[4],length(x1_bic)), col='black', pch=25,  cex=2.5, bg='black')
  cpt <- c(0, x1, len_datasub)
  print(paste0('Forecast made on ', date[end_point]))
  print(pred_grid)
  print(round(min(parameter_list[,1])*len_datasub))
  print(cpt)
  print(diff(cpt))
  
  num_segment <- length(cpt)-1

  cont <-1
  x_insample <- ((cpt[cont]+1):(cpt[cont+1]))/len_datasub
  reg_nocross <- rq.no.cross_insample(datasub[(cpt[cont]+1):(cpt[cont+1])], x=cbind(x_insample), taus=tau_sets)
  fix <- ((cpt[cont]+1):(cpt[cont+1]))+start_point-1
  j <- 1
  lines(fix, cbind(1,x_insample)%*%reg_nocross$bhat[,j], col=contcolor[j], lwd=2)
  j <- 2
  lines(fix[1:(length(fix)-3)], cbind(1,x_insample[1:(length(fix)-3)])%*%reg_nocross$bhat[,j], col=contcolor[j], lwd=2)
  j <- 3
  lines(fix[1:(length(fix)-3)], cbind(1,x_insample[1:(length(fix)-3)])%*%reg_nocross$bhat[,j], col=contcolor[j], lwd=2)

  for(cont in 2:(num_segment-1)){
    x_insample <- ((cpt[cont]+1):(cpt[cont+1]))/len_datasub
    reg_nocross <- rq.no.cross_insample(datasub[(cpt[cont]+1):(cpt[cont+1])], x=cbind(x_insample), taus=tau_sets)
    fix <- ((cpt[cont]+1):(cpt[cont+1]))+start_point-1
    for(j in 1:length(tau_sets)){
      lines(fix, cbind(1,x_insample)%*%reg_nocross$bhat[,j], col=contcolor[j], lwd=2)
    }
  }
  # plot(max_SNsweep(sweep_res), xaxt="n", xlab="", ylab="y_t", main=paste(country_name, pred_grid), cex.main=2, cex.lab=1.5)

  # Prediction based on last segment
  cont <- num_segment
  # cpt[cont] <- cpt[cont+1]-28 # 4 weeks back moving window
  len_segment <- cpt[cont+1]-cpt[cont]
  fix <- (cpt[cont]+1):cpt[cont+1]+start_point-1
  x_insample <- ((cpt[cont]+1):cpt[cont+1])/len_datasub
  x2_insample <- x_insample^2
  x3_insample <- x_insample^3
  insample_X <- cbind(1, x_insample, x2_insample, x3_insample)
  x_pred <- (len_datasub+(1:forecast_horizon))/len_datasub
  noncrossing_horizon <- forecast_horizon # for non-crossing forecasting horizon
  forecast_h <- c((len_datasub+noncrossing_horizon)/len_datasub, ((len_datasub+noncrossing_horizon)/len_datasub)^2, ((len_datasub+noncrossing_horizon)/len_datasub)^3)
  forecast_X <- cbind(1, x_pred, x_pred^2, x_pred^3)

  # model selection (via AIC/BIC of asymmetric Laplace distribution)
  loss_linear <- loss_quadratic <- c()
  # for(tau_index in 1:length(tau_sets)){
  for(tau_index in which(tau_sets==0.5)){ # Model selection based on median regression
    reg_linear_insample <- rq(datasub[(cpt[cont]+1):(cpt[cont+1])]~x_insample, tau=tau_sets[tau_index])
    reg_quadratic_insample <- rq(datasub[(cpt[cont]+1):(cpt[cont+1])]~x_insample+x2_insample, tau=tau_sets[tau_index])
    loss_linear <- c(loss_linear, log(check_loss(y=reg_linear_insample$y, est_y=reg_linear_insample$fitted.values, tau=tau_sets[tau_index])))
    loss_quadratic <- c(loss_quadratic, log(check_loss(y=reg_quadratic_insample$y, est_y=reg_quadratic_insample$fitted.values, tau=tau_sets[tau_index])))
  }
  # estimation for linear (noncrossing)
  reg_linear <- rq.no.cross_pred(y=datasub[(cpt[cont]+1):(cpt[cont+1])], x=cbind(x_insample), forecast_h[1], tau_sets)
  # estimation for quadratic (noncrossing)
  reg_quadratic <- rq.no.cross_pred(y=datasub[(cpt[cont]+1):(cpt[cont+1])], x=cbind(x_insample,x2_insample), forecast_h[1:2], tau_sets)

  # daily prediction
  data_future <- y[(reference_point-pred_grid*7)+(1:forecast_horizon)] # true result (12 days forecast)
  pred_linear_daily <- pred_quadratic_daily <- pred_selected_daily <- exp(data_future)
  for(tau_index in 1:num_tau){
    pred_linear_daily <- cbind.data.frame(pred_linear_daily, exp(forecast_X[,1:2]%*%reg_linear$bhat[,tau_index]))
    pred_quadratic_daily <- cbind.data.frame(pred_quadratic_daily, exp(forecast_X[,1:3]%*%reg_quadratic$bhat[,tau_index]))
    # selected pred_linear or pred_quadratic
    if(sum(loss_linear) <= sum(loss_quadratic)+1/2/len_segment*log(len_segment)*length(loss_quadratic)){ # BIC
      lines(c(fix, start_point-1+len_datasub+1:forecast_horizon), rbind(insample_X[,1:2],forecast_X[,1:2])%*%reg_linear$bhat[,tau_index], col=contcolor[tau_index], lty=2,lwd=3)
      pred_selected_daily <- cbind.data.frame(pred_selected_daily, exp(forecast_X[,1:2]%*%reg_linear$bhat[,tau_index]))
    }else{
      lines(c(fix, start_point-1+len_datasub+1:forecast_horizon), rbind(insample_X[,1:3],forecast_X[,1:3])%*%reg_quadratic$bhat[,tau_index], col=contcolor[tau_index], lty=4,lwd=3)
      pred_selected_daily <- cbind.data.frame(pred_selected_daily, exp(forecast_X[,1:3]%*%reg_quadratic$bhat[,tau_index]))
    }
  }
  colnames(pred_linear_daily) <- colnames(pred_quadratic_daily) <- colnames(pred_selected_daily) <-
    c('true', tau_sets)

  # Compare with CDC (one and two-week ahead)
  data_future <- y[(reference_point-pred_grid*7)+(-1:forecast_horizon)] # true result (14 days forecast)
  one_week_true <- sum(exp(data_future)[1:7])
  two_week_true <- sum(exp(data_future)[1:14])
  pred_linear_oneweek <- pred_quadratic_oneweek <- pred_selected_oneweek <- one_week_true
  pred_linear_twoweek <- pred_quadratic_twoweek <- pred_selected_twoweek <- two_week_true
  for(tau_index in 1:num_tau){
    # pred_linear
    pred_linear_oneweek <- c(pred_linear_oneweek, sum(exp(forecast_X[1:5,1:2]%*%reg_linear$bhat[,tau_index])))
    pred_linear_twoweek <- c(pred_linear_twoweek, sum(exp(forecast_X[,1:2]%*%reg_linear$bhat[,tau_index])))
    # pred quadratic
    pred_quadratic_oneweek <- c(pred_quadratic_oneweek, sum(exp(forecast_X[1:5,1:3]%*%reg_quadratic$bhat[,tau_index])))
    pred_quadratic_twoweek <- c(pred_quadratic_twoweek, sum(exp(forecast_X[,1:3]%*%reg_quadratic$bhat[,tau_index])))
    # selected pred_linear or pred_quadratic
    if(sum(loss_linear) <= sum(loss_quadratic)+1/2/len_segment*log(len_segment)*length(loss_quadratic)){ # BIC
      pred_selected_oneweek <- c(pred_selected_oneweek, sum(exp(forecast_X[1:5,1:2]%*%reg_linear$bhat[,tau_index])))
      pred_selected_twoweek <- c(pred_selected_twoweek, sum(exp(forecast_X[,1:2]%*%reg_linear$bhat[,tau_index])))
    }else{
      pred_selected_oneweek <- c(pred_selected_oneweek, sum(exp(forecast_X[1:5,1:3]%*%reg_quadratic$bhat[,tau_index])))
      pred_selected_twoweek <- c(pred_selected_twoweek, sum(exp(forecast_X[,1:3]%*%reg_quadratic$bhat[,tau_index])))
    }
  }
  # add back Saturday and Monday for one-week and two-week ahead prediction
  pred_linear_oneweek[-1] <- pred_linear_oneweek[-1]+ sum(exp(data_future)[1:2])
  pred_quadratic_oneweek[-1] <- pred_quadratic_oneweek[-1]+ sum(exp(data_future)[1:2])
  pred_selected_oneweek[-1] <- pred_selected_oneweek[-1]+ sum(exp(data_future)[1:2])
  pred_linear_twoweek[-1] <- pred_linear_twoweek[-1]+ sum(exp(data_future)[1:2])
  pred_quadratic_twoweek[-1] <- pred_quadratic_twoweek[-1]+ sum(exp(data_future)[1:2])
  pred_selected_twoweek[-1] <- pred_selected_twoweek[-1]+ sum(exp(data_future)[1:2])

  ### CDC results only available after Aug-03 (which requires pred_grid<12)
  cdc_data <- read.csv(paste0('./CDC_Forecast/', as.Date(date[reference_point-pred_grid*7], format='%Y-%b-%d'), '-all-forecasted-cases-model-data.csv'))
  cdc_data <- cdc_data[cdc_data$State=='National',]
  cdc_data <- cdc_data[order(cdc_data$model),]
  cdc_data <- cdc_data[cdc_data$target%in%c('1 wk ahead inc case', '2 wk ahead inc case'),]
  cdc_oneweek <- cdc_data[cdc_data$target%in%c('1 wk ahead inc case'),]
  cdc_oneweek <- cbind(cdc_oneweek, one_week_true)
  cdc_twoweek <- cdc_data[cdc_data$target%in%c('2 wk ahead inc case'),]
  cdc_twoweek$point <- cdc_twoweek$point + cdc_oneweek$point[cdc_oneweek$model%in%cdc_twoweek$model] # generate two-week cumulative prediction
  cdc_twoweek$quantile_0.025 <- cdc_twoweek$quantile_0.025 + cdc_oneweek$quantile_0.025[cdc_oneweek$model%in%cdc_twoweek$model] # generate two-week cumulative prediction
  cdc_twoweek$quantile_0.975 <- cdc_twoweek$quantile_0.975 + cdc_oneweek$quantile_0.975[cdc_oneweek$model%in%cdc_twoweek$model] # generate two-week cumulative prediction

  cdc_twoweek <- cbind(cdc_twoweek, two_week_true)

  cdc_oneweek_lb <- cdc_oneweek$quantile_0.025[which(cdc_oneweek$model=='Ensemble')]
  cdc_oneweek_ub <- cdc_oneweek$quantile_0.975[which(cdc_oneweek$model=='Ensemble')]
  cdc_twoweek_lb <- cdc_twoweek$quantile_0.025[which(cdc_twoweek$model=='Ensemble')]
  cdc_twoweek_ub <- cdc_twoweek$quantile_0.975[which(cdc_twoweek$model=='Ensemble')]

  tmp_prediction_result_oneweek <- c(date[reference_point-pred_grid*7], date[reference_point-pred_grid*7+5], one_week_true,
                                     round(pred_selected_oneweek[3]), round(pred_selected_oneweek[2]), round(pred_selected_oneweek[4]),
                                     pmin(pred_selected_oneweek[1]>pred_selected_oneweek[2], pred_selected_oneweek[1]<pred_selected_oneweek[4]),
                                     cdc_oneweek$point[which(cdc_oneweek$model=='Ensemble')],
                                     cdc_oneweek_lb, cdc_oneweek_ub, pmin(pred_selected_oneweek[1]>cdc_oneweek_lb, pred_selected_oneweek[1]<cdc_oneweek_ub))
  tmp_prediction_result_oneweek_relerror <- c(NA, NA, NA, (pred_selected_oneweek[3]-one_week_true)/one_week_true, NA, NA, NA,
                                              (cdc_oneweek$point[which(cdc_oneweek$model=='Ensemble')]-one_week_true)/one_week_true, NA, NA, NA)
  tmp_prediction_result_twoweek <- c(date[reference_point-pred_grid*7], date[reference_point-pred_grid*7+12], two_week_true,
                                     round(pred_selected_twoweek[3]), round(pred_selected_twoweek[2]), round(pred_selected_twoweek[4]),
                                     pmin(pred_selected_twoweek[1]>pred_selected_twoweek[2], pred_selected_twoweek[1]<pred_selected_twoweek[4]),
                                     cdc_twoweek$point[which(cdc_twoweek$model=='Ensemble')],
                                     cdc_twoweek_lb, cdc_twoweek_ub, pmin(pred_selected_twoweek[1]>cdc_twoweek_lb, pred_selected_twoweek[1]<cdc_twoweek_ub))
  tmp_prediction_result_twoweek_relerror <- c(NA, NA, NA, (pred_selected_twoweek[3]-two_week_true)/two_week_true, NA, NA, NA,
                                              (cdc_twoweek$point[which(cdc_twoweek$model=='Ensemble')]-two_week_true)/two_week_true, NA, NA, NA)

  prediction_result_oneweek <- rbind(prediction_result_oneweek, tmp_prediction_result_oneweek, round(tmp_prediction_result_oneweek_relerror*100,2))
  prediction_result_twoweek <- rbind(prediction_result_twoweek, tmp_prediction_result_twoweek, round(tmp_prediction_result_twoweek_relerror*100,2))
  prediction_result_oneweek_plot <- rbind(prediction_result_oneweek_plot, tmp_prediction_result_oneweek)
  prediction_result_twoweek_plot <- rbind(prediction_result_twoweek_plot, tmp_prediction_result_twoweek)
}
# dev.off()
options(scipen=2)
# postscript('Oneweek_fcst.eps', width=8, height=8)
prediction_result_oneweek_plot <- data.frame(prediction_result_oneweek_plot)
plot_date <- rev(format(as.Date(as.character(prediction_result_oneweek_plot$X1), format="%Y-%b-%d"), format="%b-%d"))
plot(rev(as.numeric(as.character(prediction_result_oneweek_plot$X3))),xaxt="n",
     ylim=c(200000,1100000), pch=19, col='red', main='One week ahead forecast results',
     ylab="", xlab='', cex.main=2,cex=1.5)
mtext('Cumulative New Cases', side = 2, line = 2,cex=1.5)
axis(1,seq(1,length(plot_date),by=2),plot_date[seq(1,length(plot_date),by=2)])
points(rev(as.numeric(as.character(prediction_result_oneweek_plot$X4))), pch=2, col='blue', cex=1.5)
lines(rev(as.numeric(as.character(prediction_result_oneweek_plot$X5))), col='blue',type="l")
lines(rev(as.numeric(as.character(prediction_result_oneweek_plot$X6))), col='blue',type="l")

points(rev(as.numeric(as.character(prediction_result_oneweek_plot$X8))), pch=4, cex=1.5)
lines(rev(as.numeric(as.character(prediction_result_oneweek_plot$X9))),lty=2)
lines(rev(as.numeric(as.character(prediction_result_oneweek_plot$X10))),lty=2)
legend('topleft', c('CDC', 'GOALS', 'True','CDC 95% CI','GOALS 80% CI'), pch=c(4,2,19,NA,NA),lty=c(NA,NA,NA,2,1), col=c('black', 'blue', 'red','black', 'blue'),pt.cex=1.5)
# dev.off()

# postscript('Twoweek_fcst.eps', width=8, height=8)
prediction_result_twoweek_plot <- data.frame(prediction_result_twoweek_plot)
plot_date <- rev(format(as.Date(as.character(prediction_result_twoweek_plot$X1), format="%Y-%b-%d"), format="%b-%d"))
plot(rev(as.numeric(as.character(prediction_result_twoweek_plot$X3))),xaxt="n",
     ylim=c(3e05,2.4e6), pch=19, col='red', main='Two week ahead forecast results',
     ylab="", xlab='', cex.main=2,cex=1.5)
mtext('Cumulative New Cases', side = 2, line = 2,cex=1.5)
axis(1,seq(1,length(plot_date),by=2),plot_date[seq(1,length(plot_date),by=2)])
points(rev(as.numeric(as.character(prediction_result_twoweek_plot$X4))), pch=2, col='blue', cex=1.5)
lines(rev(as.numeric(as.character(prediction_result_twoweek_plot$X5))), col='blue',type="l")
lines(rev(as.numeric(as.character(prediction_result_twoweek_plot$X6))), col='blue',type="l")

points(rev(as.numeric(as.character(prediction_result_twoweek_plot$X8))), pch=4, cex=1.5)
lines(rev(as.numeric(as.character(prediction_result_twoweek_plot$X9))),lty=2)
lines(rev(as.numeric(as.character(prediction_result_twoweek_plot$X10))),lty=2)
legend('topleft', c('CDC', 'GOALS', 'True','CDC 95% CI','GOALS 80% CI'), pch=c(4,2,19,NA,NA),lty=c(NA,NA,NA,2,1), col=c('black', 'blue', 'red','black', 'blue'),pt.cex=1.5)
# dev.off()

######## summarize results
prediction_result_oneweek <- prediction_result_oneweek # remove 1st week no true data
prediction_result_twoweek <- prediction_result_twoweek # remove 1st and 2nd week no true data
prediction_result_oneweek <- data.frame(prediction_result_oneweek)
prediction_result_twoweek <- data.frame(prediction_result_twoweek)
rownames(prediction_result_oneweek) <- rownames(prediction_result_twoweek) <- NULL
colnames(prediction_result_oneweek) <- colnames(prediction_result_twoweek) <- c('Date', 'Target', 'True', 'SN_Median', '10% Q', '90% Q', 'Covered',
                                                                                'CDC_Ensemble', '2.5% Q', '97.5% Q', 'Covered')
print(prediction_result_oneweek)
print(prediction_result_twoweek)
# print(xtable(prediction_result_oneweek), include.rownames=FALSE)
# print(xtable(prediction_result_twoweek), include.rownames=FALSE)
