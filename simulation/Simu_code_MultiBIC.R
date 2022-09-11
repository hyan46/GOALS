rm(list=ls())
source("../Util_pool.R")
source("../Util_SN.R")
source("../Util_BIC.R")
library(doParallel)
library(foreach)
library(quantreg)
# registerDoParallel(50)

### DGP
ts_settings <- list(
  list(cpt_relative=c(0,0.15,0.35,0.6,0.8,1), ###Australia setting DGP3
       b_slope=c(-23,-3,19,-14,-3), b0=6),
  list(cpt_relative=c(0,0.1,0.4,0.6,0.8,1),   ###US setting DGP2
       b_slope=c(45,-1.3,7.5,-4,6), b0=5.5),
  list(cpt_relative=c(0,0.10,0.3,0.55,0.8,1), ###UK setting  DGP1 
       b_slope=c(26,2,-7,4,11), b0=5.8))

### generating simulation date
reptime <- 1000
n_sets <- 210
sigma_sets <- c(0.2)
rho_sets <- 0.3
parameter_list <- cbind(c(0.06,rep(c(0.08,0.1,0.12,0.15),each=2)), c(0.01,rep(c(0.01,0.02),4))) # tuning parameter (epsilon, delta)
tau_sets <- c(0.1,0.5,0.9)
hetero <- c(F,T) # heteroscedastic error
null <- c(F) # null or with change-points
critical_table <- read.csv("../SN_MultipleSlope_CriticalValue.csv")

### begin simulation
for(n in n_sets){
  for(ts_settingnum in 1:length(ts_settings)){
    ts_setting <- ts_settings[[ts_settingnum]]
    cpt_relative <- ts_setting$cpt_relative
    b_slope <- ts_setting$b_slope
    b0 <- ts_setting$b0
    for(sigma in sigma_sets){
      for(rho in rho_sets){
        if(hetero){
          data <- gen_data_hetero(reptime, n, rho, sigma, cpt_relative, b_slope, b0, null=null)
        }else{
          data <- gen_data_homo(reptime, n, rho, sigma, cpt_relative, b_slope, b0, null=null)
        }
        ### Main function
        final_result <- foreach(j=1:reptime, .packages='quantreg') %dopar% {
          BIC_factor <- 0.1 # strengthend BIC with (logN)^(1+BIC_factor)
          result <- sum_90 <- sum_95 <- list()
          BIC90_individual_model  <- c()
          BIC95_individual_model  <- c()
          y <- data[j,]
          for (par in 1:nrow(parameter_list)){
            grid <- parameter_list[par,1]
            trim <- parameter_list[par,2]
            critical_sets <- subset(critical_table, epsilon==grid&delta==trim)
            critical_sets <- unlist(critical_sets)[3:5]
            grid_size <- round(n*grid) # round grid size
            trim_size <- ceiling(n*trim) # ceiling trim size
            
            SN_sweep_result <- SN_sweep_quantile_pool(data=y, grid_size=grid_size, trim_size=trim_size, tau_sets=tau_sets)
            max_SN <- max_SNsweep(SN_sweep_result)
            SN_result90 <- SN_local_peak(grid_size, max_SN, critical_value=critical_sets[1])
            SN_result95 <- SN_local_peak(grid_size, max_SN, critical_value=critical_sets[2])
            SN_result99 <- SN_local_peak(grid_size, max_SN, critical_value=critical_sets[3])
            
            # BIC for individual models
            BIC90_individual_model[par] <- BIC_cal_LeeNohPark(y, c(0,SN_result90,n), tau_sets=tau_sets, BIC_factor=BIC_factor)
            BIC95_individual_model[par] <- BIC_cal_LeeNohPark(y, c(0,SN_result95,n), tau_sets=tau_sets, BIC_factor=BIC_factor)
            
            # summary
            result[[par]] <- list(SN_result90=SN_result90, SN_result95=SN_result95, SN_result99=SN_result99)
            sum_90[[par]] <- SN_result90
            sum_95[[par]] <- SN_result95
          }

          # Best BIC among individual models
          estcp_BIC90 <- result[[which.min(BIC90_individual_model[1:par])]]$SN_result90
          estcp_BIC95 <- result[[which.min(BIC95_individual_model[1:par])]]$SN_result95
          
          # Multi-bandwidth BIC
          min_length <- round(min(parameter_list[,1])*n)
          estcp_multiband_BIC90 <- multiband_BIC(candidates=sum_90, y=y, min_length=min_length, tau_sets=tau_sets, BIC_factor=BIC_factor)
          estcp_multiband_BIC95 <- multiband_BIC(candidates=sum_95, y=y, min_length=min_length, tau_sets=tau_sets, BIC_factor=BIC_factor)
          
          # Summarize result
          list(result=result,
               BIC90_individual_model=BIC90_individual_model, 
               BIC95_individual_model=BIC95_individual_model, 
               estcp_BIC90=estcp_BIC90,  estcp_BIC95=estcp_BIC95, 
               estcp_multiband_BIC90=estcp_multiband_BIC90, 
               estcp_multiband_BIC95=estcp_multiband_BIC95)
        }
        save.image(paste0('Setting', ts_settingnum, '_quantile_pool_simu_multiBIC_', 'hetero', hetero, '_sigma', sigma, '_rho', rho, '.RData'))
      }
    }
  }
}
