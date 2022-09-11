rm(list=ls())
source("./Util_single.R")
source("../Util_pool.R")
source("../Util_SN.R")
library(doParallel)
library(foreach)
library(quantreg)
# registerDoParallel(24)

### DGP settings
code <- 'All'
ts_settings <- list(
  list(cpt_relative=c(0,0.15,0.35,0.6,0.8,1), ###Australia setting DGP3
       b_slope=c(-23,-3,19,-14,-3), b0=6),
  list(cpt_relative=c(0,0.1,0.4,0.6,0.8,1),   ###US setting DGP2
       b_slope=c(45,-1.3,7.5,-4,6), b0=5.5),
  list(cpt_relative=c(0,0.10,0.3,0.55,0.8,1), ###UK setting  DGP1 
       b_slope=c(26,2,-7,4,11), b0=5.8))

### generating simulation date
reptime <- 500
n_sets <- 210
sigma_sets <- c(0.1,0.2)
rho_sets <- c(-0.5,-0.3,0,0.3,0.5)
tau_sets <- c(0.1,0.5,0.9)
hetero_settings <- c(T, F) # heteroscedastic error
null_settings <- c(T, F) # null or with change-points


###################### begin simulation with single quantile
critical_sets <- c(65.41220,74.62814,97.90169)   ###pivotal critical value at 90% 95% 99% levels (epsilon, delta)=(0.1,0.02)
for(null in null_settings){
  for(hetero in hetero_settings){
    for(n in n_sets){
      counter <- 0
      final_result <- list()
      grid_size <- round(n*0.1) # round grid size
      trim_size <- ceiling(n*0.02) # ceiling trim size
      for(ts_setting in ts_settings){
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
            for(tau in tau_sets){
              result <- foreach(j=1:reptime, .packages='quantreg') %dopar% {
                y <- data[j,]
                SN_sweep_result <- SN_sweep_quantile_single(data=y, grid_size=grid_size, trim_size=trim_size, tau=tau)
                max_SN <- max_SNsweep(SN_sweep_result)
                SN_result90 <- SN_local_peak(grid_size, max_SN, critical_value=critical_sets[1])
                SN_result95 <- SN_local_peak(grid_size, max_SN, critical_value=critical_sets[2])
                SN_result99 <- SN_local_peak(grid_size, max_SN, critical_value=critical_sets[3])
                list(SN_result90=SN_result90, SN_result95=SN_result95, SN_result99=SN_result99)
              }
              counter <- counter + 1
              final_result[[counter]] <- list(n=n, rho=rho, sigma=sigma, tau=tau, cpt_relative=cpt_relative, result=result)
              # write.table(t(c(null,hetero,n,sigma,rho,tau,cpt_relative)), "alg5_result1.txt", append=TRUE, col.names=F)  ##for monitoring
            }
          }
        }
      }
      if(hetero){
        if(!null){
          save.image(paste0(code, '_quantile_single_simu_hetero_n', n, '.RData'))
        }else{
          save.image(paste0(code, '_quantile_single_simu_hetero_null_n', n, '.RData'))
        }
      }else{
        if(!null){
          save.image(paste0(code, '_quantile_single_simu_homo_n', n, '.RData'))
        }else{
          save.image(paste0(code, '_quantile_single_simu_homo_null', n, '.RData'))
        }
      } 
    }
  }
}


###################### begin simulation with multiple quantiles
tau_sets <- c(0.1,0.5,0.9)
critical_sets <- c(49.88713,56.58763,68.21443)   ### pivotal critical value at 90% 95% 99% levels (epsilon, delta)=(0.1,0.02)
for(null in null_settings){
  for(hetero in hetero_settings){
    for(n in n_sets){
      counter <- 0 
      final_result <- list()
      grid_size <- round(n*0.1) # round grid size
      trim_size <- ceiling(n*0.02) # ceiling trim size
      for(ts_setting in ts_settings){
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
            result <- foreach(j=1:reptime, .packages='quantreg') %dopar% {
              y <- data[j,]
              time1 <- proc.time()
              SN_sweep_result <- SN_sweep_quantile_pool(data=y, grid_size=grid_size, trim_size=trim_size, tau_sets=tau_sets)
              time2 <- proc.time()
              max_SN <- max_SNsweep(SN_sweep_result)
              SN_result90 <- SN_local_peak(grid_size, max_SN, critical_value=critical_sets[1])
              SN_result95 <- SN_local_peak(grid_size, max_SN, critical_value=critical_sets[2])
              SN_result99 <- SN_local_peak(grid_size, max_SN, critical_value=critical_sets[3])
              time <- as.numeric(time2-time1)[1]
              
              # plot(max_SN)
              # abline(h=critical_sets[1])
              list(SN_result90=SN_result90, SN_result95=SN_result95, SN_result99=SN_result99, time=time)
            }
            counter <- counter + 1
            final_result[[counter]] <- list(n=n, rho=rho, sigma=sigma, tau_sets=tau_sets, cpt_relative=cpt_relative, result=result)
            # write.table(t(c(null,hetero,n,sigma,rho,tau_sets,cpt_relative)), "alg5_result1_pool.txt", append=TRUE, col.names=F)
          }
        }
      }
      if(hetero){
        if(!null){
          save.image(paste0(code, '_quantile_pool_simu_hetero_n', n, '.RData'))
        }else{
          save.image(paste0(code, '_quantile_pool_simu_hetero_null_n', n, '.RData'))
        }
      }else{
        if(!null){
          save.image(paste0(code, '_quantile_pool_simu_homo_n', n, '.RData'))
        }else{
          save.image(paste0(code, '_quantile_pool_simu_homo_null', n, '.RData'))
        }
      }
    }
  }
}
