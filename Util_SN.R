####generate autogressive models
MAR <- function(n, reptime, rho){   
  inter <- matrix(0, n+30, reptime)
  epsilon <- matrix(rnorm((n+30)*reptime,0,1),(n+30),reptime)
  for (j in 1:(n+29)){
    inter[j+1,1:reptime] <- epsilon[j+1,1:reptime]+rho*inter[j,1:reptime]
  }
  return(inter[31:(n+30),1:reptime])
}

####generate the piecewise linear trend location mmodel
LinearTrend <- function(n, reptime, cp_sets, rho, lines_sets, sigma=1){    
  no_seg <- length(cp_sets)-1
  inter <- MAR(n, reptime, rho)*sqrt(1-rho^2)*sigma # standardized
  for(index in 1:no_seg){ # Mean shift
    a <- lines_sets[[index]][1]
    b <- lines_sets[[index]][2]
    tau1 <- cp_sets[index]+1
    tau2 <- cp_sets[index+1]
    inter[tau1:tau2,] <- inter[tau1:tau2,] + a + b*(tau1:tau2)
  }
  return(inter)
}

####homoscedastic
gen_data_homo <- function(reptime, n, rho, sigma, cpt_relative, b_slope, b0, null=F){
  data <- matrix(0,nrow=reptime,ncol=n) 
  for (i in 1:reptime) {
    set.seed(i)
    u <- rep(0,n)
    u[1] <- rnorm(1)*sigma # initial
    epsilon <- rnorm(n,0,1)*sigma*sqrt(1-rho^2)
    for (t in 2:n) {
      u[t] <- rho*u[t-1]+epsilon[t]
    }
    y <- rep(0,n)
    intercept <- b0
    for (j in 1:(length(cpt_relative)-1)){
      y[(round(n*cpt_relative[j])+1):round(n*cpt_relative[j+1])] <- intercept+b_slope[j]*(((round(n*cpt_relative[j])+1):round(n*cpt_relative[j+1]))/n-cpt_relative[j])#+u[1:(n*0.12)]
      intercept <- intercept+b_slope[j]*(cpt_relative[j+1]-cpt_relative[j])
    }
    if(!null){
      data[i,] <- y+u
    }else{
      data[i,] <- u
    }
  }
  return(data)
}

############### heteroscedastic
gen_data_hetero <- function(reptime, n, rho, sigma, cpt_relative, b_slope, b0, null=F){
  data <- matrix(0,nrow=reptime,ncol=n) 
  for (i in 1:reptime) {
    set.seed(i)
    u <- rep(0,n)
    u[1] <- rnorm(1)*sigma # initial
    epsilon <- rnorm(n,0,1)*sigma*sqrt(1-rho^2)
    for (t in 2:n) {
      u[t] <- rho*u[t-1]+epsilon[t]
    }
    u <- u*(1+0.3*(1:n)/n)###this part for location scale
    y <- rep(0,n)
    intercept <- b0
    for (j in 1:(length(cpt_relative)-1)){
      y[(round(n*cpt_relative[j])+1):round(n*cpt_relative[j+1])] <- intercept+b_slope[j]*(((round(n*cpt_relative[j])+1):round(n*cpt_relative[j+1]))/n-cpt_relative[j])#+u[1:(n*0.12)]
      intercept <- intercept+b_slope[j]*(cpt_relative[j+1]-cpt_relative[j])
    }
    if(!null){
      data[i,] <- y+u
    }else{
      data[i,] <- u
    }
  }
  return(data)
}

###Global testing part find the global testing statistics at each potential location
max_SNsweep <- function(SN_sweep_result) {     
  max_matrix <- function(matrix_data) {
    if (is.null(dim(matrix_data))) {
      return(0)
    } else{
      return(max(matrix_data[, 1]))
    }
  }
  return(sapply(SN_sweep_result, max_matrix))
}

###local scanning part
SN_local_peak <- function(grid_size, max_sweep_result, critical_value){   
  candidate <- NULL
  over_peak <- max_sweep_result>critical_value
  search_set <- which(over_peak==1)
  for (i in search_set)
  {
    sub_set <- (i-grid_size+1):(i+grid_size)
    if (max_sweep_result[i]==max(max_sweep_result[sub_set]))
    {
      candidate <- c(candidate,i)
    }
  }
  return(unique(sort(candidate)))
}

#### Function to calculate the pairwise quantile regression
SN_pairwise_rq <- function(data, trim_size, total_length, tau_sets){
  total_combn <- combn(total_length, 2)
  dist <- total_combn[2,]-total_combn[1,]
  total_combn_selected <- total_combn[, dist>trim_size]
  num_tau <- length(tau_sets)
  time_index <- (1:total_length)/total_length
  rq_pairwise_matrix_intercept <- rq_pairwise_matrix_slope <- list()
  for(tau_index in 1:num_tau){
    rq_pairwise_matrix_intercept[[tau_index]] <- rq_pairwise_matrix_slope[[tau_index]] <- matrix(0, total_length, total_length)
  }
  
  for(pair_index in 1:dim(total_combn_selected)[2]){
    pair_pos <- total_combn_selected[,pair_index]
    loc1 <- pair_pos[1]; loc2 <- pair_pos[2]
    tmp_rq <- rq(data[loc1:loc2]~time_index[loc1:loc2], tau=tau_sets)
    if(num_tau==1){
      rq_pairwise_matrix_intercept[[1]][loc1,loc2] <- rq_pairwise_matrix_intercept[[1]][loc2,loc1] <- tmp_rq$coefficients[1]
      rq_pairwise_matrix_slope[[1]][loc1,loc2] <- rq_pairwise_matrix_slope[[1]][loc2,loc1] <- tmp_rq$coefficients[2]
    }else{
      for(tau_index in 1:num_tau){
        rq_pairwise_matrix_intercept[[tau_index]][loc1,loc2] <- tmp_rq$coefficients[1,tau_index]
        rq_pairwise_matrix_slope[[tau_index]][loc1,loc2] <- tmp_rq$coefficients[2,tau_index]
      }
    }
  }
  return(list(rq_pairwise_matrix_intercept=rq_pairwise_matrix_intercept, rq_pairwise_matrix_slope=rq_pairwise_matrix_slope))
}