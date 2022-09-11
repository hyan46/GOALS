SN_sweep_quantile_single <- function(data, grid_size, trim_size, tau){
  # Define cumsum function
  cumsum_quantile_contrast_single <- function(ts_length, type="L", start_t, trim_size){
    n <- ts_length
    end_t <- start_t+n-1
    a_result <- b_result <- matrix(0, ncol=1, nrow=n)
    if(type=="L"){
      for (i in 2:(n-2)){
        a_result[i,1] <- rq_pairwise_matrix_intercept[[1]][start_t,(start_t+i-1)]-rq_pairwise_matrix_intercept[[1]][(start_t+i),end_t]
        b_result[i,1] <- rq_pairwise_matrix_slope[[1]][start_t,(start_t+i-1)]-rq_pairwise_matrix_slope[[1]][(start_t+i),end_t]
      }
    }
    if(type=="R"){
      for (i in 3:(n-1)){
        a_result[i,1] <- rq_pairwise_matrix_intercept[[1]][(start_t+i-1),end_t]-rq_pairwise_matrix_intercept[[1]][start_t,(start_t+i-2)]
        b_result[i,1] <- rq_pairwise_matrix_slope[[1]][(start_t+i-1),end_t]-rq_pairwise_matrix_slope[[1]][start_t,(start_t+i-2)]
      }
    }
    
    if(type=='L'){
      b_result[2:(2+trim_size-1),] <- 0
      b_result[(n-2-trim_size+1):(n-2),] <- 0
      a_result[2:(2+trim_size-1),] <- 0
      a_result[(n-2-trim_size+1):(n-2),] <- 0
    }else{
      b_result[3:(3+trim_size-1),] <- 0
      b_result[(n-1-trim_size+1):(n-1),] <- 0
      a_result[3:(3+trim_size-1),] <- 0
      a_result[(n-1-trim_size+1):(n-1),] <- 0
    }
    return(cbind(a_result, b_result))
  }
  
  # Define SN_test function
  SN_test_quantile_single <- function(ts_length, k, pre_grid_position, trim_size){
    n <- ts_length
    start_t <- pre_grid_position
    end_t <- start_t+n-1
    
    D1 <- c(rq_pairwise_matrix_intercept[[1]][start_t,(start_t+k-1)], rq_pairwise_matrix_slope[[1]][start_t,(start_t+k-1)])
    D2 <- c(rq_pairwise_matrix_intercept[[1]][(start_t+k),end_t], rq_pairwise_matrix_slope[[1]][(start_t+k),end_t])
    
    D <- k*(n-k)/n^1.5*(D1-D2)
    inter1 <- inter2 <- NULL
    inter1 <- cumsum_quantile_contrast_single(ts_length=k, 'L', pre_grid_position, trim_size)
    inter2 <- cumsum_quantile_contrast_single(ts_length=n-k, 'R', pre_grid_position+k, trim_size)
    
    multiplier1 <- ((1:k)*((k-1):0))^2/n^2/k^2
    multiplier2 <- ((0:(n-k-1))*((n-k):1))^2/n^2/(n-k)^2
    
    no_stat <- 2
    M1 <- M2 <- matrix(0, no_stat, no_stat)
    for(index1 in 1:(no_stat-1)){
      for(index2 in (index1+1):no_stat){
        M1[index1, index2] <- M1[index2, index1] <- sum(inter1[,index1]*inter1[,index2]*multiplier1)
        M2[index1, index2] <- M2[index2, index1] <- sum(inter2[,index1]*inter2[,index2]*multiplier2)
      }
    }
    for(index1 in 1:no_stat){
      M1[index1, index1] <- sum(inter1[,index1]^2*multiplier1)
      M2[index1, index1] <- sum(inter2[,index1]^2*multiplier2)
    }
    test_SN <- t(D)%*%solve(M1+M2)%*%D
    return(test_SN)
  }
  
  #### SN test starts
  n <- length(data)
  # Calculate all the pairwise rq results
  rq_pairwise_result <- SN_pairwise_rq(data=data, trim_size=trim_size, total_length=n, tau_sets=tau)
  rq_pairwise_matrix_intercept <- rq_pairwise_result$rq_pairwise_matrix_intercept
  rq_pairwise_matrix_slope <- rq_pairwise_result$rq_pairwise_matrix_slope
  
  # SN test for single quantile
  substat <- list()
  substat[1:(grid_size-1)] <- NA
  for (k in grid_size:(n-grid_size)){
    pre_grid_no <- floor(k/grid_size)
    post_grid_no <- floor((n-k)/grid_size)
    pre_grid_sets <- k-(pre_grid_no:1)*grid_size+1
    post_grid_sets <- k+(1:post_grid_no)*grid_size
    sn_grid_stat <- c()
    for(pre_grid_position in pre_grid_sets){
      for(post_grid_position in post_grid_sets){
        sn_grid_stat <- rbind(sn_grid_stat,
                              c(SN_test_quantile_single(ts_length=post_grid_position-pre_grid_position+1, k=k-pre_grid_position+1, pre_grid_position, trim_size),
                                pre_grid_position, k, post_grid_position))
      }
    }
    substat[[k]] <- sn_grid_stat
  }
  substat[(n-grid_size+1):n] <- NA
  return(substat)
}
