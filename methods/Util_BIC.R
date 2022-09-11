check_loss <- function(y, est_y, tau){
  # Function to calculate check loss
  u <- y-est_y
  check <- sum(u*(tau-(u<0)))
  return(check)
}

BIC_cal_LeeNohPark <- function(y, est_cp, tau_sets, BIC_factor=0){
  # Function to calculate the BIC (Lee et al., 2014) of a given est_cp estimation,
  n <- length(y)
  no_seg <- length(est_cp)-1
  total_BIC <- 0
  for(tau in tau_sets){
    check_losses <- 0
    for(seg_index in 1:no_seg){
      start <- est_cp[seg_index]+1
      end <- est_cp[seg_index+1]
      x_insample <- (start:end)/n
      reg_linear_insample <- rq(y[start:end]~x_insample, tau=tau)
      check_losses <- check_losses + check_loss(y=reg_linear_insample$y, est_y=reg_linear_insample$fitted.values, tau=tau)
    }
    BIC <- log(check_losses) + 1/2*((2+1)*no_seg)*(log(n))^(1+BIC_factor)/n
    total_BIC <- total_BIC + BIC
  }
  return(total_BIC)
}

multiband_BIC <- function(candidates, y, min_length, tau_sets, BIC_factor=0, max_num_cand=12){
  # Function to heuristcally find the best combination of change-points across several bandwidth
  n <- length(y)
  
  # individual model
  BIC_result_individual <- c()
  num_cand <- length(candidates)
  for(cand_index in 1:num_cand){
    BIC_result_individual <- c(BIC_result_individual, BIC_cal_LeeNohPark(y, c(0, candidates[[cand_index]], n), tau_sets, BIC_factor))
  }
  
  # filtering + exhaustive search
  BIC_result_multiband <- c()
  candidates_cp <- unique(sort(unlist(candidates)))
  num_para_list <- length(candidates) # num of combination of (grid, trim)
  while(length(candidates_cp)>max_num_cand & num_para_list>1){
    # only select top-performing (grid, trim) by BIC to reduce consideration
    candidates_cp <- unique(sort(unlist(candidates[which(BIC_result_individual<sort(BIC_result_individual)[num_para_list])])))
    num_para_list <- num_para_list-1
  }
  num_cand <- 2^length(candidates_cp)
  poss_combination <- expand.grid(rep(list(0:1), length(candidates_cp)))
  for(cand_index in 1:num_cand){
    tmp_cand_cp <- c(0,candidates_cp[as.logical(poss_combination[cand_index,])],n)
    if(min(diff(tmp_cand_cp))<min_length){
      BIC_result_multiband <- c(BIC_result_multiband, Inf)
    }else{
      BIC_result_multiband <- c(BIC_result_multiband, BIC_cal_LeeNohPark(y, tmp_cand_cp, tau_sets, BIC_factor))
    }
  }
  
  # stepwise selection
  previous_cp <- c(0,n)
  previous_BIC_result <- BIC_cal_LeeNohPark(y, previous_cp, tau_sets, BIC_factor)
  current_candidates_cp <- unique(sort(unlist(candidates)))
  continue <- T
  while(continue){
    continue <- F
    num_cand <- length(current_candidates_cp)
    tmp_BIC_result <- c()
    for(cand_index in 1:num_cand){
      current_cp <- sort(c(current_candidates_cp[cand_index], previous_cp))
      if(min(diff(current_cp))<min_length){
        tmp_BIC_result <- c(tmp_BIC_result, Inf)
      }else{
        tmp_BIC_result <- c(tmp_BIC_result, BIC_cal_LeeNohPark(y, current_cp, tau_sets, BIC_factor))
      }
    }
    current_BIC_result <- min(tmp_BIC_result)
    if(current_BIC_result<previous_BIC_result){
      current_candidates_cp <- current_candidates_cp[-which.min(tmp_BIC_result)]
      previous_cp <- sort(c(current_candidates_cp[which.min(tmp_BIC_result)], previous_cp))
      previous_BIC_result <- current_BIC_result
      if(length(current_candidates_cp)>0){
        continue <- T
      }
    }
  }
  
  # combine together
  winner <- which.min(c(min(BIC_result_individual), min(BIC_result_multiband), previous_BIC_result))
  min_BIC <- min(c(min(BIC_result_individual), min(BIC_result_multiband), previous_BIC_result))
  if(winner==1){
    est_cp <- candidates[[which.min(BIC_result_individual)]]
  }
  if(winner==2){
    est_cp <- candidates_cp[as.logical(poss_combination[which.min(BIC_result_multiband),])]
  }
  if(winner==3){
    est_cp <- previous_cp[-c(1,length(previous_cp))]
  }
  return(list(est_cp=est_cp, min_BIC=min_BIC))
}


