rm(list=ls())
library(doParallel)
library(foreach)
library(quantreg)
library(MASS)
library(latex2exp)
registerDoParallel(6)
source("../Util_pool.R")
source("../Util_SN.R")
source("../Util_BIC.R")
source("./Util_noncrossing.R")

######## multiple grid-trim
parameter_list <- cbind(c(0.06,rep(c(0.08,0.1,0.12,0.15),each=2)), c(0.01,rep(c(0.01,0.02),4)))
critical_table <- read.csv("../SN_MultipleSlope_CriticalValue.csv")
tau_sets <- c(0.1,0.5,0.9)
num_tau <- length(tau_sets)

######## in-sample change-point analysis
full_data <- read.csv("./owid-covid-data.csv")
last_day <- '2020-11-07' # last day of analysis
Nov07data <- subset(full_data,date==last_day)
Nov07data <- Nov07data[order(Nov07data$total_cases,decreasing = T),]


#### 35 countries, used in Figure 2
G20_country <- c("Argentina","Australia","Brazil","Canada","China","France","Germany","India","Indonesia","Italy",
                 "Japan","Mexico","Russia","Saudi Arabia","South Africa","South Korea","Turkey","United Kingdom","United States")
top_country <- Nov07data$location[2:31]
country_list <- union(top_country,G20_country)

######## SNCP at different (trim, grid)
country_result <- foreach(country_index=1:length(country_list), .packages='quantreg') %dopar% {
  country_name <- country_list[country_index]
  country_data <- subset(full_data, (location==country_name)&(total_cases>1000)&(new_cases>=1))
  last_day_reference_point <- max(which(as.Date(country_data$date)<=as.Date(last_day))) # some country does not have last day data
  y <- log(country_data$new_cases)[1:last_day_reference_point]
  y[is.na(y)] <- 0
  n <- length(y)
  ### Sensitivity analysis
  result <- list()
  for (par in 1:nrow(parameter_list)){
    grid <- parameter_list[par,1]
    trim <- parameter_list[par,2]
    critical_sets <- subset(critical_table, epsilon==grid&delta==trim)
    critical_sets <- unlist(critical_sets)[3:5]
    grid_size <- floor(n*grid) # round grid size
    trim_size <- ceiling(n*trim) # ceiling trim size
    SN_sweep_result <- tryCatch({SN_sweep_quantile_pool(data=y, grid_size=grid_size, trim_size=trim_size, tau_sets=tau_sets)},
                                error=function(cond){return(NULL)})
    if(is.null(SN_sweep_result)){
      SN_result90 <- SN_result95 <- NULL
    }else{
      max_SN <- max_SNsweep(SN_sweep_result)
      SN_result90 <- SN_local_peak(grid_size, max_SN, critical_value=critical_sets[1])
      SN_result95 <- SN_local_peak(grid_size, max_SN, critical_value=critical_sets[2])
    }
    result[[par]] <- list(grid=grid, trim=trim, SN_result90=SN_result90, SN_result95=SN_result95, SN_sweep_result=SN_sweep_result)
  }
  list(result=result)
}

######## SNCP selected by  BIC
max_num_cand <- 12 # maximum number of candidates to consider in multi-bandwidth BIC
country_result_BIC <- foreach(country_index=1:length(country_list), .packages='quantreg') %dopar% {
  country_name <- country_list[country_index]
  country_data <- subset(full_data, (location==country_name)&(total_cases>1000)&(new_cases>=1))
  last_day_reference_point <- max(which(as.Date(country_data$date)<=as.Date(last_day))) # some country does not have last day data
  y <- log(country_data$new_cases)[1:last_day_reference_point]
  y[is.na(y)] <- 0
  n <- length(y)
  # BIC selection
  BIC_factor <- 0.1 # strengthend BIC with (logN)^(1+BIC_factor)
  sum_90 <- list()
  BIC90_individual_model_mtau <- c()
  for (par in 1:nrow(parameter_list)){
    SN_result90 <- country_result[[country_index]]$result[[par]]$SN_result90
    # BIC for individual models at tau=(0.1,0.5,0.9)
    BIC90_individual_model_mtau[par] <- BIC_cal_LeeNohPark(y, c(0,SN_result90,n), tau_sets=tau_sets, BIC_factor=BIC_factor)
    sum_90[[par]] <- SN_result90
  }
  # Best BIC among individual models
  estcp_BIC90_mtau <- sum_90[[which.min(BIC90_individual_model_mtau)]]
  
  # Multi-bandwidth BIC at tau=0.5
  min_length <- round(min(parameter_list[,1])*n)
  # Multi-bandwidth BIC at tau=(0.1,0.5,0.9)
  estcp_multiband_BIC90_mtau <- multiband_BIC(candidates=sum_90, y=y, min_length=min_length, tau_sets=tau_sets, BIC_factor=BIC_factor,
                                              max_num_cand=max_num_cand)

  # Summary
  # Summarize result
  list(BIC_factor=BIC_factor, tau_sets=tau_sets,
       # multiple tau at tau=(0.1,0.5,0.9)
       BIC90_individual_model_mtau=BIC90_individual_model_mtau, 
       estcp_BIC90_mtau=estcp_BIC90_mtau, 
       estcp_multiband_BIC90_mtau=estcp_multiband_BIC90_mtau)
}

########## GOALS with noncrossing constraint; Fig1 & Fig S2
contcolor <- c("black","red","#3E9F05")
shortlist <- c("United States","United Kingdom","India","Brazil","France","Russia","Spain","Argentina","South Africa","Australia")
for(country_index in 1:length(shortlist)){
  country_name <- shortlist[country_index]
  country_data <- subset(full_data, (location==country_name)&(total_cases>1000)&(new_cases>=1))
  last_day_reference_point <- max(which(as.Date(country_data$date)<=as.Date(last_day))) # some country does not have last day data
  y <- log(country_data$new_cases)[1:last_day_reference_point]
  y[is.na(y)] <- 0
  n <- length(y)
  date <- format(as.Date(country_data$date),format="%b-%d")[1:n]
  print(country_name)
  print(c(date[1], date[n]))
  print(n)
  
  # plot estimated change-points (grid, trim) = (0.1, 0.02)
  plot(y, xaxt="n", xlab="",ylab="", main=country_name, cex.main=2, cex.lab=1.5, ylim=c(min(y)-1, max(y)+1))
  mtext(TeX("Y_t"), side = 2, line = 2,cex=1.5)
  axis(1,c(1,floor(n/5)*(1:4),n),date[c(1,floor(n/5)*(1:4),n)])
  fit <- list()
  sweep_res <- country_result[[which(country_list==country_name)]]$result[[5]]$SN_sweep_result # tuning parameter=(0.1,0.02)
  sweep <- sweep_res
  critical_sets <- unlist(subset(critical_table, epsilon==0.1&delta==0.02))[3:5]
  x1 <- SN_local_peak(grid_size=floor(0.1*n), max_sweep_result=max_SNsweep(sweep), critical_value=critical_sets[1])
  x1_localmax <- SN_local_peak(grid_size=floor(0.1*n), max_sweep_result=max_SNsweep(sweep), critical_value=0)
  
  
  
  cpt <- c(0,x1,n)
  for (i in 2:(length(cpt)-1)){
    buffer <- 2
    x <- ((cpt[i]+1):(cpt[i+1]-buffer))/n
    reg_nocross <- rq.no.cross_insample(y[(cpt[i]+1):(cpt[i+1]-buffer)], x=cbind(x), taus=tau_sets)
    fix <- ((cpt[i]+1):(cpt[i+1]-buffer))
    for (j in 1:length(tau_sets)){
      lines(fix, cbind(1,x)%*%reg_nocross$bhat[,j], col=contcolor[j], lwd=2)
    }
    text(cpt[i]+10,max(y)+0.7,round(reg_nocross$bhat[2,2]/n,3),cex=1.15)
  }
  mylims <- par("usr")
  
  # plot estimated change-points using multi-scanning
  x1_bic <- country_result_BIC[[which(country_list==country_name)]]$estcp_multiband_BIC90_mtau$est_cp
  points(x1_bic, rep(mylims[4],length(x1_bic)), col='black', pch=25,  cex=2.5, bg='black')
  if (length(cpt)>=3){
    for (i in 1:(length(cpt)-1)){
      # points(cpt[i],y[cpt[i]],pch=16,cex=1.5)
      # text(cpt[i]+8,y[cpt[i]],date[cpt[i]],cex=1.2)
      abline(v=cpt[i])
      text(cpt[i]+10,min(y)-1,date[cpt[i]+1],cex=1)
    }
  }
  
  x1_localmax_ignored <- x1_localmax[!x1_localmax%in%x1]
  for(tmp_x1 in x1_localmax_ignored){
    abline(v=tmp_x1, lty=2)
    text(tmp_x1+10,min(y)-1,date[tmp_x1],cex=1)
  }
  
  #Plot the SN test statistics   ##Fig S2
  plot(max_SNsweep(sweep), xaxt="n", xlab="", ylab="", main=country_name, cex.main=2)
  mtext(TeX("$T_{n,\\epsilon,\\delta}(k)$"), side = 2, line = 2,cex=1.5)
  axis(1,c(1,floor(n/5)*(1:4),n),date[c(1,floor(n/5)*(1:4),n)])
  abline(v=x1_localmax)
  abline(h=critical_sets[1], col='blue', lty=2)
}
# dev.off()

######### Clustering analysis; Fig2 & Fig S5
selected_country_list <- country_list
selected_country_result <- country_result[which(country_list%in%selected_country_list)]
selected_num_country <- length(selected_country_list)
cp_results <- list()
for(country_index in 1:selected_num_country){
  country_name <- selected_country_list[country_index]
  country_data <- subset(full_data, (location==country_name)&(total_cases>1000)&(new_cases>=1))
  continent_index <- country_data$continent[1]
  country_abb <- country_data$iso_code[1]
  last_day_reference_point <- max(which(as.Date(country_data$date)<=as.Date(last_day))) # some country does not have last day data
  y <- log(country_data$new_cases)[1:last_day_reference_point]
  y[is.na(y)] <- 0
  n <- length(y)
  date <- format(as.Date(country_data$date),format="%b-%d")[1:n]
  
  fit <- list()
  sweep_res <- country_result[[which(country_list==country_name)]]$result[[5]]$SN_sweep_result
  sweep <- sweep_res
  critical_sets <- unlist(subset(critical_table, epsilon==0.1&delta==0.02))[3:5]
  x1 <- SN_local_peak(grid_size=floor(0.1*n), max_sweep_result=max_SNsweep(sweep), critical_value=critical_sets[1])
  # x1 <- country_result_BIC[[which(country_list==country_name)]]$estcp_multiband_BIC90_mtau$est_cp
  
  cpt <- c(0,x1,n)
  tmp_slope10 <- tmp_slope50 <- tmp_slope90 <- c()
  for (i in 1:(length(cpt)-1)){
    x <- ((cpt[i]+1):(cpt[i+1]))/n
    reg_nocross <- rq.no.cross_insample(y[(cpt[i]+1):(cpt[i+1])], x=cbind(x), taus=tau_sets)
    reg_nocross_slope <- rq.no.cross_insample(y[(cpt[i]+1):(cpt[i+1])], x=cbind(x*n), taus=tau_sets)
    tmp_slope10 <- c(tmp_slope10, reg_nocross_slope$bhat[2,which(tau_sets==0.1)])
    tmp_slope50 <- c(tmp_slope50, reg_nocross_slope$bhat[2,which(tau_sets==0.5)])
    tmp_slope90 <- c(tmp_slope90, reg_nocross_slope$bhat[2,which(tau_sets==0.9)])
  }
  # Detected change-point for each country
  cp_results[[country_index]] <- list(country_name=country_name, y=y, date=date, cpt=cpt, slope10=rep(tmp_slope10, diff(cpt)),
                                      slope50=rep(tmp_slope50, diff(cpt)), slope90=rep(tmp_slope90, diff(cpt)),continent=continent_index,code=country_abb)
}

### Slope analysis
mean_slope <- last_slope <- max_slope <- median_slope <- c()
continent_index <- c()
country_abb <- c()
for(country_index in 1:selected_num_country){
  mean_slope <- c(mean_slope, mean(cp_results[[country_index]]$slope50))
  median_slope <- c(median_slope, median(cp_results[[country_index]]$slope50))
  last_slope <- c(last_slope, cp_results[[country_index]]$slope50[length(cp_results[[country_index]]$slope50)])
  max_slope <- c(max_slope, max(cp_results[[country_index]]$slope50))
  continent_index <- c(continent_index, as.character(cp_results[[country_index]]$continent))
  country_abb <- c(country_abb, as.character(cp_results[[country_index]]$code))
}

continent_index[continent_index%in%c("Oceania","Asia","Africa")] <- "Asia+Others"
continent_index[continent_index%in%c("South America")] <- "Latin America"
continent_index[10] <- "Latin America"
continent_index[1] <- "USA+CAN"
continent_index[32] <- "USA+CAN"
continent_index <- factor(continent_index,levels=c("Europe","USA+CAN","Latin America","Asia+Others"))
slopes <- data.frame(last_slope, max_slope, median_slope, continent_index)
boxplot(last_slope~continent_index, data=slopes, main="Boxplot of Current Slopes", ylab="", xlab=NULL, cex.main=2, cex.lab=1.5)

### Time series clustering
# calculation of the correlation matrix
cor_matrix <- matrix(1, selected_num_country, selected_num_country)
colnames(cor_matrix) <- rownames(cor_matrix) <- selected_country_list
for(country_index1 in 1:(selected_num_country-1)){
  for(country_index2 in (country_index1+1):selected_num_country){
    common_date1 <- cp_results[[country_index1]]$date%in%intersect(cp_results[[country_index1]]$date, cp_results[[country_index2]]$date)
    common_date2 <- cp_results[[country_index2]]$date%in%intersect(cp_results[[country_index1]]$date, cp_results[[country_index2]]$date)
    # tmp_slope1 <- cp_results[[country_index1]]$slope10[common_date1]
    # tmp_slope2 <- cp_results[[country_index2]]$slope10[common_date2]
    tmp_slope1 <- cp_results[[country_index1]]$slope50[common_date1]
    tmp_slope2 <- cp_results[[country_index2]]$slope50[common_date2]
    # tmp_slope1 <- cp_results[[country_index1]]$slope90[common_date1]
    # tmp_slope2 <- cp_results[[country_index2]]$slope90[common_date2]
    tmp_cor <- cor(tmp_slope1, tmp_slope2)
    cor_matrix[country_index1, country_index2] <- cor_matrix[country_index2, country_index1] <- tmp_cor
  }
}
cor_matrix[is.na(cor_matrix)] <- -1


### Use the correlation matrix as a distance matrix
dist_matrix <- as.dist(1-cor_matrix)
contcolor <- c("black","red","#0000FF","#3E9F05")
contpch <- c(2,3,4,19)
continent_num <- factor(continent_index)
levels(continent_num) <- c("1","2","3","4")
continent_num <- as.numeric(as.character(continent_num))

### MDS visulization
fit <- sammon(dist_matrix, k=2)

x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="", ylab="",type="n", main=TeX("\\textbf{Sammon MDS}, $\\tau$=0.5"),ylim=c(-0.85,1.05),cex.main=2, cex.lab=1.5)
mtext("Coordinate 2", side = 2, line = 2,cex=1.5)
mtext("Coordinate 1", side = 1, line = 2,cex=1.5)

for (i in 1:length(x)){
  points(x[i],y[i],pch=contpch[continent_num[i]],col=contcolor[continent_num[i]],cex=1.4)
  text(x[i],y[i]-0.05, labels=country_abb[i],col=contcolor[continent_num[i]], cex=1.2)
}
legend("topright",cex=1,c("Europe", "USA+CAN","Latin America","Asia+Others"),col=contcolor,pch=contpch)

