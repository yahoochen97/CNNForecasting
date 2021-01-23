library(dplyr)
# define horizons in the loyo process
horizons = c('0',
             '7',
             '14',
             '21',
             '28',
             '42',
             '56')

# the optimal index of hyperparameters in the loyo process
# optimal index can be obtained with loocv_nlZs.R
best_gp_idx = read.csv(paste("results/", "GP", "_opthyp.csv", sep=''))
best_gp_idx = best_gp_idx$opt_idx

best_lm_idx = read.csv(paste("results/", "LM", "_opthyp.csv", sep=''))
best_lm_idx = best_lm_idx$opt_idx

# PVALUES = c()
# 
# gp_avg_nlz = c()
# lm_avg_nlz = c()
# 
# for (a in 1:length(horizons)) {
#   
#   # load the prior files
#   gp_nlz = paste('results/stan_NLZ', 'GP', '_' , test_year, 'day', horizons[a], '_', best_gp_idx[a] ,'.csv',sep='')
#   lm_nlz = paste('results/stan_NLZ', 'LM', '_' , test_year, 'day', horizons[a], '_', best_lm_idx[a] ,'.csv',sep='')
#   gp_data <- read.csv(gp_nlz)
#   lm_data <- read.csv(lm_nlz)
#   
#   gp_data = gp_data$nlz
#   lm_data = lm_data$nlz
#   
#   print(strtoi(horizons[a]))
#   TEST = t.test(gp_data, lm_data, paired = TRUE)
#   PVALUES = c(PVALUES, TEST$p.value)
#   
#   gp_avg_nlz = c(gp_avg_nlz, (gp_data))
#   lm_avg_nlz = c(lm_avg_nlz, (lm_data))
# }
# 
# results <- data.frame(horizons, PVALUES)
# write.csv(results, "results/pvalues.csv")


generate_DR_results <- function(TYPE, best_idx, test_years, horizons){
  RMSE = c() 
  ACCURACY = c()
  ENTROPY = c()
  RATIO = c()
  LL = c()
  CORRELATION = c()
  
  for (a in length(horizons):1) {
    RMSE_Y = c() 
    ACCURACY_Y = c()
    ENTROPY_Y = c()
    RATIO_Y = c()
    NLZ_Y = c()
    CORRELATION_Y = c()
    for(test_year in test_years){
      # load result files
      data = paste('results/stan_LOO', TYPE, '_' , test_year, 'day', horizons[a], '_', best_idx[a] ,'.csv',sep='')
      data <- read.csv(data)
      CORRELATION_Y = c(CORRELATION_Y, cor(data$vote, data$median))
      RMSE_Y = c(RMSE_Y, mean(data$rmse))
      
      correct_prediction = 0
      print(horizons[a])
      for(state in unique(data$state)){
        tmp = data[data$cycle==test_year & data$state==state,c("win","vote")]
        if( which(max(tmp$win)==tmp$win)== which(max(tmp$vote)==tmp$vote)){
          correct_prediction = correct_prediction + 1
        }
        else{
          print(test_year)
          print(state)
        }
      }
      
      ACCURACY_Y = c(ACCURACY_Y, correct_prediction/length(unique(data$state)))
      ENTROPY_Y = c(ENTROPY_Y, mean(log(data$win+1e-10)*data$winners))
      RATIO_Y = c(RATIO_Y, sum((data$vote>data$lower95) & (data$vote<data$upper95))/nrow(data))
      
      nlz = paste('results/stan_NLZ', TYPE, '_' , test_year, 'day', horizons[a], '_', best_idx[a] ,'.csv',sep='')
      nlz <- read.csv(nlz)
      nlz = nlz$nlz
      NLZ_Y = c(NLZ_Y, mean(nlz))
    }
    RMSE = c(RMSE, mean(RMSE_Y)) 
    ACCURACY = c(ACCURACY, mean(ACCURACY_Y))
    ENTROPY = c(ENTROPY, mean(ENTROPY_Y))
    RATIO = c(RATIO, mean((RATIO_Y)))
    LL = c(LL, -mean(NLZ_Y))
    CORRELATION = c(CORRELATION, mean(CORRELATION_Y))
  }
  
  results <- data.frame(CORRELATION,
                        RMSE,
                        ACCURACY,
                        ENTROPY,
                        RATIO,
                        LL)
  
  results = results %>% mutate_if(is.numeric, round, digits = 4)
  if(length(test_years)==1){
    test_year = test_years[1]
  }
  else{
    test_year = "92-16"
  }
  results = t(results)
  colnames(results) = rev(horizons)
  write.csv((results), paste("results/DR_", TYPE,test_year,"results.csv",sep=""))
}



generate_Prior_results <- function(TYPE, best_idx, test_years, horizons){
  RMSE = c() 
  ACCURACY = c()
  ENTROPY = c()
  RATIO = c()
  LL = c()
  CORRELATION = c()
  
  for (a in length(horizons):1) {
    RMSE_Y = c() 
    ACCURACY_Y = c()
    ENTROPY_Y = c()
    RATIO_Y = c()
    NLZ_Y = c()
    CORRELATION_Y = c()
    for(test_year in test_years){
      # load result files
      data = paste('results/LOO', TYPE, '_' , test_year, 'day', horizons[a], '_', best_idx[a] ,'.csv',sep='')
      data <- read.csv(data)
      data <- data[data$cycle==test_year,]
      data$vote = data$vote / 100
      CORRELATION_Y = c(CORRELATION_Y, cor(data$vote, data$posteriormean))
      # prior is Gaussian, close form is available
      RMSE_Y = c(RMSE_Y, sqrt(mean((data$posteriorstd)^2+(data$posteriormean-data$vote)^2)))
      
      winning_probs = c()
      winners = c()
      nlz = c()
      correct_prediction = 0
      
      for(state in unique(data$state)){
        tmp = data[data$cycle==test_year & data$state==state,]
        n = nrow(tmp)
        m = 10000
        samples = c()
        for(k in 1:n){
          samples = c(samples, (rnorm(m, mean=tmp$posteriormean[k],sd = tmp$posteriorstd[k])))
          nlz = c(nlz, -log(dnorm(tmp$vote[k],mean=tmp$posteriormean[k],sd = tmp$posteriorstd[k])))
        }
        samples <- matrix(samples, nrow = n, byrow = TRUE)
        win_rates = rep(0, n)
        for(k in 1:ncol(samples)){
          idx = which.max(samples[,k])
          win_rates[idx] = win_rates[idx] + 1
        }
        win_rates = win_rates / sum(win_rates)
        winning_probs = c(winning_probs, win_rates)
        winner = rep(0, n)
        winner[which(max(tmp$vote)==tmp$vote)] = 1
        winners = c(winners, winner)
        if( which(max(win_rates)==win_rates)== which(max(tmp$vote)==tmp$vote)){
          correct_prediction = correct_prediction + 1
        }
        
      }
      
      ACCURACY_Y = c(ACCURACY_Y, correct_prediction/length(unique(data$state)))
      ENTROPY_Y = c(ENTROPY_Y, mean(log(winning_probs+1e-10)*winners))
      l = data$posteriormean - 1.96*data$posteriorstd
      u = data$posteriormean + 1.96*data$posteriorstd
      RATIO_Y = c(RATIO_Y, sum((data$vote>l) & (data$vote<u))/nrow(data))
      NLZ_Y = c(NLZ_Y, mean(nlz))
    }
    RMSE = c(RMSE, mean(RMSE_Y)) 
    ACCURACY = c(ACCURACY, mean(ACCURACY_Y))
    ENTROPY = c(ENTROPY, mean(ENTROPY_Y))
    RATIO = c(RATIO, mean((RATIO_Y)))
    LL = c(LL, -mean(NLZ_Y))
    CORRELATION = c(CORRELATION, mean(CORRELATION_Y))
  }
  
  results <- data.frame(CORRELATION,
                        RMSE,
                        ACCURACY,
                        ENTROPY,
                        RATIO,
                        LL)
  
  results = results %>% mutate_if(is.numeric, round, digits = 4)
  if(length(test_years)==1){
    test_year = test_years[1]
  }
  else{
    test_year = "92-16"
  }
  results = t(results)
  colnames(results) = rev(horizons)
  write.csv((results), paste("results/Prior_", TYPE,test_year,"results.csv",sep=""))
}


generate_BRW_results <- function(test_years, horizons){
  RMSE = c() 
  ACCURACY = c()
  ENTROPY = c()
  RATIO = c()
  LL = c()
  CORRELATION = c()
  
  for (a in length(horizons):1) {
    RMSE_Y = c() 
    ACCURACY_Y = c()
    ENTROPY_Y = c()
    RATIO_Y = c()
    NLZ_Y = c()
    CORRELATION_Y = c()
    for(test_year in test_years){
      # load result files
      data = paste('results/brw_' , test_year, 'day_', horizons[a], '.csv',sep='')
      data <- read.csv(data)
      CORRELATION_Y = c(CORRELATION_Y, cor(data$vote, data$pmean))
      RMSE_Y = c(RMSE_Y, data$rmse)
      
      winning_probs = c()
      winners = c()
      correct_prediction = 0
      
      for(state in unique(data$state)){
        tmp = data[data$cycle==test_year & data$state==state,]
        n = nrow(tmp)
        winning_probs = c(winning_probs, tmp$win)
        winner = rep(0, n)
        winner[which(max(tmp$vote)==tmp$vote)] = 1
        winners = c(winners, winner)
        if( which(max(tmp$win)==tmp$win)==which(max(tmp$vote)==tmp$vote)){
          correct_prediction = correct_prediction + 1
        }
      }
      
      ACCURACY_Y = c(ACCURACY_Y, correct_prediction/length(unique(data$state)))
      ENTROPY_Y = c(ENTROPY_Y, mean(log(winning_probs+1e-10)*winners))
      RATIO_Y = c(RATIO_Y, sum((data$vote>data$lower95) & (data$vote<data$upper95))/nrow(data))
      nlz = data$nlz
      nlz = nlz[nlz<=1000]
      # nlz = -dnorm(data$vote, data$pmean, (data$upper95-data$lower95)/4, log=TRUE)
      NLZ_Y = c(NLZ_Y, mean(nlz, na.rm = TRUE))
    }
    RMSE = c(RMSE, mean(RMSE_Y)) 
    ACCURACY = c(ACCURACY, mean(ACCURACY_Y))
    ENTROPY = c(ENTROPY, mean(ENTROPY_Y))
    RATIO = c(RATIO, mean((RATIO_Y)))
    LL = c(LL, -mean(NLZ_Y))
    CORRELATION = c(CORRELATION, mean(CORRELATION_Y))
  }
  
  results <- data.frame(CORRELATION,
                        RMSE,
                        ACCURACY,
                        ENTROPY,
                        RATIO,
                        LL)
  
  results = results %>% mutate_if(is.numeric, round, digits = 3)
  if(length(test_years)==1){
    test_year = test_years[1]
  }
  else{
    test_year = "92-16"
  }
  results = t(results)
  colnames(results) = rev(horizons)
  write.csv((results), paste("results/brw_",test_year,"results.csv",sep=""))
}


test_years = c(1992,1994,1996,1998,2000,2002,2004,2006,2008,2010,2012,2014,2016)

# generate_DR_results('GP',best_gp_idx, c(2018), horizons)

# generate_DR_results('LM',best_lm_idx, c(2018), horizons)

# generate_DR_results('GP',best_gp_idx, test_years , horizons)

# generate_DR_results('LM',best_lm_idx, test_years, horizons)

# generate_Prior_results('GP',best_gp_idx, c(2018), horizons)

# generate_Prior_results('LM',best_lm_idx, c(2018), horizons)

# generate_Prior_results('GP',best_gp_idx, test_years, horizons)
 
# generate_Prior_results('LM',best_lm_idx, test_years, horizons)

# generate_BRW_results(test_years, horizons)

# generate_BRW_results(c(2018), horizons)
