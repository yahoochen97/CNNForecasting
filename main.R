# setwd('/Users/yahoo/Documents/WashU/CSE515T/Code/Gaussian Process/')

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, use default 2020 with GP model
if (length(args)==0) {
  test_year = 2020
  # define model type: gp prior or lm prior
  TYPE = 'GP'
}
if (length(args)==1){
  # define the test year
  test_year = as.double(args[1])
  TYPE = 'GP'
}
if (length(args)==2){
  # define the test year
  test_year = as.double(args[1])
  TYPE = args[2]
}

# load packages
library(rstan)
# library(bayesplot)
# library(MCMCpack)
# library(dplyr)
# library(ggridges)
# library(ggplot2)
# library(grid)

PLOT = FALSE

# define horizons
horizons = c('0',
               '7',
               '14',
               '21',
               '28',
               '42',
               '56')

# the optimal index of hyperparameters in the loyo process
# optimal index can be obtained with loocv_nlZs.R

# best_cv_idx = read.csv(paste("results/", TYPE, "_opthyp.csv", sep=''));
# best_cv_idx = best_cv_idx$opt_idx

# R&R
horizon = args[1]
cv_year = args[2]
TYPE = args[3]

best_cv_idx = read.csv(paste("results/GP_opthyp.csv", sep=''))
IDX = best_cv_idx$opt_idx[best_cv_idx$horizons==str2lang(horizon)]

test_year=2016
a = 1

# 1:length(horizons)
for (b in (IDX):(IDX)) {
  
  # load the prior files
  # input_file = paste('results/LOO', TYPE, '_' , test_year, 'day', horizons[a], '_', best_cv_idx[a] ,'.csv',sep='')
  # output_file = paste('results/stan_LOO', TYPE, '_' , test_year, 'day', horizons[a], '_', best_cv_idx[a] ,'.csv',sep='')
  input_file = paste('results/GP_', TYPE, '_' , cv_year, 'day', horizon, '_', b ,'.csv',sep='')
  output_file = paste('nlZs/stan_GP_', TYPE, '_' , cv_year, 'day', horizon, '_', b,'.csv',sep='')
  
  data <- read.csv(input_file)
  print(input_file)
  
  # remove unlike candidates of races with >4 candidates
  data <- data[data$cycle!=2016 | data$state!='Louisiana' | data$candidate!='Flemsing',]
  data <- data[data$cycle!=2020 | data$state!='Georgia' | data$candidate!='Loeffler',]
  data <- data[data$cycle!=2020 | data$state!='Georgia' | data$candidate!='Tarver',]
  
  # data %>%
  #   group_by(cycle, state) %>%
  #   summarise(count=n()) %>%
  #   filter(count >=4)
  
  # split training and testing data
  data_test <- data[(data$cycle==test_year),]
  data <- data[(data$cycle!=test_year & data$cycle!=2018 & data$cycle!=2020),]

  cycles <- unique(data$cycle)
  states <- union(unique(data$state), unique(data_test$state))

  # maximal number of candidates allowed in the model  
  C <- 4
  
  # define meta variables
  metadata <- list()
  year_idx <- c()
  mu <- list()
  sigma <- list()
  y <- list()
  nc <- c()
  pvi <- list()
  party <- list()
  experienced <- list()
  counter <- 0
  
  # iterate over races
  for (cycle in cycles) {
    for (state in states) {
      # obtain priors and fundamentals
      pmu = data[data$cycle==cycle & data$state==state,c("posteriormean")]
      pstd = data[data$cycle==cycle & data$state==state,c("posteriorstd")]
      vote = data[data$cycle==cycle & data$state==state,c("vote")]
      pvi_ = data[data$cycle==cycle & data$state==state,c("pvi")]
      party_ = data[data$cycle==cycle & data$state==state,c("party")]
      experienced_ = data[data$cycle==cycle & data$state==state,c("experienced")]
      
      # not every state has election in all election cycles
      # proceed only it has
      if(length(pmu)){
        counter <- counter + 1
        metadata[[counter]] = c(cycle, state)
        mu[[counter]] = pmu
        sigma[[counter]] = pstd
        y[[counter]] = vote / 100
        pvi[[counter]] = pvi_
        party[[counter]] = party_
        experienced[[counter]] = experienced_
        nc = c(nc, length(vote))
        year_idx = c(year_idx, (cycle-1990)/2)
      }
    }
  }
  
  # build stan data
  stan_mu <- matrix(0,counter,C)
  stan_sigma <- matrix(0.0001,counter,C)
  stan_y <- matrix(0,counter,C)
  stan_pvi <- matrix(0,counter,C)
  stan_party <- matrix(0,counter,C)
  stan_experienced <- matrix(0,counter,C)
  
  for (i in 1:counter) {
    stan_mu[i,1:nc[i]] <- mu[[i]]
    stan_sigma[i,1:nc[i]] = sigma[[i]]
    stan_y[i,1:nc[i]] = y[[i]]
    if(sum(stan_y[i,1:nc[i]])!=0){
      stan_y[i,1:nc[i]] = stan_y[i,1:nc[i]]/sum(stan_y[i,1:nc[i]])
    }
    else{
      stan_y[i,1:nc[i]] = 1/nc[i]
    }
    
    stan_pvi[i,1:nc[i]] = pvi[[i]]
    stan_party[i,1:nc[i]] = party[[i]]
    stan_experienced[i, 1:nc[i]] = experienced[[i]]
  }
  
  # split data into categories based on number of candidates
  idx2 <- c()
  idx3 <- c()
  idx4 <- c()
  
  for (i in 1:counter) {
    tmp = sum(stan_mu[i,]!=0)
    if(tmp==2) idx2 = c(idx2, i)
    if(tmp==3) idx3 = c(idx3, i)
    if(tmp==4) idx4 = c(idx4, i)
  }
  
  # test data
  test_metadata <- list()
  test_year_idx <- c()
  test_mu <- list()
  test_sigma <- list()
  test_y <- list()
  
  test_nc <- c()
  test_pvi <- list()
  test_party <- list()
  test_experienced <- list()
  test_counter <- 0
  
  # iterate over races
  for (cycle in unique(data_test$cycle)){
    for (state in states) {
      pmu = data_test[(data_test$state==state & data_test$cycle==cycle),c("posteriormean")]
      pstd = data_test[data_test$state==state & data_test$cycle==cycle,c("posteriorstd")]
      vote = data_test[data_test$state==state & data_test$cycle==cycle,c("vote")]
      pvi_ = data_test[data_test$state==state & data_test$cycle==cycle,c("pvi")]
      party_ = data_test[data_test$state==state & data_test$cycle==cycle,c("party")]
      experienced_ = data_test[data_test$state==state & data_test$cycle==cycle,c("experienced")]
      if(length(pmu)){
        test_counter = test_counter + 1
        test_metadata[[test_counter]] = c(cycle,state)
        test_mu[[test_counter]] = pmu
        test_sigma[[test_counter]] = pstd
        test_y[[test_counter]] = vote / 100
        test_pvi[[test_counter]] = pvi_
        test_party[[test_counter]] = party_
        test_experienced[[test_counter]] = experienced_
        test_nc = c(test_nc, length(vote))
        test_year_idx = c(test_year_idx, (cycle-1990)/2)
      }
    }
  }
  
  # build stan data
  test_stan_mu <- matrix(0,test_counter,C)
  test_stan_sigma <- matrix(0.0001,test_counter,C)
  test_stan_y <- matrix(0,test_counter,C)
  test_stan_pvi <- matrix(0,test_counter,C)
  test_stan_party <- matrix(0,test_counter,C)
  test_stan_experienced <- matrix(0,test_counter,C)
  
  for (i in 1:test_counter) {
    test_stan_mu[i,1:test_nc[i]] = test_mu[[i]]
    
    test_stan_sigma[i,1:test_nc[i]] = test_sigma[[i]]
    test_stan_y[i,1:test_nc[i]] = test_y[[i]]
    if(sum(test_stan_y[i,1:test_nc[i]])!=0){
      test_stan_y[i,1:test_nc[i]] = test_stan_y[i,1:test_nc[i]]/sum(test_stan_y[i,1:test_nc[i]])
    }
    else{
      test_stan_y[i,1:test_nc[i]] = 1/test_nc[i]
    }
    test_stan_pvi[i,1:test_nc[i]] = test_pvi[[i]]
    test_stan_party[i,1:test_nc[i]] = test_party[[i]]
    test_stan_experienced[i, 1:test_nc[i]] = test_experienced[[i]]
  }
  
  
  test_idx2 <- c()
  test_idx3 <- c()
  test_idx4 <- c()
  
  for (i in 1:test_counter) {
    tmp = sum(test_stan_mu[i,]!=0)
    if(tmp==2) test_idx2 = c(test_idx2, i)
    if(tmp==3) test_idx3 = c(test_idx3, i)
    if(tmp==4) test_idx4 = c(test_idx4, i)
  }
  
  # define stan data structure
  stan_data <- list(N2 = length(idx2), 
                    mu2 = stan_mu[idx2,1:2], 
                    sigma2 = stan_sigma[idx2,1:2],
                    y2 = stan_y[idx2,1:2],
                    pvi2 = stan_pvi[idx2,1:2],
                    party2 = stan_party[idx2,1:2],
                    experienced2 = stan_experienced[idx2,1:2],
                    year_idx2 = year_idx[idx2],
                    N3 = length(idx3), 
                    mu3 = stan_mu[idx3,1:3], 
                    sigma3 = stan_sigma[idx3,1:3],
                    y3 = stan_y[idx3,1:3],
                    pvi3 = stan_pvi[idx3,1:3],
                    party3 = stan_party[idx3,1:3],
                    experienced3 = stan_experienced[idx3,1:3],
                    year_idx3 = year_idx[idx3],
                    N4 = length(idx4), 
                    mu4 = matrix(stan_mu[idx4,1:4],ncol=4,byrow = FALSE), 
                    sigma4 = matrix(stan_sigma[idx4,1:4],ncol=4,byrow = FALSE),
                    y4 = matrix(stan_y[idx4,1:4],ncol=4,byrow = FALSE),
                    pvi4 = matrix(stan_pvi[idx4,1:4],ncol=4,byrow = FALSE),
                    party4 = matrix(stan_party[idx4,1:4],ncol=4,byrow = FALSE),
                    experienced4 = matrix(stan_experienced[idx4,1:4],ncol=4,byrow = FALSE),
                    year_idx4 = array(year_idx[idx4]),
                    test_N2 = length(test_idx2), 
                    test_mu2 = test_stan_mu[test_idx2,1:2], 
                    test_sigma2 = test_stan_sigma[test_idx2,1:2],
                    test_pvi2 = test_stan_pvi[test_idx2,1:2],
                    test_party2 = test_stan_party[test_idx2,1:2],
                    test_experienced2 = test_stan_experienced[test_idx2,1:2],
                    test_year_idx2 = test_year_idx[test_idx2],
                    test_f2 = test_stan_y[test_idx2,1:2],
                    test_N3 = length(test_idx3), 
                    test_mu3 = matrix(test_stan_mu[test_idx3,1:3],ncol=3,byrow = FALSE), 
                    test_sigma3 = matrix(test_stan_sigma[test_idx3,1:3],ncol=3,byrow = FALSE),
                    test_pvi3 = matrix(test_stan_pvi[test_idx3,1:3],ncol=3,byrow = FALSE),
                    test_party3 = matrix(test_stan_party[test_idx3,1:3],ncol=3,byrow = FALSE),
                    test_experienced3 = matrix(test_stan_experienced[test_idx3,1:3],ncol=3,byrow = FALSE),
                    test_year_idx3 = array(test_year_idx[test_idx3]),
                    test_f3 = matrix(test_stan_y[test_idx3,1:3],ncol=3,byrow = FALSE),
                    test_N4 = length(test_idx4), 
                    test_mu4 = matrix(test_stan_mu[test_idx4,1:4],ncol=4,byrow = FALSE), 
                    test_sigma4 = matrix(test_stan_sigma[test_idx4,1:4],ncol=4,byrow = FALSE),
                    test_pvi4 = matrix(test_stan_pvi[test_idx4,1:4],ncol=4,byrow = FALSE),
                    test_party4 = matrix(test_stan_party[test_idx4,1:4],ncol=4,byrow = FALSE),
                    test_experienced4 = matrix(test_stan_experienced[test_idx4,1:4],ncol=4,byrow = FALSE),
                    test_year_idx4 = array(test_year_idx[test_idx4]),
                    test_f4 = matrix(test_stan_y[test_idx4,1:4],ncol=4,byrow = FALSE),
                    max_year_idx = max(c(year_idx, test_year_idx)))
  
  # define stan model
  # model <- stan_model("model.stan")
  
  # train stan model
  fit <- stan(file = "model.stan",
              data = stan_data, 
              warmup = 1000, 
              iter = 5000, 
              chains = 3, 
              cores = 3, 
              thin = 4,
              control=list(adapt_delta=.99, max_treedepth = 15),
              seed = a,
              refresh=0
  )
  # saveRDS(fit, file = paste("models/",TYPE, "_", test_year, "day_", horizons[a] ,"_fit.rds",sep=''))
  
  saveRDS(fit, file = paste("models/GP_",TYPE, "_", test_year, "day_", horizon ,"_fit.rds",sep=''))
  
  # print(summary(fit,c('alpha','beta','ppb','eb','year_sig'))$summary)
  fit_params <- as.data.frame(fit)
  
  CYCLE <- c()
  STATE <- c()
  CANDIDATE <- c()
  POSTERIORMEAN <- c()
  POSTERIORSTD <- c()
  PMEAN <- c()
  PSTD <- c()
  VOTE <- c()
  NORM_VOTE <- c()
  RMSE = c()
  LOWER95 <- c()
  UPPER95 <- c()
  WIN <- c()
  WINNERS <- c()
  MEDIAN <- c()
  NLZ <- c()
  PARTY <- c()
  
  correct_predictions <- 0
  Nout_test <- 0
  Nout <- 0
  
  COLNAMES = c('Posterior_Vote','Party','State','Type','VOTE')
  posteriors = data.frame(matrix(ncol = length(COLNAMES), nrow = 0))
  colnames(posteriors) = COLNAMES

  for(i in 1:length(test_idx2)) {
    cycle = test_metadata[[test_idx2[i]]][1]
    state = test_metadata[[test_idx2[i]]][2]
    pmu = data_test[data_test$state==state & data_test$cycle==cycle,c("posteriormean")]
    pstd = data_test[data_test$state==state & data_test$cycle==cycle,c("posteriorstd")]
    vote = test_stan_y[test_idx2[i],1:test_nc[test_idx2[i]]]
    party = data_test[data_test$state==state & data_test$cycle==cycle,c("party")]
    candidates = data_test[data_test$state==state & data_test$cycle==cycle,c("candidate")]
    preds= c()
    for(j in 1:2){
      tmp = paste('test_y2[',i,',',j,']',sep='')
      pred = fit_params[[tmp]]
      preds = c(preds, pred)
      u=quantile(pred,probs=c(0.975),names = FALSE)
      l=quantile(pred,probs=c(0.025),names = FALSE)
      m = mean(pred)
      s = sd(pred)
      
      one_posterior = data.frame(matrix(ncol = length(COLNAMES), nrow = length(pred)))
      colnames(one_posterior) = COLNAMES
      
      one_posterior$Posterior_Vote = pred*100
      one_posterior$VOTE = vote[j]*100
      if (state=='MinnesotaS'){
        one_posterior$State='MNS'
      }
      else if(state=='MississippiS'){
        one_posterior$State='MSS'
      }
      else{
        one_posterior$State = state.abb[match(state,state.name)]
      }
      if(party[j]==-1){
        one_posterior$Party = 'REP'
        RepWin = mean(one_posterior$Posterior_Vote>50)
      }
      else{
        one_posterior$Party = 'DEM'
        RepWin = mean(one_posterior$Posterior_Vote<50)
      }
      
      # if(state=="California" & j==2){
      #   one_posterior$Party = 'REP'
      #   RepWin = mean(one_posterior$Posterior_Vote>50)
      # }
      
      if(RepWin>0.95){
        one_posterior$Type = "Safe R"
      }
      else if (RepWin>0.7) {
        one_posterior$Type = "Likely R"
      }
      else if (RepWin>0.55) {
        one_posterior$Type = "Lean R"
      }
      else if (RepWin>0.45) {
        one_posterior$Type = "Toss-up"
      }
      else if (RepWin>0.3) {
        one_posterior$Type = "Lean D"
      }
      else if (RepWin>0.05) {
        one_posterior$Type = "Likely D"
      }
      else {
        one_posterior$Type = "Safe D"
      }

      posteriors = rbind(posteriors, one_posterior)
      
      if (test_stan_y[test_idx2[i],j]>u | test_stan_y[test_idx2[i],j]<l){
        Nout_test = Nout_test + 1
      }

      CYCLE <- c(CYCLE, cycle)
      STATE <- c(STATE,state)
      CANDIDATE <- c(CANDIDATE,as.character(candidates[j]))
      POSTERIORMEAN <- c(POSTERIORMEAN,pmu[j])
      POSTERIORSTD <- c(POSTERIORSTD,pstd[j])
      PMEAN <- c(PMEAN, m)
      PSTD <- c(PSTD, s)
      RMSE = c(RMSE, sqrt(mean((pred - rep(vote[j],length(pred)))^2)))
      VOTE <- c(VOTE, vote[j])
      MEDIAN <- c(MEDIAN, median(pred))
      LOWER95 <- c(LOWER95, l)
      UPPER95 <- c(UPPER95, u)
      PARTY <- c(PARTY, party[j])
    }
    NLZ <- c(NLZ, -log(mean(exp(fit_params[[paste('test_ll2[',i,']',sep='')]]))))
    preds <- matrix(preds, nrow = 2, byrow = TRUE)
    win_rates = rep(0, 2)
    for(k in 1:ncol(preds)){
      idx = which.max(preds[,k])
      win_rates[idx] = win_rates[idx] + 1
    }
    win_rates = win_rates / sum(win_rates)
    
    WIN <- c(WIN, win_rates)
    winners = rep(0, 2)
    winners[which.max(vote)] = 1
    WINNERS  <- c(WINNERS, winners)
    if (which.max(win_rates)==which.max(vote)){
      correct_predictions = correct_predictions + 1
    }
    else{
      print("Wrong prediction:")
      print(test_metadata[[test_idx2[i]]])
    }
  }

  if(length(test_idx3)){
    for(i in 1:length(test_idx3)) {
      cycle = test_metadata[[test_idx3[i]]][1]
      state = test_metadata[[test_idx3[i]]][2]
      pmu = data_test[data_test$state==state & data_test$cycle==cycle,c("posteriormean")]
      pstd = data_test[data_test$state==state & data_test$cycle==cycle,c("posteriorstd")]
      vote = test_stan_y[test_idx3[i],1:test_nc[test_idx3[i]]]
      party = data_test[data_test$state==state & data_test$cycle==cycle,c("party")]
      candidates = data_test[data_test$state==state & data_test$cycle==cycle,c("candidate")]
      # preds <- sample_posterior(stan_mu[i,], stan_sigma[i,], nc[i], gs=10, ds=1000, fit_params=fit_params)
      preds= c()
      for(j in 1:3){
        tmp = paste('test_y3[',i,',',j,']',sep='')
        pred = fit_params[[tmp]]
        preds = c(preds, pred)
        u=quantile(pred,probs=c(0.975),names = FALSE)
        l=quantile(pred,probs=c(0.025),names = FALSE)
        m = mean(pred)
        s = sd(pred)
        
        # if (state=='Maine' & j==3){
        #   next
        # }
        
        if(party[j]!=0){
          one_posterior = data.frame(matrix(ncol = length(COLNAMES), nrow = length(pred)))
          colnames(one_posterior) = COLNAMES
          
          one_posterior$Posterior_Vote = pred*100
          one_posterior$VOTE = vote[j]*100
          if (state=='MinnesotaS'){
            one_posterior$State='MNS'
          }
          else if(state=='MississippiS'){
            one_posterior$State='MSS'
          }
          else{
            one_posterior$State = state.abb[match(state,state.name)]
          }
          if(party[j]==-1){
            one_posterior$Party = 'REP'
            RepWin = mean(one_posterior$Posterior_Vote>50)
          }
          else{
            one_posterior$Party = 'DEM'
            RepWin = mean(one_posterior$Posterior_Vote<50)
          }
          if(RepWin>0.95){
            one_posterior$Type = "Safe R"
          }
          else if (RepWin>0.7) {
            one_posterior$Type = "Likely R"
          }
          else if (RepWin>0.55) {
            one_posterior$Type = "Lean R"
          }
          else if (RepWin>0.45) {
            one_posterior$Type = "Toss-up"
          }
          else if (RepWin>0.3) {
            one_posterior$Type = "Lean D"
          }
          else if (RepWin>0.05) {
            one_posterior$Type = "Likely D"
          }
          else {
            one_posterior$Type = "Safe D"
          }
          posteriors = rbind(posteriors, one_posterior)
        }
        
        
        if (test_stan_y[test_idx3[i],j]>u | test_stan_y[test_idx3[i],j]<l){
          Nout_test = Nout_test + 1
        }
        CYCLE <- c(CYCLE, cycle)
        STATE <- c(STATE,state)
        CANDIDATE <- c(CANDIDATE,as.character(candidates[j]))
        POSTERIORMEAN <- c(POSTERIORMEAN,pmu[j])
        POSTERIORSTD <- c(POSTERIORSTD,pstd[j])
        PMEAN <- c(PMEAN, m)
        PSTD <- c(PSTD, s)
        RMSE = c(RMSE, sqrt(mean((pred - rep(vote[j],length(pred)))^2)))
        VOTE <- c(VOTE, vote[j])
        MEDIAN <- c(MEDIAN, median(pred))
        LOWER95 <- c(LOWER95, l)
        UPPER95 <- c(UPPER95, u)
        PARTY <- c(PARTY, party[j])
      }
      NLZ <- c(NLZ, -log(mean(exp(fit_params[[paste('test_ll3[',i,']',sep='')]]))))
      preds <- matrix(preds, nrow = 3, byrow = TRUE)
      win_rates = rep(0, 3)
      for(k in 1:ncol(preds)){
        idx = which.max(preds[,k])
        win_rates[idx] = win_rates[idx] + 1
      }
      win_rates = win_rates / sum(win_rates)
      WIN <- c(WIN, win_rates)
      winners = rep(0, 3)
      winners[which.max(vote)] = 1
      WINNERS  <- c(WINNERS, winners)
      if (which.max(win_rates)==which.max(vote)){
        correct_predictions = correct_predictions + 1
      }
      else{
        print("Wrong prediction:")
        print(test_metadata[[test_idx3[i]]])
      }
    }
  }

  # if (length(test_idx4)){
  #   for(i in 1:length(test_idx4)) {
  #     cycle = test_metadata[[test_idx4[i]]][1]
  #     state = test_metadata[[test_idx4[i]]][2]
  #     pmu = data_test[data_test$state==state & data_test$cycle==cycle ,c("posteriormean")]
  #     pstd = data_test[data_test$state==state & data_test$cycle==cycle ,c("posteriorstd")]
  #     vote = test_stan_y[test_idx4[i],1:test_nc[test_idx4[i]]]
  #     party = data_test[data_test$state==state & data_test$cycle==cycle,c("party")]
  #     candidates = data_test[data_test$state==state & data_test$cycle==cycle,c("candidate")]
  #     # preds <- sample_posterior(stan_mu[i,], stan_sigma[i,], nc[i], gs=10, ds=1000, fit_params=fit_params)
  #     preds= c()
  #     for(j in 1:4){
  #       tmp = paste('test_y4[',i,',',j,']',sep='')
  #       pred = fit_params[[tmp]]
  #       preds = c(preds, pred)
  #       u=quantile(pred,probs=c(0.975),names = FALSE)
  #       l=quantile(pred,probs=c(0.025),names = FALSE)
  #       m = mean(pred)
  #       s = sd(pred)
  #       if (test_stan_y[test_idx4[i],j]>u | test_stan_y[test_idx4[i],j]<l){
  #         Nout_test = Nout_test + 1
  #       }
  #       CYCLE <- c(CYCLE, cycle)
  #       STATE <- c(STATE, state)
  #       CANDIDATE <- c(CANDIDATE,as.character(candidates[j]))
  #       POSTERIORMEAN <- c(POSTERIORMEAN,pmu[j])
  #       POSTERIORSTD <- c(POSTERIORSTD,pstd[j])
  #       PMEAN <- c(PMEAN, m)
  #       PSTD <- c(PSTD, s)
  #       VOTE <- c(VOTE, vote[j])
  #       RMSE = c(RMSE, sqrt(mean((pred - rep(vote[j],length(pred)))^2)))
  #       MEDIAN <- c(MEDIAN, median(pred))
  #       LOWER95 <- c(LOWER95, l)
  #       UPPER95 <- c(UPPER95, u)
  #     }
  #     NLZ <- c(NLZ, -log(mean(exp(fit_params[[paste('test_ll4[',i,']',sep='')]]))))
  #     preds <- matrix(preds, nrow = 4, byrow = TRUE)
  #     win_rates = rep(0, 4)
  #     for(k in 1:ncol(preds)){
  #       idx = which.max(preds[,k])
  #       win_rates[idx] = win_rates[idx] + 1
  #     }
  #     win_rates = win_rates / sum(win_rates)
  #     WIN <- c(WIN, win_rates)
  #     winners = rep(0, 4)
  #     winners[which.max(vote)] = 1
  #     WINNERS  <- c(WINNERS, winners)
  #     if (which.max(win_rates)==which.max(vote)){
  #       correct_predictions = correct_predictions + 1
  #     }
  #     else{
  #       print("Wrong prediction:")
  #       print(test_metadata[[test_idx4[i]]])
  #     }
  #   }
  # }
  
  # LEVELS = posteriors[posteriors$Party=='DEM',] %>% 
  #   dplyr::group_by(State) %>% 
  #   dplyr::mutate(tmp = mean(Posterior_Vote)) %>% 
  #   dplyr::select(State, tmp) %>%
  #   dplyr::distinct(State, tmp) %>%
  #   dplyr::arrange((tmp)) %>%
  #   dplyr::select(State)
  # 
  # LEVELS = LEVELS$State
  # 
  # DEMVOTE2018 = c()
  # REPVOTE2018 = c()
  # DEMUPPER = c()
  # DEMLOWER = c()
  # DEMMEDIAN = c()
  # REPUPPER = c()
  # REPLOWER = c()
  # REPMEDIAN = c()
  # for (i in length(LEVELS):1){
  #   DEMVOTE2018 = c(DEMVOTE2018,posteriors[(posteriors$State)==LEVELS[i] & posteriors$Party=='DEM','VOTE'][1])
  #   REPVOTE2018 = c(REPVOTE2018,posteriors[(posteriors$State)==LEVELS[i] & posteriors$Party=='REP','VOTE'][1])
  #   pred = posteriors[(posteriors$State)==LEVELS[i] & posteriors$Party=='DEM','Posterior_Vote']
  #   u=quantile(pred,probs=c(0.975),names = FALSE)
  #   m=quantile(pred,probs=c(0.5),names = FALSE)
  #   l=quantile(pred,probs=c(0.025),names = FALSE)
  #   DEMUPPER = c(DEMUPPER, u)
  #   DEMLOWER = c(DEMLOWER, l)
  #   DEMMEDIAN = c(DEMMEDIAN, m)
  #   pred = posteriors[(posteriors$State)==LEVELS[i] & posteriors$Party=='REP','Posterior_Vote']
  #   u=quantile(pred,probs=c(0.975),names = FALSE)
  #   m=quantile(pred,probs=c(0.5),names = FALSE)
  #   l=quantile(pred,probs=c(0.025),names = FALSE)
  #   REPUPPER = c(REPUPPER, u)
  #   REPLOWER = c(REPLOWER, l)
  #   REPMEDIAN = c(REPMEDIAN, m)  
  # }
  # VOTE2018 = data.frame(c(rev(as.character(LEVELS)),rev(as.character(LEVELS))),
  #                       c(DEMVOTE2018, REPVOTE2018),
  #                       c(rep('DEM',length(DEMVOTE2018)),rep('REP',length(REPVOTE2018))),
  #                       c(DEMUPPER, REPUPPER),
  #                       c(DEMLOWER, REPLOWER),
  #                       c(DEMMEDIAN, REPMEDIAN)
  #                       )
  # colnames(VOTE2018) = c("state","vote","party", 'upper','lower','m')
  # if(PLOT){
  #   ggplot(VOTE2018, aes(fill=party, x=state, y=vote)) + 
  #     geom_bar(position="stack", stat="identity") +
  #     scale_x_discrete(limits = rev(LEVELS)) + 
  #     coord_flip() + 
  #     scale_fill_manual("Party", values = c("REP" = "red", "DEM" = "blue")) +
  #     ggtitle("Actual vote share for major party candidates") +
  #     theme(plot.title = element_text(hjust=0.5),
  #           panel.background = element_rect(fill = 'white', colour = 'white'))
  # }
  # 
  # posteriors$State <- factor(posteriors$State, levels = LEVELS)
  # 
  # LIKENAMES = c('Safe R', 'Likely R','Lean R', 'Toss-up', 'Lean D','Likely D','Safe D')
  # posteriors$Type <- factor(posteriors$Type, levels = LIKENAMES)
  # 
  # LEVELS = posteriors[posteriors$Party=='DEM',] %>% 
  #   dplyr::group_by(State) %>% 
  #   dplyr::mutate(tmp = mean(Posterior_Vote)) %>% 
  #   dplyr::select(State, tmp) %>%
  #   dplyr::distinct(State, tmp) %>%
  #   dplyr::arrange((tmp))
  # 
  # STATE_COLORS = c()
  # incumbents = read.csv("data/incumbent2020.csv")
  # for (i in nrow(LEVELS):1) {
  #   if(incumbents[as.character(incumbents$State)==LEVELS$State[i],'Incumbent']=='REP'){
  #     STATE_COLORS = c(STATE_COLORS, 'red')
  #   }
  #   else{
  #     STATE_COLORS = c(STATE_COLORS, 'blue')
  #   }
  # }
  # 
  # if(PLOT){
  #   p = ggplot(posteriors, aes(x = Posterior_Vote, y = reorder(State, desc(State)), color = Party, fill = Party)) +
  #     geom_density_ridges(alpha=0.6) +
  #     scale_y_discrete(expand = c(0, 0), name = "") +
  #     scale_x_continuous(expand = c(0, 0), breaks = c(0,20,40,60,80,100),
  #                        name = "Posterior Vote (%)") +
  #     theme(panel.grid.major.x = element_line(color = "gray"),
  #           # explicitly set the horizontal lines (or they will disappear too)
  #           panel.grid.major.y = element_line(size=.4, color="grey" ) ,
  #           axis.text.y = element_text(colour =STATE_COLORS)) +
  #     # geom_hline(yintercept=0:length(STATE_COLORS), linetype="solid", color = "grey", size=0.2) + 
  #     scale_fill_manual(values = c("blue","red"), labels = c("DEM","REP")) +
  #     scale_color_manual(values = c(NA,NA), guide = "none") +
  #     coord_cartesian(xlim = c(0, 100), clip='on') +
  #     guides(fill = guide_legend(
  #       override.aes = list(
  #         fill = c("blue","red"),
  #         color = NA, point_color = NA)
  #     )
  #     ) +
  #     geom_hline(yintercept=sum(LEVELS$tmp>=50)+1, linetype="dashed", color = "black") + 
  #     ggtitle("Posterior predictive density of vote share for major party candidates") +
  #     theme(plot.title = element_text(hjust=0.5),
  #           panel.background = element_rect(fill = 'white', colour = 'white'))
  #   
  # }
  # 
  # LEVELS = posteriors[posteriors$Party=='DEM',] %>% 
  #   dplyr::group_by(State) %>% 
  #   dplyr::mutate(tmp = mean(Posterior_Vote)) %>% 
  #   dplyr::select(State, tmp) %>%
  #   dplyr::distinct(State, tmp) %>%
  #   dplyr::arrange((tmp))
  # LEVELS = LEVELS$State
  # 
  # STATE_COLORS = c()
  # for (s in rev(sort(as.character(LEVELS)))) {
  #   if(s=="AZ" | s=="NV"){
  #     STATE_COLORS = c(STATE_COLORS, 'red')
  #   }
  #   else{
  #     STATE_COLORS = c(STATE_COLORS, 'black')
  #   }
  # }
  # 
  # LeonLine = data.frame(y1=VOTE2018[VOTE2018$state=="CA" & VOTE2018$party=='REP',"lower"],
  #                       y2=VOTE2018[VOTE2018$state=="CA" & VOTE2018$party=='REP',"upper"],
  #                       x1=length(STATE_COLORS)-1,
  #                       x2=length(STATE_COLORS)-1,
  #                       party="DEM")
  # if(PLOT){
  #   ggplot(VOTE2018, aes(y = m, x = reorder(state, desc(state)), fill = party)) + 
  #     geom_errorbar(size = 0.75, aes(ymin=lower, ymax=upper, color=party), width=0.1, position=position_dodge(width=0.5)) +
  #     geom_segment(data = LeonLine, aes(x = x1, y = y1, xend = x2, yend = y2, colour = "blue"),position=position_dodge(width=0.5)) +
  #     scale_color_manual(name="party", values = alpha(c("blue","blue","red"), 0.5),
  #                        labels=c("DEM","DEM","REP"))+
  #     geom_point(size=1, mapping = aes(x = state, y = m, fill=party, color=party),position=position_dodge(width=0.5)) + 
  #     scale_x_discrete(expand = c(0, 0), name = "") +
  #     scale_y_continuous(expand = c(0, 0),
  #                        name = "Posterior Vote (%)") +
  #     coord_flip() + 
  #     theme(panel.grid.major.x = element_line(color = "gray"),
  #           panel.grid.major.y = element_line(size=.4, color="grey" ) ,
  #           axis.text.y = element_text(colour =STATE_COLORS)) +
  #     geom_point(data = VOTE2018, mapping = aes(x = state, y = vote, fill=party, color=party),
  #                size=4,shape='x', position=position_dodge(width=0.5)) + 
  #     ggtitle("Posterior credible intervals of vote share for major party candidates") +
  #     theme(plot.title = element_text(hjust=0.5),
  #           panel.background = element_rect(fill = 'white', colour = 'white'))
  # }
  
  # write results to csv
  result <- data.frame(CYCLE,
                       STATE,
                       CANDIDATE,
                       POSTERIORMEAN,
                       POSTERIORSTD,
                       LOWER95,
                       UPPER95,
                       MEDIAN,
                       WIN,
                       VOTE,
                       PMEAN,
                       PARTY)
  
  names(result) <- tolower(names(result))
  
  write.csv(result,output_file)
  
  # output_file = paste('results/stan_NLZ', TYPE, '_' , test_year, 'day', horizons[a], '_', best_cv_idx[a] ,'.csv',sep='')
  output_file = paste('results/stan_NLZGP', TYPE, '_' , test_year, 'day', horizon, '_', b ,'.csv',sep='')
  
  result <- data.frame(NLZ)

  names(result) <- tolower(names(result))

  write.csv(result,output_file)
  
  # print("Test")
  # 
  # print(paste("Correct predictions: ",correct_predictions))
  # 
  # print(paste("Accuracy: ",correct_predictions/length(test_metadata)))
  # 
  # print(paste("Correlation: ",cor(PMEAN, VOTE)))
  # 
  # print(paste("RSME: ",mean(RMSE)))
  # 
  # print(paste("Ratio in 95% : ",1-Nout_test/length(PMEAN)))
  # 
  # print(paste("Predictive averaged nlZ: ",mean(NLZ)))

  # print(paste("Mean of predictive std: ",mean(PSTD)))
  # 
  # print(paste("Median of predictive std: ",median(PSTD)))
  # 
  # print(paste("Std of predictive std: ",sd(PSTD)))
}

# save the r session
# save.image(file = paste('models/LOOCV_',TYPE ,'.RData',sep=''))

# saveRDS(fit_objs, file = paste("models/",TYPE, "_", test_year, "_fit_objs.rds",sep=''))

# fit_objs = readRDS("models/GP_2018_fit_objs.rds")
# fit = fitobjs[[1]]
# mcmc_trace(fit, par=c("alpha","beta","ppb","eb"))

# fit_objs = readRDS("models/GP_2018_fit_objs.rds")
# HORIZON = c()
# MEAN = c()
# SE_MEAN = c()
# SD = c()
# L = c()
# U = c()
# N_EFF = c()
# RHAT = c()
# PARMS = c("alpha","beta","ppb","eb","year_sig")
# for (a in 1:length(horizons)) {
#   fit = fit_objs[[a]]
#   tmp = summary(fit, par=PARMS)
#   tmp = tmp$summary
#   for (p in PARMS) {
#     HORIZON = c(HORIZON, a)
#     MEAN = c(MEAN, tmp[p, 'mean'])
#     SE_MEAN = c(SE_MEAN, tmp[p, 'se_mean'])
#     SD = c(SD, tmp[p, 'sd'])
#     L = c(L, tmp[p, '2.5%'])
#     U = c(U, tmp[p, '97.5%'])
#     N_EFF = c(N_EFF, tmp[p, 'n_eff'])
#     RHAT = c(RHAT,tmp[p, 'Rhat'])
#   }
# }

# MEAN = round(MEAN, 2)
# SE_MEAN = round(SE_MEAN, 3)
# SD = round(SD,2)
# L = round(L,2)
# U = round(U,2)
# N_EFF = round(N_EFF,1)
# RHAT = round(RHAT, 3)
# 
# results <- data.frame(HORIZON,
#                       MEAN,
#                       SE_MEAN,
#                       SD,
#                       L,
#                       U,
#                       N_EFF,
#                       RHAT)
# 
# write.csv(results,"results/stan.csv")
