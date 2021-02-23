horizons = c('0',
             '7',
             '14',
             '21',
             '28',
             '42',
              '56')

TYPE='incumbent_GP'

best_cv_idx = read.csv(paste("results/GP_opthyp.csv", sep=''));
best_cv_idx = best_cv_idx$opt_idx

# fit_objs = readRDS(file = paste("models/",TYPE, "_", test_year, "_fit_objs.rds",sep=''))

test_years = c(1992,1994,1996,1998,2000,2002,2004,2006,2008,2010,2012,2014,2016,2018)

test_years = c(2018, 2020)

for (a in 1:length(horizons)) {
  for (test_year in test_years) {
    fit = readRDS(file = paste("RR/models/",TYPE, "_", test_year,"day_", horizons[a] , "_fit.rds",sep=''))
    fit_params <- as.data.frame(fit)
    
    # load the prior files
    input_file = paste('RR/results/', TYPE, '_' , test_year, 'day', horizons[a], '_', best_cv_idx[a] ,'.csv',sep='')
    output_file = paste('RR/results/stan_LOO', TYPE, '_' , test_year, 'day', horizons[a], '_', best_cv_idx[a] ,'.csv',sep='')
    data <- read.csv(input_file)
    print(input_file)
    
    # remove unlike candidates of races with >4 candidates
    data <- data[data$cycle!=2016 | data$state!='Louisiana' | data$candidate!='Flemsing',]
    data <- data[data$cycle!=2020 | data$state!='Georgia' | data$candidate!='Loeffler',]
    data <- data[data$cycle!=2020 | data$state!='Georgia' | data$candidate!='Tarver',]
    
    # split training and testing data
    data_test <- data[(data$cycle==test_year),]
    data <- data[(data$cycle!=test_year & data$cycle!=2018 & data$cycle!=2020),]
    
    states <- union(unique(data$state), unique(data_test$state))
    
    # test data
    test_metadata <- list()
    test_y <- list()
    test_counter <- 0
    
    # iterate over races
    for (cycle in unique(data_test$cycle)){
      for (state in states) {
        vote = data_test[data_test$state==state & data_test$cycle==cycle,c("vote")]
        if(length(vote)){
          test_counter = test_counter + 1
          test_metadata[[test_counter]] = c(cycle,state)
          test_y[[test_counter]] = vote / sum(vote)
        }
      }
    }
    
    test_idx2 <- c()
    test_idx3 <- c()
    test_idx4 <- c()
    
    for (i in 1:test_counter) {
      tmp = length(test_y[[i]])
      if(tmp==2) test_idx2 = c(test_idx2, i)
      if(tmp==3) test_idx3 = c(test_idx3, i)
      if(tmp==4) test_idx4 = c(test_idx4, i)
    }
    
    CYCLE <- c()
    STATE <- c()
    CANDIDATE <- c()
    VOTE <- c()
    RMSE = c()
    LOWER95 <- c()
    UPPER95 <- c()
    MEDIAN <- c()
    WIN <- c()
    WINNERS <- c()
    NLZ <- c()
    PSTD <- c()
    
    correct_predictions <- 0
    Nout_test <- 0
    Nout <- 0
    
    for(i in 1:length(test_idx2)) {
      cycle = test_metadata[[test_idx2[i]]][1]
      state = test_metadata[[test_idx2[i]]][2]
      vote = test_y[[test_idx2[i]]]
      candidates = data_test[data_test$state==state & data_test$cycle==cycle,c("candidate")]
      preds= c()
      for(j in 1:2){
        tmp = paste('test_y2[',i,',',j,']',sep='')
        pred = fit_params[[tmp]]
        preds = c(preds, pred)
        u=quantile(pred,probs=c(0.975),names = FALSE)
        l=quantile(pred,probs=c(0.025),names = FALSE)
        
        if (vote[j]>u | vote[j]<l){
          Nout_test = Nout_test + 1
        }
        
        CYCLE <- c(CYCLE, cycle)
        STATE <- c(STATE,state)
        CANDIDATE <- c(CANDIDATE,as.character(candidates[j]))
        RMSE = c(RMSE, sqrt(mean((pred - rep(vote[j],length(pred)))^2)))
        VOTE <- c(VOTE, vote[j])
        MEDIAN <- c(MEDIAN, median(pred))
        LOWER95 <- c(LOWER95, l)
        UPPER95 <- c(UPPER95, u)
        PSTD <- c(PSTD, sd(pred))
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
        vote = test_y[[test_idx3[i]]]
        candidates = data_test[data_test$state==state & data_test$cycle==cycle,c("candidate")]
        preds= c()
        for(j in 1:3){
          tmp = paste('test_y3[',i,',',j,']',sep='')
          pred = fit_params[[tmp]]
          preds = c(preds, pred)
          u=quantile(pred,probs=c(0.975),names = FALSE)
          l=quantile(pred,probs=c(0.025),names = FALSE)
          if (vote[j]>u | vote[j]<l){
            Nout_test = Nout_test + 1
          }
          CYCLE <- c(CYCLE, cycle)
          STATE <- c(STATE,state)
          CANDIDATE <- c(CANDIDATE,as.character(candidates[j]))
          RMSE = c(RMSE, sqrt(mean((pred - rep(vote[j],length(pred)))^2)))
          VOTE <- c(VOTE, vote[j])
          MEDIAN <- c(MEDIAN, median(pred))
          LOWER95 <- c(LOWER95, l)
          UPPER95 <- c(UPPER95, u)
          PSTD <- c(PSTD, sd(pred))
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
    
    if (length(test_idx4)){
      for(i in 1:length(test_idx4)) {
        cycle = test_metadata[[test_idx4[i]]][1]
        state = test_metadata[[test_idx4[i]]][2]
        vote = test_y[[test_idx4[i]]]
        party = data_test[data_test$state==state & data_test$cycle==cycle,c("party")]
        candidates = data_test[data_test$state==state & data_test$cycle==cycle,c("candidate")]
        preds= c()
        for(j in 1:4){
          tmp = paste('test_y4[',i,',',j,']',sep='')
          pred = fit_params[[tmp]]
          preds = c(preds, pred)
          u=quantile(pred,probs=c(0.975),names = FALSE)
          l=quantile(pred,probs=c(0.025),names = FALSE)
          if (vote[j]>u | vote[j]<l){
            Nout_test = Nout_test + 1
          }
          CYCLE <- c(CYCLE, cycle)
          STATE <- c(STATE, state)
          CANDIDATE <- c(CANDIDATE,as.character(candidates[j]))
          VOTE <- c(VOTE, vote[j])
          RMSE = c(RMSE, sqrt(mean((pred - rep(vote[j],length(pred)))^2)))
          MEDIAN <- c(MEDIAN, median(pred))
          LOWER95 <- c(LOWER95, l)
          UPPER95 <- c(UPPER95, u)
          PSTD <- c(PSTD, sd(pred))
        }
        NLZ <- c(NLZ, -log(mean(exp(fit_params[[paste('test_ll4[',i,']',sep='')]]))))
        preds <- matrix(preds, nrow = 4, byrow = TRUE)
        win_rates = rep(0, 4)
        for(k in 1:ncol(preds)){
          idx = which.max(preds[,k])
          win_rates[idx] = win_rates[idx] + 1
        }
        win_rates = win_rates / sum(win_rates)
        WIN <- c(WIN, win_rates)
        winners = rep(0, 4)
        winners[which.max(vote)] = 1
        WINNERS  <- c(WINNERS, winners)
        if (which.max(win_rates)==which.max(vote)){
          correct_predictions = correct_predictions + 1
        }
        else{
          print("Wrong prediction:")
          print(test_metadata[[test_idx4[i]]])
        }
      }
    }
    
    # write results to csv
    result <- data.frame(CYCLE,
                         STATE,
                         CANDIDATE,
                         LOWER95,
                         UPPER95,
                         MEDIAN,
                         WIN,
                         WINNERS,
                         RMSE,
                         VOTE)
    print(output_file)
    
    names(result) <- tolower(names(result))
    
    write.csv(result,output_file)
    
    output_file = paste('RR/results/stan_NLZ', TYPE, '_' , test_year, 'day', horizons[a], '_', best_cv_idx[a] ,'.csv',sep='')
    
    result <- data.frame(NLZ)
    
    names(result) <- tolower(names(result))
    
    write.csv(result,output_file)
    
    # print(paste("Correct predictions: ",correct_predictions))
    # 
    # print(paste("Accuracy: ",correct_predictions/length(test_metadata)))
    # 
    # print(paste("RSME: ",mean(RMSE)))
    # 
    # print(paste("Ratio in 95% : ",1-Nout_test/length(RMSE)))
    # 
    # print(paste("Predictive averaged nlZ: ",mean(NLZ)))
  }
}
