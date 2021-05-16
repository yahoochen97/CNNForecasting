library(ggridges)
library(ggplot2)
library(grid)
library(dplyr)

setwd('/Users/yahoo/Documents/WashU/CSE515T/Code/Gaussian Process')

horizons = c('0',
             '7',
             '14',
             '21',
             '28',
             '42',
             '56')

TYPE='GP'

best_cv_idx = read.csv(paste("results/", TYPE, "_opthyp.csv", sep=''));
best_cv_idx = best_cv_idx$opt_idx

test_year = 2016
STATE = "Arizona"
CANDIDATE = "McCain"


COLNAMES = c('horizon','Posterior_Vote')
HARRIS = data.frame(matrix(ncol = length(COLNAMES), nrow = 0))
colnames(HARRIS) = COLNAMES

for (a in length(horizons):1) {
    fit = readRDS(file = paste("models/",TYPE, "_", test_year,"day_", horizons[a] , "_fit.rds",sep=''))
    fit_params <- as.data.frame(fit)
    
    # load the prior files
    input_file = paste('results/LOO', TYPE, '_' , test_year, 'day', horizons[a], '_', best_cv_idx[a] ,'.csv',sep='')
    output_file = paste('results/stan_LOO', TYPE, '_' , test_year, 'day', horizons[a], '_', best_cv_idx[a] ,'.csv',sep='')
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
    
    i = 1
    for(tmp in 1:length(test_idx2)){
      if (test_metadata[[test_idx2[tmp]]][2]==STATE){
         i = tmp
      }
    }
    cycle = test_metadata[[test_idx2[i]]][1]
    vote = test_y[[test_idx2[i]]]
    candidates = data_test[data_test$state==STATE & data_test$cycle==test_year,c("candidate")]
    j = 1
    for(tmp in 1:length(candidates)){
      if (candidates[tmp]==CANDIDATE){
        j = tmp
      }
    }
    
    tmp = paste('test_y2[',i,',',j,']',sep='')
    pred = fit_params[[tmp]]
    
    harris = data.frame(matrix(ncol = length(COLNAMES), nrow = length(pred)))
    colnames(harris) = COLNAMES
    
    harris$Posterior_Vote = pred*100
    harris$horizon = as.numeric(horizons[a])
    HARRIS = rbind(HARRIS, harris)
}

v = 100*vote[j]
v = round(v, 2)
LABELS = rev(c("Election Day", "One week left", "Two weeks left",
           "Three weeks left", "For weeks left", "Six weeks left", "Eight weeks left"))

ggplot(HARRIS, aes(x = Posterior_Vote, y = reorder(horizon, desc(horizon)))) +
  geom_density_ridges(alpha=0.6) +
  scale_y_discrete(expand = c(0, 0), name = "Horizon",labels=LABELS) +
  scale_x_continuous(expand = c(0, 0), breaks = c(40,50,60,v, 70,80),
                     name = "Posterior Vote (%)") +
  geom_vline(xintercept=v, colour="blue") +
  theme(panel.grid.major.x = element_line(color = "gray"),
        panel.grid.major.y = element_line(size=.2, color="grey" )) +
  coord_cartesian(xlim = c(40, 80), clip='on') +
  theme(plot.title = element_text(hjust=0.5),
        panel.background = element_rect(fill = 'white', colour = 'white'))
