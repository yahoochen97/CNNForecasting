setwd('/Users/yahoo/Documents/WashU/CSE515T/Code/Gaussian Process/')

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is argument: if not, use default 2018
if (length(args)==0) {
  horizon = '0'
  test_year = 2018
}

horizon = args[1]
test_year = as.double(args[2])

# horizons = c('0',
#              '7',
#              '14',
#              '21',
#              '28',
#              '42',
#              '56')

library(rstan)

output_path = 'results/brw_'

computeh = function(df, test_year){
  if(test_year==2018){
    data = df[df$cycle!=test_year, ]
  }
  else{
    data = df[df$cycle!=test_year & df$cycle!=2018, ]
  }
  model <- lm(vote ~ pvi + experienced + Democrat + Republican + pvi*Republican, data) 
  return(model)
}

# read data
df <- read.csv("data/CNNdata1992-2016.csv")
colnames(df)[colnames(df) == "Percentage_of_Vote_won_x"] = 'vote'
df$vote = df$vote/100

polls = unique(df[,c("cycle", "state", "samplesize", "daysLeft","pollster")])
for (i in 1:nrow(polls)) {
  poll = polls[i,]
  cycle = poll$cycle[1]
  state = poll$state[1]
  samplesize = poll$samplesize[1]
  daysLeft = poll$daysLeft[1]
  pollster = poll$pollster[1]
  mask = df$cycle==cycle & df$state==state & df$samplesize==samplesize & df$daysLeft==daysLeft & df$pollster==pollster
  if(daysLeft==1){
    mask = df$cycle==cycle & df$state==state & df$samplesize==samplesize & df$daysLeft==daysLeft
  }
  if(sum(mask)==1 || sum(mask)!=length(unique(df[df$state==state&df$cycle==cycle,"Candidateidentifier"]))){
    df = df[!mask,]
    next
  }
  if(length(unique(df[mask,c("Candidateidentifier")]))!=length(df[mask,c("samplesize")])){
    print(df[mask,])
  }
  df[mask,c("samplesize")] = sum(df[mask,c("numberSupport")])
  df[mask, c("vote")] = df[mask, c("vote")]/sum(df[mask, c("vote")])
}

df2018 <- read.csv("data/CNNdata2018.csv")
colnames(df2018)[colnames(df2018) == "Percentage.of.Vote.won.x"] = 'vote'
df2018$vote = df2018$vote/100
df2018 = subset(df2018, select=-c(candidate_name))

polls = unique(df2018[,c("cycle", "state", "samplesize", "daysLeft","pollster")])
for (i in 1:nrow(polls)) {
  poll = polls[i,]
  cycle = poll$cycle[1]
  state = poll$state[1]
  samplesize = poll$samplesize[1]
  daysLeft = poll$daysLeft[1]
  pollster = poll$pollster[1]
  mask = df2018$cycle==cycle & df2018$state==state & df2018$samplesize==samplesize & df2018$daysLeft==daysLeft & df2018$pollster==pollster
  if(daysLeft==1){
    mask = df2018$cycle==cycle & df2018$state==state & df2018$samplesize==samplesize & df2018$daysLeft==daysLeft
  }
  if(sum(mask)==1 || sum(mask)!=length(unique(df2018[df2018$state==state&df2018$cycle==cycle,"Candidateidentifier"]))){
    df2018 = df2018[!mask,]
    next
  }
  df2018[mask, c("vote")] = df2018[mask, c("vote")]/sum(df2018[mask, c("vote")])
  df2018[mask,c("samplesize")] = sum(df2018[mask,c("numberSupport")])
}

df = rbind(df, df2018)

df$polls <- df$numberSupport/df$samplesize

priorModel = computeh(df, test_year)

df <- df[df$cycle==test_year, ]


horizons = c(horizon)

start.time <- Sys.time()

for (horizon in horizons) {
  output_file <- paste(output_path, test_year, "day_", horizon,'.csv',sep='')
  # slicing parameter
  W = 4
  precision = 1/0.1^2
  
  # iterate over candidate
  states = unique(df$state)
  I = 0
  NS <- c()
  YS <- c()
  CANDIDATES <- c()
  DAYS <- c()
  HS <- c()
  TAUS <- c()
  VS<- c()
  K = 0
  for(state in states){
    cs = unique(df[df$state==state, c('Candidateidentifier')])
    
    for (c in cs[1:(length(cs)-1)]) {
      # only iterate first C-1 candidates 
      data = df[as.character(df$Candidateidentifier)==c,]
      vote = data$vote[1]
      republican = data$Republican[1]
      democratic = data$Democrat[1]
      
      pvi = data$pvi[1]
      experienced = data$experienced[1]
      h = predict(priorModel, data[1,])[[1]]
  
      data = data[data$daysLeft<=-as.numeric(horizon),]
      n_poll = nrow(data)
      
      if(n_poll>0 && max(data$daysLeft)<0){
        # when there is available polls
        I = I + 1
        ns = as.array(data$samplesize)
        ys = as.array(data$numberSupport)
        days = as.array(ceiling(-data$daysLeft/W))
        k = length(ns)
        NS <- c(NS, ns)
        YS <- c(YS, ys)
        CANDIDATES <- c(CANDIDATES, rep(I, k))
        HS <- c(HS, h)
        TAUS <- c(TAUS, precision)
        VS<- c(VS, vote)
        DAYS <- c(DAYS, days)
        K = K + k
      }
    }
  }
  
  J = max(DAYS)
  DAYS = J + 1 - DAYS
  
  stan_data <- list(I=I,
                    J=J,
                    K=K,
                    n=NS,
                    y=YS,
                    candidate=CANDIDATES,
                    day=DAYS,
                    h=HS,
                    tau=TAUS,
                    v=VS)
  
  # define stan model
  # model <- stan_model("dynamicBayesian.stan")
  
  # train stan model
  fit <- stan(file = "dynamicBayesian.stan",
              data = stan_data, 
              warmup = 500, 
              iter = 5000, 
              chains = 1, 
              cores = 1, 
              thin = 4,
              control=list(adapt_delta=.98, max_treedepth = 15),
              seed = 1,
              refresh=0
  )
  
  fit_params <- as.data.frame(fit)
  
  saveRDS(fit, file = paste("models/brw_", test_year, "day_", horizon , "_fit.rds",sep=''))
  
  CYCLE <- c()
  STATE <- c()
  CANDIDATE <- c()
  
  PRIOR <- c()
  VOTE <- c()
  RMSE <- c()
  
  WIN <- c()
  WINNERS <- c()
  MEDIAN <- c()
  NLZ <- c()
  
  PMEAN <- c()
  LOWER95 <- c()
  UPPER95 <- c()
  
  i = 0
  
  for(state in states){
    cs = unique(df[df$state==state, c('Candidateidentifier')])
    preds = 0
    PREDS = c()
    votes = c()
    
    
    for (c_idx in 1:length(cs)) {
      c = as.character(cs[c_idx])
      data = df[as.character(df$Candidateidentifier)==c,]
      vote = data$vote[1]
      republican = data$Republican[1]
      democratic = data$Democrat[1]
      
      pvi = data$pvi[1]
      experienced = data$experienced[1]
      h = predict(priorModel, data[1,])[[1]]
      
      CYCLE <- c(CYCLE, test_year)
      STATE <- c(STATE, state)
      CANDIDATE <- c(CANDIDATE, substr(c,7,100))
      PRIOR <- c(PRIOR, h)
      VOTE <- c(VOTE, vote)
      votes = c(votes,vote)
      
      COMPUTE_LL = TRUE
      if(c_idx==length(cs)){
        # for the final candidate
        # apply the sum to 1 contraints
        pred = 1 - preds
        PREDS = c(PREDS, pred)
      }
      else{
        data = data[data$daysLeft<=-as.numeric(horizon),]
        n_poll = nrow(data)
  
        if(n_poll>0 && max(data$daysLeft)<0){
          # when there is available poll
          # use stan posteriors
          i = i + 1
          pred = fit_params[[paste("beta[",as.character(i) ,",", J, "]", sep="")]]
          pred = 1/(1+exp(-pred))
        }
        else{
          # otherwise use prior
          pred = rnorm(10000, mean = h, sd = 1/sqrt(precision))
          COMPUTE_LL = FALSE
        }
        preds = preds + pred
        PREDS = c(PREDS, pred)
      }
      
      m = as.numeric(quantile(pred, 0.5, na.rm = TRUE))
      u = as.numeric(quantile(pred, 0.975, na.rm = TRUE))
      l = as.numeric(quantile(pred, 0.025, na.rm = TRUE))
      PMEAN = c(PMEAN, mean(pred))
      LOWER95 = c(LOWER95, l)
      UPPER95 = c(UPPER95, u)
      MEDIAN <- c(MEDIAN, m)
      RMSE <- c(RMSE, sqrt(mean((pred - rep(vote,length(pred)))^2)))
      
      if(COMPUTE_LL){
        ll = c()
        for(t in 1:length(pred)){
          # logit normal pdf
          ll = c(ll, log(dnorm(log(vote/(1-vote)), log(pred[t]/(1-pred[t])) ,fit_params$sigma_beta[t])/vote/(1-vote)))
        }
        NLZ <- c(NLZ, -log(mean(exp(ll))))
      }
      else{
        NLZ <- c(NLZ,  -dnorm(vote,h, 1/sqrt(precision), log=TRUE))
      }
    }
    
    PREDS  <- matrix(PREDS , nrow = length(cs), byrow = TRUE)
    win_rates = rep(0, length(cs))
    for(k in 1:ncol(PREDS)){
      idx = which.max(PREDS [,k])
      win_rates[idx] = win_rates[idx] + 1
    }
    win_rates = win_rates / sum(win_rates)
    WIN <- c(WIN, win_rates)
    winners = rep(0, length(cs))
    winners[which.max(votes)] = 1
    WINNERS  <- c(WINNERS, winners)
  }
  
  
  for (i in 1:I){
    ms = c()
    upper = c()
    lower = c()
    for (j in 1:J){
      pred = fit_params[[paste("beta[",as.character(i) ,",", as.character(j) , "]", sep="")]]
      m = as.numeric(quantile(pred, 0.5, na.rm = TRUE))
      ms = c(ms, m)
      u = as.numeric(quantile(pred, 0.975, na.rm = TRUE))
      upper = c(upper, u)
      l = as.numeric(quantile(pred, 0.025, na.rm = TRUE))
      lower = c(lower, l)
    }
    PMEAN = c(PMEAN, 1/(1+exp(-m)))
    LOWER95 = c(LOWER95, 1/(1+exp(-l)))
    UPPER95 = c(UPPER95, 1/(1+exp(-u)))
    jpeg(paste("plots/brw/",as.character(i), ".jpg", sep=""), width = 1000, height = 1000)
    plot(-(J:1),1/(1+exp(-ms)), type="l",ylim=c(0,1))
    lines(-(J:1),1/(1+exp(-upper)),type="l")
    lines(-(J:1),1/(1+exp(-lower)),type="l")
    points(-(J+1-DAYS[CANDIDATES==i]), (YS[CANDIDATES==i]/NS[CANDIDATES==i]))
    points(0, VS[i], col="blue",pch=20, cex=2)
    dev.off()
  }
  
  
  # write results to csv
  result <- data.frame(CYCLE,
                       STATE,
                       CANDIDATE,
                       PRIOR,
                       PMEAN,
                       VOTE,
                       LOWER95,
                       UPPER95,
                       MEDIAN,
                       WIN,
                       WINNERS,
                       RMSE,
                       NLZ)
  
  names(result) <- tolower(names(result))
  
  write.csv(result,output_file)
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken