library(rstan)
library(MCMCpack)

# loading data
data = read.csv("results/forecast1992-2016.csv")
data2016 = data[data$cycle==2016,]
data = data[data$cycle!=2016,]
cycles = unique(data$cycle)
states = unique(data$state)

# define variables
metadata = list()
mu = list()
sigma = list()
y = list()
nc = c()
counter = 0

# iterate over races
for (cycle in cycles) {
  for (state in states) {
    pmu = data[data$cycle==cycle & data$state==state,c("posteriormean")]
    pstd = data[data$cycle==cycle & data$state==state,c("posteriorstd")]
    vote = data[data$cycle==cycle & data$state==state,c("vote")]
    if(length(pmu)){
      counter = counter + 1
      metadata[[counter]] = c(cycle, state)
      mu[[counter]] = pmu
      sigma[[counter]] = pstd
      y[[counter]] = vote / 100
      nc = c(nc, length(vote))
    }
  }
}

# build stan data
stan_mu = matrix(0,counter,4)
stan_sigma = matrix(0.0001,counter,4)
stan_y = matrix(0.0001,counter,4)

for (i in 1:counter) {
  stan_mu[i,1:nc[i]] = mu[[i]]
  stan_sigma[i,1:nc[i]] = sigma[[i]]
  stan_y[i,1:nc[i]] = y[[i]]
  stan_y[i,] = stan_y[i,]/sum(stan_y[i,])
}

# define stan data structure
stan_data = list(N = counter, 
                 mu = stan_mu, 
                 sigma = stan_sigma,
                 nc = nc,
                 y = stan_y)

# define stan model
model <- stan_model("model.stan")

# train stan model
fit <- stan(file = "model.stan",
            data = stan_data, 
            warmup = 500, 
            iter = 3000, 
            chains = 3, 
            cores = 3, 
            thin = 4,
            control=list(adapt_delta=.95)
            )

# summary(fit)
fit_params = extract(fit)

# test data
metadata = list()
mu = list()
sigma = list()
y = list()
nc = c()
counter = 0

# iterate over races
for (state in states) {
  pmu = data2016[data2016$state==state,c("posteriormean")]
  pstd = data2016[data2016$state==state,c("posteriorstd")]
  vote = data2016[data2016$state==state,c("vote")]
  if(length(pmu)){
    counter = counter + 1
    metadata[[counter]] = c(state)
    mu[[counter]] = pmu
    sigma[[counter]] = pstd
    y[[counter]] = vote / 100
    nc = c(nc, length(vote))
  }
}

# build stan data
stan_mu = matrix(0,counter,4)
stan_sigma = matrix(0.0001,counter,4)
stan_y = matrix(0.0001,counter,4)

for (i in 1:counter) {
  stan_mu[i,1:nc[i]] = mu[[i]]
  stan_sigma[i,1:nc[i]] = sigma[[i]]
  stan_y[i,1:nc[i]] = y[[i]]
  stan_y[i,] = stan_y[i,]/sum(stan_y[i,])
}

# prediction

sample_posterior <- function(mus, sigmas, nc, gs=1, ds=1, fit_params=fit_params){
  n = length(fit_params$alpha)
  preds = matrix(0, nc, n*gs*ds)
  count = 1
  for(k in 1:n){
    alpha = fit_params$alpha[i]
    beta = fit_params$beta[i]
    # g: monte carlo normal
    for(g in 1:gs){
      p = rep(0, 4)
      for(j in 1:4){
        if(j<=nc){
          gamma = rnorm(1, mean = mus[j], sd = sigmas[j])
          gamma = min(max(gamma, 0),1)
          p[j] = alpha + beta*gamma
        }
        else{
          p[j]=0.0001;
        }
      }
      pred = rdirichlet(ds, p[1:nc])
      preds[,count:(count+ds-1)] = t(pred)
      count = count + ds
    }
  }
  return(preds)
}

# within 95% CI
flags = matrix(0, counter, 4)

CYCLE = c()
STATE = c()
CANDIDATE = c()
POSTERIORMEAN = c()
POSTERIORSTD = c()
VOTE = c()
LOWER95 = c()
UPPER95 = c()

for(i in 1:counter) {
  state = metadata[[i]]
  pmu = data2016[data2016$state==state,c("posteriormean")]
  pstd = data2016[data2016$state==state,c("posteriorstd")]
  vote = data2016[data2016$state==state,c("vote")]
  candidates = data2016[data2016$state==state,c("candidate")]
  preds <- sample_posterior(stan_mu[i,], stan_sigma[i,], nc[i], gs=100, ds=100, fit_params=fit_params)
  for(j in 1:nc[i]){
    u=quantile(preds[j,],probs=c(0.975),names = FALSE)
    l=quantile(preds[j,],probs=c(0.025),names = FALSE)
    if (stan_y[i,j]<=u & stan_y[i,j]>=l){
      flags[i,j] = 1
    }
    print(state)
    print(as.character(candidates[j]))
    CYCLE = c(CYCLE, 2016)
    STATE = c(STATE,state)
    CANDIDATE = c(CANDIDATE,as.character(candidates[j]))
    POSTERIORMEAN = c(POSTERIORMEAN,pmu[j])
    POSTERIORSTD = c(POSTERIORSTD,pstd[j])
    VOTE = c(VOTE, vote[j])
    LOWER95 = c(LOWER95, l)
    UPPER95 = c(UPPER95, u)
  }
}

result = data.frame(CYCLE,
                    STATE,
                    CANDIDATE,
                    POSTERIORMEAN,
                    POSTERIORSTD,
                    VOTE,
                    LOWER95,
                    UPPER95)

names(result) <- tolower(names(result))

write.csv(result,'results/stan_prediction.csv')