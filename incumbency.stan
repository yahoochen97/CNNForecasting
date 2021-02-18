
data {
  int<lower=0> N2; // number of races
  real mu2[N2, 2]; // vector of prior Gaussian means
  real sigma2[N2, 2]; // vector of prior Gaussian stds
  // int<lower=0> nc[N]; // number of candidates each race
  vector[2] pvi2[N2];
  vector[2] party2[N2];
  vector[2] experienced2[N2];
  vector[2] incumbency2[N2];
  vector<lower=0, upper=1>[2] y2[N2]; // vector of true votings
  int year_idx2[N2];
  
  int<lower=0> N3; // number of races
  real mu3[N3, 3]; // vector of prior Gaussian means
  real sigma3[N3, 3]; // vector of prior Gaussian stds
  vector[3] pvi3[N3];
  vector[3] party3[N3];
  vector[3] experienced3[N3];
  vector[3] incumbency3[N3];
  vector<lower=0, upper=1>[3] y3[N3]; // vector of true votings
  int year_idx3[N3];
  
  int<lower=0> N4; // number of races
  real mu4[N4, 4]; // vector of prior Gaussian means
  real sigma4[N4, 4]; // vector of prior Gaussian stds
  vector[4] pvi4[N4];
  vector[4] party4[N4];
  vector[4] experienced4[N4];
  vector[4] incumbency4[N4];
  vector<lower=0, upper=1>[4] y4[N4]; // vector of true votings
  int year_idx4[N4];
  
  int<lower=0> test_N2; // number of races
  real test_mu2[test_N2, 2]; // vector of prior Gaussian means
  real test_sigma2[test_N2, 2]; // vector of prior Gaussian stds
  vector[2] test_pvi2[test_N2];
  vector[2] test_party2[test_N2];
  vector[2] test_experienced2[test_N2];
  vector[2] test_incumbency2[test_N2];
  int test_year_idx2[test_N2];
  vector<lower=0, upper=1>[2] test_f2[test_N2]; // vector of true votings
  
  int<lower=0> test_N3; // number of races
  real test_mu3[test_N3, 3]; // vector of prior Gaussian means
  real test_sigma3[test_N3, 3]; // vector of prior Gaussian stds
  vector[3] test_pvi3[test_N3];
  vector[3] test_party3[test_N3];
  vector[3] test_experienced3[test_N3];
  vector[3] test_incumbency3[test_N3];
  int test_year_idx3[test_N3];
  vector<lower=0, upper=1>[3] test_f3[test_N3]; 
  
  int<lower=0> test_N4; // number of races
  real test_mu4[test_N4, 4]; // vector of prior Gaussian means
  real test_sigma4[test_N4, 4]; // vector of prior Gaussian stds
  vector[4] test_pvi4[test_N4];
  vector[4] test_party4[test_N4];
  vector[4] test_experienced4[test_N4];
  vector[4] test_incumbency4[test_N4];
  // int<lower=0> test_nc[test_N]; // number of candidates each race
  int test_year_idx4[test_N4];
  vector<lower=0, upper=1>[4] test_f4[test_N4]; 
  int max_year_idx;
  
  
}

parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0> ppb; // party*pvi b
  real<lower=0> eb; // experienced b
  real<lower=0> ib; // incumbency b
  real<lower=0> year_sig; // The SD for the year-level random effects
  vector[max_year_idx] year_re; // Year level random effects in the vote share equation
  vector<lower=0, upper=1>[2] gamma2[N2]; // vector of underlying true support rates
  vector<lower=0, upper=1>[3] gamma3[N3];
  vector<lower=0, upper=1>[4] gamma4[N4];
  
  
}

transformed parameters{
  vector<lower=0>[2] p2[N2]; // vector of underlying dirichlet parameters
  vector<lower=0>[3] p3[N3];
  vector<lower=0>[4] p4[N4];
  
  for(i in 1:N2){
    for(j in 1:2){
        // a linear model for the concentration parameters in the Dirichlet distribution
        p2[i,j] = beta*gamma2[i,j] + 
            party2[i,j]*pvi2[i,j]*ppb + experienced2[i,j]*eb + incumbency2[i,j]*ib
            + party2[i,j]*year_re[year_idx2[i]];
        // numerical stablity
        if(p2[i,j]<0.0001){
          p2[i,j] = 0.0001;
        }
    }
  }
  
  for(i in 1:N3){
    for(j in 1:3){
        p3[i,j] = beta*gamma3[i,j] + 
            party3[i,j]*pvi3[i,j]*ppb + experienced3[i,j]*eb + incumbency3[i,j]*ib
            + party3[i,j]*year_re[year_idx3[i]];
        if(p3[i,j]<0.0001){
          p3[i,j] = 0.0001;
        }
    }
  }
  
  for(i in 1:N4){
    for(j in 1:4){
        p4[i,j] = beta*gamma4[i,j] + 
            party4[i,j]*pvi4[i,j]*ppb + experienced4[i,j]*eb + incumbency4[i,j]*ib
            + party4[i,j]*year_re[year_idx4[i]];
        if(p4[i,j]<0.0001){
          p4[i,j] = 0.0001;
        }
    }
  }
  
  
}

model {
  // very flat priors on linear model parameters
  alpha ~ normal(0, 100);
  beta ~ normal(0, 100);
  ppb ~ normal(0, 10);
  eb ~ normal(0, 10);
  ib ~ normal(0,10);
  year_sig ~ gamma(1, .5); // Prior on the SD for the year-level RE
  year_re ~ normal(0, year_sig); // Prior for the year-level RE
  // generate underlying true support rates
  for(i in 1:N2){
    gamma2[i] ~ normal(mu2[i], sigma2[i]);
  }
  
  for(i in 1:N3){
    gamma3[i] ~ normal(mu3[i], sigma3[i]);
  }
  
  for(i in 1:N4){
    gamma4[i] ~ normal(mu4[i], sigma4[i]);
  }

  // generate voting rates from support rates
  for(i in 1:N2){
    y2[i] ~ dirichlet(alpha+p2[i]);
  }
  
  for(i in 1:N3){
    y3[i] ~ dirichlet(alpha+p3[i]);
  }
  
  for(i in 1:N4){
    y4[i] ~ dirichlet(alpha+p4[i]);
  }
  
  
}

generated quantities {
  vector[2] test_y2[test_N2]; // vector of true votings
  vector<lower=0>[2] test_p2[test_N2]; // vector of underlying dirichlet parameters
  vector<lower=0, upper=1>[2] test_gamma2[test_N2]; // vector of underlying true support rates
  
  vector[3] test_y3[test_N3]; // vector of true votings
  vector<lower=0>[3] test_p3[test_N3]; // vector of underlying dirichlet parameters
  vector<lower=0, upper=1>[3] test_gamma3[test_N3]; // vector of underlying true support rates
  
  vector[4] test_y4[test_N4]; // vector of true votings
  vector<lower=0>[4] test_p4[test_N4]; // vector of underlying dirichlet parameters
  vector<lower=0, upper=1>[4] test_gamma4[test_N4]; // vector of underlying true support rates
  
  vector[2] rep_y2[N2]; // vector of forecasting votings
  
  vector[3] rep_y3[N3]; // vector of forecasting votings
  
  vector[4] rep_y4[N4]; // vector of forecasting votings
  
  vector[N2] ll2; // vector of in-sample model log likelihood

  vector[N3] ll3;

  vector[N4] ll4;
  
  vector[test_N2] test_ll2;  // vector of out-of-sample model log likelihood

  vector[test_N3] test_ll3;

  vector[test_N4] test_ll4;
  
  // generate in-sample forecasting
  for (i in 1:N2)
    rep_y2[i] = dirichlet_rng(alpha+p2[i]);
    
  for (i in 1:N3)
    rep_y3[i] = dirichlet_rng(alpha+p3[i]);
    
  for (i in 1:N4)
    rep_y4[i] = dirichlet_rng(alpha+p4[i]);
    
  // obtain in-sample log likelihood
  for (i in 1:N2)
    ll2[i] = dirichlet_lpdf(y2[i]|alpha+p2[i]);

  for (i in 1:N3)
    ll3[i] = dirichlet_lpdf(y3[i]|alpha+p3[i]);

  for (i in 1:N4)
    ll4[i] = dirichlet_lpdf(y4[i]|alpha+p4[i]);
  
  // generate samples for voter preference priors as 
  // one of the independent variables for the concentration
  // parameters for test data
  for(i in 1:test_N2){
    for(j in 1:2){
      test_gamma2[i,j] = normal_rng(test_mu2[i,j], test_sigma2[i,j]);
      if(test_gamma2[i,j]<0){
        test_gamma2[i,j] = 0;
      }
      else if(test_gamma2[i,j]>1){
        test_gamma2[i,j]=1;
      }
    }
  }
  
  for(i in 1:test_N3){
    for(j in 1:3){
      test_gamma3[i,j] = normal_rng(test_mu3[i,j], test_sigma3[i,j]);
      if(test_gamma3[i,j]<0){
        test_gamma3[i,j] = 0;
      }
      else if(test_gamma3[i,j]>1){
        test_gamma3[i,j]=1;
      }
    }
  }
  
  for(i in 1:test_N4){
    for(j in 1:4){
      test_gamma4[i,j] = normal_rng(test_mu4[i,j], test_sigma4[i,j]);
      if(test_gamma4[i,j]<0){
        test_gamma4[i,j] = 0;
      }
      else if(test_gamma4[i,j]>1){
        test_gamma4[i,j]=1;
      }
    }
  }
  
  // the linear model for the concentration parameters in the Dirichlet distribution
  for(i in 1:test_N2){
    for(j in 1:2){
      test_p2[i,j] = beta*test_gamma2[i,j] + 
            test_party2[i,j]*test_pvi2[i,j]*ppb + test_experienced2[i,j]*eb + test_incumbency2[i,j]*ib + 
            test_party2[i,j]*year_re[test_year_idx2[i]];
      if(test_p2[i,j]<0.0001){
          test_p2[i,j] = 0.0001;
        }
    }
  }
  
  for(i in 1:test_N3){
    for(j in 1:3){
      test_p3[i,j] = beta*test_gamma3[i,j] + 
            test_party3[i,j]*test_pvi3[i,j]*ppb + test_experienced3[i,j]*eb +test_incumbency3[i,j]*ib + 
            test_party3[i,j]*year_re[test_year_idx3[i]];
      if(test_p3[i,j]<0.0001){
          test_p3[i,j] = 0.0001;
        }
    }
  }
  
  for(i in 1:test_N4){
    for(j in 1:4){
      test_p4[i,j] = beta*test_gamma4[i,j] + 
            test_party4[i,j]*test_pvi4[i,j]*ppb + test_experienced4[i,j]*eb +test_incumbency4[i,j]*ib + 
            test_party4[i,j]*year_re[test_year_idx4[i]];
      if(test_p4[i,j]<0.0001){
          test_p4[i,j] = 0.0001;
        }
    }
  }
  
  // generate out-of-sample forecasting
  for (i in 1:test_N2)
    test_y2[i] = dirichlet_rng(alpha+test_p2[i]);
    
  for (i in 1:test_N3)
    test_y3[i] = dirichlet_rng(alpha+test_p3[i]);
    
  for (i in 1:test_N4)
    test_y4[i] = dirichlet_rng(alpha+test_p4[i]);
    
  // obtain out-of-sample log likelihood
  for (i in 1:test_N2)
    test_ll2[i] = dirichlet_lpdf(test_f2[i]|alpha+test_p2[i]);

  for (i in 1:test_N3)
    test_ll3[i] = dirichlet_lpdf(test_f3[i]|alpha+test_p3[i]);

  for (i in 1:test_N4)
    test_ll4[i] = dirichlet_lpdf(test_f4[i]|alpha+test_p4[i]);


}
