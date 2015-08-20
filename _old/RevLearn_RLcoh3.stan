data {
  int<lower=1> nSubjects;                              // number of subjects
  int<lower=1> nTrials;                                // number of trials 
  int<lower=1,upper=2> choice1[nSubjects,nTrials];     // 1st choices, 1 or 2
  int<lower=1,upper=2> choice2[nSubjects,nTrials];     // 2nd choices, 1 or 2
  real<lower=-1,upper=1> reward[nSubjects,nTrials];    // outcome, 1 or -1
  real<lower=0,upper=1>  with[nSubjects,nTrials];      // No. of with / 4
  real<lower=0,upper=1>  against[nSubjects,nTrials];   // No. of against / 4
}

// transformed data{
//   real L;
//   real U;
// 
//   L <- -10.0;
//   U <-  10.0;
// }

parameters {
  // group-level parameters
  real lr_mu_pr;    // lr_mu before probit
  real tau_mu_pr;   // tau_mu before probit  
  //real<lower=L,upper=U> coha_mu;
  //real<lower=L,upper=U> cohw_mu;
  real coha_mu;
  real cohw_mu;

  real<lower=0> lr_sd;
  real<lower=0> tau_sd;
  real<lower=0.01, upper=2> coha_sd;
  real<lower=0.01, upper=2> cohw_sd;
  
  // subject-level row parameters, follows norm(0,1), for later Matt Trick
  vector[nSubjects] lr_row;
  vector[nSubjects] tau_row;
  
  //vector<lower=L,upper=U>[nSubjects] coha;
  //vector<lower=L,upper=U>[nSubjects] cohw;
  vector[nSubjects] coha;
  vector[nSubjects] cohw;
}

transformed parameters {
  // subject-level parameters
  vector<lower=0,upper=1>[nSubjects] lr;
  vector<lower=0,upper=10>[nSubjects] tau;

  // Matt Trick, note that the input of Phi_approx must be 'real' rather than 'vector'
  for (s in 1:nSubjects) {
    lr[s] <- Phi_approx( lr_mu_pr + lr_sd * lr_row[s] );
    tau[s] <- 10 * Phi_approx( tau_mu_pr + tau_sd * tau_row[s] );
  }
}

model {
  // define the value and pe vectors
  vector[2] v[nTrials+1];    //  values
  vector[nTrials] pe;         // prediction errors
  vector[nTrials] penc;       // pe for the non-chosen choice

  // hyperparameters
  lr_mu_pr   ~ normal(0,1);
  tau_mu_pr ~ normal(0,1);
  coha_mu    ~ normal(0,1);// T[L,U];
  cohw_mu    ~ normal(0,1);// T[L,U];

  lr_sd   ~ cauchy(0,5);
  tau_sd ~ cauchy(0,5);
  
  // Matt Trick
  lr_row   ~ normal(0,1);
  tau_row ~ normal(0,1);

  coha ~ normal(coha_mu,coha_sd);
  cohw ~ normal(cohw_mu,cohw_sd);

  // subject loop and trial loop
  for (s in 1:nSubjects) {
    //coha[s] ~ normal(coha_mu,coha_sd) T[L,U];
    //cohw[s] ~ normal(cohw_mu,cohw_sd) T[L,U];

    v[1] <- rep_vector(0.0,2);
    
    for (t in 1:nTrials) {
      
      v[t][3-choice1[s,t]] <- v[t][3-choice1[s,t]] + coha[s] * against[s,t]*4;
      v[t][choice1[s,t]]   <- v[t][choice1[s,t]]   + cohw[s] * with[s,t]*4;
      
      choice2[s,t] ~ categorical_logit( tau[s] * v[t] );
      pe[t]   <-  reward[s,t] - v[t][choice2[s,t]];
      penc[t] <- -reward[s,t] - v[t][3-choice2[s,t]];

      v[t+1][choice2[s,t]]   <- v[t][choice2[s,t]]   + lr[s] * pe[t]; 
      v[t+1][3-choice2[s,t]] <- v[t][3-choice2[s,t]] + lr[s] * penc[t];
    }
  }
}

generated quantities {
  real<lower=0,upper=1> lr_mu; 
  real<lower=0,upper=10> tau_mu;
  
  real log_lik[nSubjects]; 
  vector[2] v2[nTrials+1];
  vector[nTrials] pe2;
  vector[nTrials] penc2;
  
  lr_mu <- Phi_approx(lr_mu_pr);
  tau_mu <- 10 * Phi_approx(tau_mu_pr);
  
  for (s in 1:nSubjects) {
    log_lik[s] <- 0;
    v2[1] <- rep_vector(0.0,2);
    
    for (t in 1:nTrials) {
      v2[t][3-choice1[s,t]] <- v2[t][3-choice1[s,t]] + coha[s] * against[s,t]*4;
      v2[t][choice1[s,t]]   <- v2[t][choice1[s,t]]   + cohw[s] * with[s,t]*4;
      
      log_lik[s] <- log_lik[s] + categorical_logit_log(choice2[s,t], tau[s] * v2[t]);
      
      pe2[t]   <-  reward[s,t] - v2[t][choice2[s,t]];
      penc2[t] <- -reward[s,t] - v2[t][3-choice2[s,t]];

      v2[t+1][choice2[s,t]]   <- v2[t][choice2[s,t]]   + lr[s] * pe2[t];
      v2[t+1][3-choice2[s,t]] <- v2[t][3-choice2[s,t]] + lr[s] * penc2[t];
    }
  }
}