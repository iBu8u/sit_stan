data {
  int<lower=1> nSubjects;                    // number of subjects
  int<lower=1> nTrials;                      // number of trials 
  int<lower=1,upper=2> choice1[nTrials];     // 1st choices, 1 or 2
  int<lower=1,upper=2> choice2[nTrials];     // 2nd choices, 1 or 2
  real<lower=-1,upper=1> reward[nTrials];    // outcome, 1 or -1
  real<lower=0,upper=1>  with[nTrials];      // No. of with / 4
  real<lower=0,upper=1>  against[nTrials];   // No. of against / 4
  int<lower=0> obs[nSubjects+1];   
  int<lower=0> cumobs[nSubjects+1];
}

transformed data{
  int L;
  int U;
  int maxTrials;
  
  L <- -10;
  U <-  10;
  maxTrials <- 100;
}

parameters {
  // group-level parameters
  real lr_mu_pr;    // lr_mu before probit
  real tau_mu_pr;   // tau_mu before probit  
  real<lower=L,upper=U> coha_mu;
  real<lower=L,upper=U> cohw_mu;

  real<lower=0> lr_sd;
  real<lower=0> tau_sd;
  real<lower=0.01, upper=2> coha_sd;
  real<lower=0.01, upper=2> cohw_sd;
  
  // subject-level row parameters, follows norm(0,1), for later Matt Trick
  vector[nSubjects] lr_row;
  vector[nSubjects] tau_row;
  
  vector<lower=L,upper=U>[nSubjects] coha;
  vector<lower=L,upper=U>[nSubjects] cohw;
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
  vector[2] v[maxTrials+1];    //  values
  vector[maxTrials] pe;         // prediction errors
  vector[maxTrials] penc;       // pe for the non-chosen choice

  // hyperparameters
  lr_mu_pr   ~ normal(0,1);
  tau_mu_pr ~ normal(0,1);
  coha_mu    ~ normal(0,1) T[L,U];
  cohw_mu    ~ normal(0,1) T[L,U];

  lr_sd   ~ cauchy(0,5);
  tau_sd ~ cauchy(0,5);
  
  // Matt Trick
  lr_row   ~ normal(0,1);
  tau_row ~ normal(0,1);

  // subject loop and trial loop
  for (s in 1:nSubjects) {
    coha[s] ~ normal(coha_mu,coha_sd) T[L,U];
    cohw[s] ~ normal(cohw_mu,cohw_sd) T[L,U];
    v[1] <- rep_vector(0.0,2);
    
    for (t in 1:obs[s+1] ) {
      
      v[t][3-choice1[t+cumobs[s]]] <- v[t][3-choice1[t+cumobs[s]]] + coha[s] * against[t+cumobs[s]] * 4;
      v[t][choice1[t+cumobs[s]]]   <- v[t][choice1[t+cumobs[s]]]   + cohw[s] * with[t+cumobs[s]] * 4;
      
      choice2[t + cumobs[s]] ~ categorical_logit( tau[s] * v[t] );
      
      pe[t]   <-  reward[t+cumobs[s]] - v[t][choice2[t+cumobs[s]]];
      penc[t] <- -reward[t+cumobs[s]] - v[t][3-choice2[t+cumobs[s]]];

      v[t+1][choice2[t+cumobs[s]]]   <- v[t][choice2[t+cumobs[s]]]   + lr[s] * pe[t]; // overwrite chosen value with pe update
      v[t+1][3-choice2[t+cumobs[s]]] <- v[t][3-choice2[t+cumobs[s]]] + lr[s] * penc[t]; // overwrite non-chosen value with penc update
    }
  }
}

generated quantities {
  real<lower=0,upper=1> lr_mu; 
  real<lower=0,upper=10> tau_mu;
  
  real log_lik[nSubjects]; 
  vector[2] v2[maxTrials+1];
  vector[maxTrials] pe2;
  vector[maxTrials] penc2;
  
  lr_mu <- Phi_approx(lr_mu_pr);
  tau_mu <- 10 * Phi_approx(tau_mu_pr);
  
  for (s in 1:nSubjects) {
    log_lik[s] <- 0;
    v2[1] <- rep_vector(0.0,2);
    
    for (t in 1:obs[s+1] ) {
      v2[t][3-choice1[t+cumobs[s]]] <- v2[t][3-choice1[t+cumobs[s]]] + coha[s] * against[t+cumobs[s]] * 4;
      v2[t][choice1[t+cumobs[s]]]   <- v2[t][choice1[t+cumobs[s]]]   + cohw[s] * with[t+cumobs[s]] * 4;
      
      log_lik[s] <- log_lik[s] + categorical_logit_log(choice2[t+cumobs[s]], tau[s] * v2[t]);
      
      pe2[t]   <-  reward[t+cumobs[s]] - v2[t][choice2[t+cumobs[s]]];
      penc2[t] <- -reward[t+cumobs[s]] - v2[t][3-choice2[t+cumobs[s]]];

      v2[t+1][choice2[t+cumobs[s]]]   <- v2[t][choice2[t+cumobs[s]]]   + lr[s] * pe2[t]; // overwrite chosen value with pe update
      v2[t+1][3-choice2[t+cumobs[s]]] <- v2[t][3-choice2[t+cumobs[s]]] + lr[s] * penc2[t]; // overwrite non-chosen value with penc update
    }
  }
}