# this mode tries coha/cohw with matt trick without truncation
data {
  int<lower=1> nSubjects;                              // number of subjects
  int<lower=1> nTrials;                                // number of trials 
  int<lower=1,upper=2> choice1[nSubjects,nTrials];     // 1st choices, 1 or 2
  int<lower=1,upper=2> choice2[nSubjects,nTrials];     // 2nd choices, 1 or 2
  real<lower=-1,upper=1> reward[nSubjects,nTrials];    // outcome, 1 or -1
  real<lower=0,upper=1>  with[nSubjects,nTrials];      // No. of with / 4
  real<lower=0,upper=1>  against[nSubjects,nTrials];   // No. of against / 4
}

parameters {
  // group-level parameters
  real lr_mu_pr;    // lr_mu before probit
  real tau_mu_pr;   // tau_mu before probit  
  real coha_mu;
  real cohw_mu;

  real<lower=0> lr_sd;
  real<lower=0> tau_sd;
  real<lower=0, upper=3> coha_sd;
  real<lower=0, upper=3> cohw_sd;
  
  // subject-level raw parameters, follows norm(0,1), for later Matt Trick
  vector[nSubjects] lr_raw;
  vector[nSubjects] tau_raw;
  vector[nSubjects] coha_raw;
  vector[nSubjects] cohw_raw;
}

transformed parameters {
  // subject-level parameters
  vector<lower=0,upper=1>[nSubjects] lr;
  vector<lower=0,upper=2.5>[nSubjects] tau;
  vector[nSubjects] coha;
  vector[nSubjects] cohw;

  // Matt Trick, note that the input of Phi_approx must be 'real' rather than 'vector'
  for (s in 1:nSubjects) {
    lr[s]   <- Phi_approx( lr_mu_pr + lr_sd * lr_raw[s] );
    tau[s]  <- Phi_approx( tau_mu_pr + tau_sd * tau_raw[s] ) * 2.5;
  }
  coha <- coha_mu + coha_sd * coha_raw;
  cohw <- cohw_mu + cohw_sd * cohw_raw;
}

model {
  // define the value and pe vectors
  vector[2] v[nTrials+1];    //  values
  vector[nTrials] pe;         // prediction errors
  vector[nTrials] penc;       // pe for the non-chosen choice

  // hyperparameters
  lr_mu_pr  ~ normal(0,1);
  tau_mu_pr ~ normal(0,1);
  coha_mu   ~ normal(0,0.5); 
  cohw_mu   ~ normal(0,0.5);

  lr_sd   ~ cauchy(0,5);
  tau_sd  ~ cauchy(0,5);
//   coha_sd ~ cauchy(0,3);
//   cohw_sd ~ cauchy(0,3);
  
  // Matt Trick
  lr_raw   ~ normal(0,1);
  tau_raw  ~ normal(0,1);
  coha_raw ~ normal(0,0.5);
  cohw_raw ~ normal(0,0.5);
  
  // subject loop and trial loop
  for (s in 1:nSubjects) {
    v[1] <- rep_vector(0.0,2);
    
    //print("s=", s, "=======================================================================");
    //print("lr=", lr[s], "tau=", tau[s], ",coha=", coha[s], ", cohw=", cohw[s]);
    
    for (t in 1:nTrials) {
      
      //print("t=", t, "---------------------------------------------------------------------");
      //print("choice1=", choice1[s,t], ", choice2=", choice2[s,t], ", with=", with[s,t], ", against=", against[s,t]);
      
      //print("before social", v[t]);
      
      v[t][3-choice1[s,t]] <- v[t][3-choice1[s,t]] + coha[s] * against[s,t];
      v[t][choice1[s,t]]   <- v[t][choice1[s,t]]   + cohw[s] * with[s,t];
      
      //print("after social", v[t], ", coha[s]=", coha[s], ", cohw[s]=", cohw[s], ", tau[s]=", tau[s], ", tau*v", tau[s]*v[t]);
      
      choice2[s,t] ~ categorical_logit( tau[s] * v[t] );
      
      pe[t]   <-  reward[s,t] - v[t][choice2[s,t]];
      penc[t] <- -reward[s,t] - v[t][3-choice2[s,t]];
      
      //print("pe=", pe[t], ", penc=", penc[t]);

      v[t+1][choice2[s,t]]   <- v[t][choice2[s,t]]   + lr[s] * pe[t]; 
      v[t+1][3-choice2[s,t]] <- v[t][3-choice2[s,t]] + lr[s] * penc[t];
      
      //print("after update", v[t], v[t+1]);
    }
  }
}

generated quantities {
  real<lower=0,upper=1> lr_mu; 
  real<lower=0,upper=2.5> tau_mu;

  real log_lik[nSubjects]; 
  vector[2] v2[nTrials+1];
  vector[nTrials] pe2;
  vector[nTrials] penc2;
  
  lr_mu  <- Phi_approx(lr_mu_pr);
  tau_mu <- Phi_approx(tau_mu_pr) * 2.5;


  for (s in 1:nSubjects) {
    log_lik[s] <- 0;
    v2[1] <- rep_vector(0.0,2);
    
    for (t in 1:nTrials) {
      v2[t][3-choice1[s,t]] <- v2[t][3-choice1[s,t]] + coha[s] * against[s,t];
      v2[t][choice1[s,t]]   <- v2[t][choice1[s,t]]   + cohw[s] * with[s,t];
      
      log_lik[s] <- log_lik[s] + categorical_logit_log(choice2[s,t], tau[s] * v2[t]);
      
      pe2[t]   <-  reward[s,t] - v2[t][choice2[s,t]];
      penc2[t] <- -reward[s,t] - v2[t][3-choice2[s,t]];

      v2[t+1][choice2[s,t]]   <- v2[t][choice2[s,t]]   + lr[s] * pe2[t];
      v2[t+1][3-choice2[s,t]] <- v2[t][3-choice2[s,t]] + lr[s] * penc2[t];
    }
  }
}