data {
  int<lower=1> nSubjects;                              // number of subjects
  int<lower=1> nTrials;                                // number of trials 
  int<lower=1,upper=2> choice1[nSubjects,nTrials];     // 1st choices, 1 or 2
  int<lower=1,upper=2> choice2[nSubjects,nTrials];     // 2nd choices, 1 or 2
  real<lower=-1,upper=1> reward[nSubjects,nTrials];    // outcome, 1 or -1
  real<lower=0,upper=1>  with[nSubjects,nTrials];      // No. of with / 4
  real<lower=0,upper=1>  against[nSubjects,nTrials];   // No. of against / 4
}

transformed data {
  vector[2] initV;  // initial values for V
  initV <- rep_vector(0.0,2);    
}

parameters {
  // group-level parameters
  real lr_mu_pr;    // lr_mu before probit
  real tau_mu_pr;   // tau_mu before probit  
  real disc_mu_pr;  // discounting gamma, before probit
  real cra_mu;
  real crw_mu;

  real<lower=0> lr_sd;
  real<lower=0> tau_sd;
  real<lower=0> disc_sd;
  real<lower=0> cra_sd;
  real<lower=0> crw_sd;
  
  // subject-level row parameters, follows norm(0,1), for later Matt Trick
  vector[nSubjects] lr_row;
  vector[nSubjects] tau_row;
  vector[nSubjects] disc_row;
  vector[nSubjects] cra_row;
  vector[nSubjects] crw_row;
}

transformed parameters {
  // subject-level parameters
  vector<lower=0,upper=1>[nSubjects] lr;
  vector<lower=0,upper=10>[nSubjects] tau;
  vector<lower=0,upper=1>[nSubjects] disc;
  vector[nSubjects] cra;
  vector[nSubjects] crw; //?LZ -- how to use truncate in Matt Trick

  // Matt Trick, note that the input of Phi_approx must be 'real' rather than 'vector'
  for (s in 1:nSubjects) {
    lr[s]   <- Phi_approx( lr_mu_pr + lr_sd * lr_row[s] );
    tau[s]  <- 10 * Phi_approx( tau_mu_pr + tau_sd * tau_row[s] );
    disc[s] <- Phi_approx( disc_mu_pr + disc_sd * disc_row[s] );
  }
  cra <- cra_mu + cra_sd * cra_row; // vectorization
  crw <- crw_mu + crw_sd * crw_row;
}

model {
  // define the value and pe vectors
  vector[2] v[nTrials+1];     // values
  vector[nTrials] pe;         // prediction errors
  vector[nTrials] penc;       // pe for the non-chosen choice
  real cr[nSubjects,nTrials]; // cumulative reward

  // hyperparameters
  lr_mu_pr   ~ normal(0,1);
  tau_mu_pr ~ normal(0,1);
  coha_mu    ~ normal(0,1);
  cohw_mu    ~ normal(0,1);

  lr_sd   ~ cauchy(0,5);
  tau_sd ~ cauchy(0,5);
  coha_sd ~ cauchy(0,5);
  cohw_sd ~ cauchy(0,5);
  
  // Matt Trick
  lr_row   ~ normal(0,1);
  tau_row ~ normal(0,1);
  coha_row ~ normal(0,1);
  cohw_row ~ normal(0,1);
  
  // subject loop and trial loop
  for (s in 1:nSubjects) {

    v[1] <- initV;
    
    for (t in 1:nTrials) {
      
      //* re-weight value after have seen the group decisions */
      // also try multiply with the coha/cohw reweight --LZ
      v[t][3-choice1[s,t]] <- v[t][3-choice1[s,t]] + coha[s] * against[s,t];
      v[t][choice1[s,t]]   <- v[t][choice1[s,t]]   + cohw[s] * with[s,t];

      //* compute action probs using built-in softmax function and related to choice data */
      choice2[s,t] ~ categorical_logit( tau[s] * v[t] );

      //* prediction error */
      pe[t]   <-  reward[s,t] - v[t][choice2[s,t]];
      penc[t] <- -reward[s,t] - v[t][3-choice2[s,t]];

      //* value updating (learning) */
      v[t+1][choice2[s,t]]   <- v[t][choice2[s,t]]   + lr[s] * pe[t];
      v[t+1][3-choice2[s,t]] <- v[t][3-choice2[s,t]] + lr[s] * penc[t];
    }
  }
}

generated quantities {
  real<lower=0,upper=1> lr_mu; 
  real<lower=0,upper=10> tau_mu;
  real<lower=0,upper=1> disc_mu;
  
  lr_mu   <- Phi_approx(lr_mu_pr);
  tau_mu  <- 10 * Phi_approx(tau_mu_pr);
  disc_mu <- Phi_approx( disc_mu_pr );
}