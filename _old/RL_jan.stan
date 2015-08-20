data {
  
  int nSubjects;
  int nTrials;
  row_vector<lower=1,upper=2>[nTrials] choice[nSubjects]; // 2nd choices here
  row_vector<lower=-1,upper=1>[nTrials] reward[nSubjects];

}

transformed data {

  real A; // prior for alpha parameter
  real B; // prior for beta parameter
  row_vector[2] initV; // initial values for v

  A <- 1.0;
  B <- 1.0;
  initV <- rep_row_vector( 0.0 , 2 );

}

parameters {

  real<lower=0,upper=1> lr_mu;
  real<lower=0,upper=1> lr_kappa0;
  real<lower=0,upper=1> temp_mu;
  real<lower=0,upper=1> temp_kappa0;

  vector<lower=0,upper=1>[nSubjects] lr;
  vector<lower=0,upper=1>[nSubjects] temp0;

}

transformed parameters {

  real<lower=0> lr_kappa;
  real<lower=0> lr_a;
  real<lower=0> lr_b;

  real<lower=0> temp_kappa;
  real<lower=0> temp_a;
  real<lower=0> temp_b;

  vector<lower=0,upper=1>[nSubjects] temp;

  lr_kappa <- 3.0 / (lr_kappa0 * lr_kappa0) - 1.0;
  lr_a <- lr_mu * lr_kappa;
  lr_b <- (1-lr_mu) * lr_kappa;

  temp_kappa <- 3.0 / (temp_kappa0 * temp_kappa0) - 1.0;
  temp_a <- temp_mu * temp_kappa;
  temp_b <- (1-temp_mu) * temp_kappa;

  for (s in 1:nSubjects) temp[s] <- 20.0 * temp0[s];

}

model {

  row_vector[2] v[nTrials+1]; // values
  vector[nTrials] pe; // prediction errors

  lr_mu ~ beta(A,B);
  lr_kappa0 ~ beta(A,B);

  temp_mu ~ beta(A,B);
  temp_kappa0 ~ beta(A,B);

  for (s in 1:nSubjects) {

    lr[s] ~ beta(lr_a,lr_b);
    temp0[s] ~ beta(temp_a,temp_b);

    v[1] <- initV;

    for (t in 1:nTrials) {

      /* compute action probs using built-in softmax function and related to choice data */
      choice[s,t] ~ categorical(softmax(temp[s] * v[t]));

      /* prediction error */
      pe[t] <- reward[t] - v[t,choice[s,t]];

      /* value updating (learning) */
      v[t+1] <- v[t]; // make a copy of current value into t+1
      v[t+1,choice[s,t]] <- v[t,choice[s,t]] + lr[s] * pe[t]; // overwrite chosen value with pe update

    }
  }
}