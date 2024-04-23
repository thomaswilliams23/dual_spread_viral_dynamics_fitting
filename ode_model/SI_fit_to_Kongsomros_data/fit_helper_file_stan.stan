// TEMPLATE FOR m LATENT and n INFECTED STAGES

// ensure stage numbers passed in agree with those declared in the function!



functions {
  real[] staged_TEIV_model(real t, real[] y, real[] theta, real[] x_r, int[] x_i){
    // y[1,2...2+m-1,2+m...2+m+n-1,2+m+n,2+m+n+1...2+m+n+1+n-1,2+m+n+1+n] = 
    //     [T, E_r_1, ... E_r_m, I_r_1, ... I_r_n, I_r_dead, I_d_1, ..., I_d_n, V]
    // theta[1,2,3,4,5,6] = [alpha, beta, gamma, delta, p, c]
    
    
    
    //INPUT NUMBER STAGES HERE
    int latent_stages = 3;
    int inf_stages = 1;
    //////////////////////////


    real dydt[3 + latent_stages + 2*inf_stages];

    int target_ind = 1;
    int latent_start = 2;
    int inf_r_start = latent_start + latent_stages;
    int inf_dead_ind = latent_start + latent_stages + inf_stages;
    int inf_d_start = latent_start + latent_stages + inf_stages + 1;
    int virus_ind = latent_start + latent_stages + inf_stages + 1 + inf_stages;
    
    
    // count all infected cells
    real all_inf_cells = 0;
    for (cell_ind in inf_r_start:(inf_r_start+inf_stages-1)){
      all_inf_cells += y[cell_ind];
    }
    for (cell_ind in inf_d_start:(inf_d_start+inf_stages-1)){
      all_inf_cells += y[cell_ind];
    }
    
    
    //target cells
    dydt[target_ind] = -theta[1]*y[target_ind]*(all_inf_cells) - theta[2]*y[target_ind]*y[virus_ind];
    
    //eclipse cells
    dydt[latent_start] = theta[1]*y[target_ind]*(all_inf_cells) + theta[2]*y[target_ind]*y[virus_ind] - latent_stages * theta[3]*y[latent_start];
    for (latent in (latent_start+1):(latent_start+latent_stages-1)){
      dydt[latent] = latent_stages * theta[3] * (y[latent-1] - y[latent]);
    }
    
    //infected (recipient) cells
    dydt[inf_r_start] = latent_stages * theta[3]*y[inf_r_start - 1] - inf_stages * theta[4]*y[inf_r_start];
    for (inf_cell in (inf_r_start+1):(inf_r_start+inf_stages-1)){
      dydt[inf_cell] = inf_stages * theta[4] * (y[inf_cell-1] - y[inf_cell]);
    }
    
    //dead infected cells
    dydt[inf_dead_ind] = inf_stages * theta[4]*y[inf_dead_ind-1];
    
    //infected (donor) cells
    dydt[inf_d_start] = - inf_stages * theta[4]*y[inf_d_start];
    for (inf_cell in (inf_d_start+1):(inf_d_start+inf_stages-1)){
      dydt[inf_cell] = inf_stages * theta[4] * (y[inf_cell-1] - y[inf_cell]);
    }
    
    //virus
    dydt[virus_ind] = exp(theta[5])*all_inf_cells - theta[6]*y[virus_ind];                        //CAUTION!!! CHECK IF USING LOG_P!!
    
    
    
    return dydt;} 
    
}




data {
  int<lower=0> N;
  real time_data[N];
  int<lower=0> num_infected_data[N];
  real T0;
  real E0;
  real I_r_viable_0;
  real I_r_dead_0;
  real I_d_0;
  real V0;
  real t0;
  
  int<lower=0> num_components;                           //ENSURE THIS AGREES WITH VALUES DECLARED ABOVE
  int<lower=0> outind1;
  int<lower=0> outind2;
  
  int target_ind;
  int latent_start;
  int inf_r_start;
  int inf_dead_ind;
  int inf_d_start;
  int virus_ind;
  
  int num_samples;
}




transformed data {
  real x_r[0];
  int x_i[0];
}




parameters {
  real<lower=0> theta[6];
}




transformed parameters {
  real model_output[N,num_components];
  real Y0[num_components];
  real pred_prop_infected[N];
  
  Y0[target_ind] = T0;
  Y0[latent_start] = E0;
  for (i in (latent_start+1):(inf_r_start-1)){
    Y0[i]=0;
  }
  Y0[inf_r_start] = I_r_viable_0;
    for (i in (inf_r_start+1):(inf_dead_ind-1)){
    Y0[i]=0;
  }
  Y0[inf_dead_ind] = I_r_dead_0;
  Y0[inf_d_start] = I_d_0;
  for (i in (inf_d_start+1):(virus_ind-1)){
    Y0[i]=0;
  }
  Y0[virus_ind] = V0;
  
  model_output = integrate_ode_bdf(staged_TEIV_model, Y0, t0, time_data, theta, x_r, x_i);
  
  for (i in 1:N){
    real prop_inf_local = 0;
    for (outind in outind1:outind2){
       prop_inf_local = prop_inf_local + model_output[i,outind]/T0;
    }
    pred_prop_infected[i] = prop_inf_local;
  }
}






model { 
  // priors     (using parameters from Baccam as first guess)
  theta[2] ~ lognormal(log(1e-3*2.167e-3),1);  //beta
  theta[4] ~ lognormal(log(0.217),1);          //delta
  theta[6] ~ lognormal(log(0.3125),1);         //c
  
  //WEAK PRIORS
  theta[1] ~ normal(0, 4);     //alpha
  theta[3] ~ normal(0, 10);    //gamma
  theta[5] ~ normal(0, 10);    //log(p)
  
  
  // BINOMIAL LIKELIHOOD
  num_infected_data ~ binomial(num_samples, pred_prop_infected);
}


