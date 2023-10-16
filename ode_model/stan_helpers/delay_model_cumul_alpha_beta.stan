// TEMPLATE FOR m LATENT and n INFECTED STAGES

// ensure stage numbers passed in agree with those declared in the function!



functions {
  real[] staged_TEIV_model(real t, real[] y, real[] theta, real[] x_r, int[] x_i){
    // y[1,2...2+m-1,2+m...2+m+n-1,2+m+n,2+m+n+1...2+m+n+1+n-1,2+m+n+1+n] = 
    //     [T, E_r_1, ... E_r_m, I_r_1, ... I_r_n, I_r_dead, I_d_1, ..., I_d_n, V]
    // theta[1,2,3,4,5,6] = [alpha, beta, gamma, delta, p, c]
    
    
    
    //HARD CODED number of delay stages and fixed parameters other than alpha and beta
    int latent_stages = 3;
    int inf_stages = 1;

    real gamma_param = 3.366934e-01;
    real delta_param = 8.256588e-02;
    real log_p_param = 1.409457e+01;
    real c_param = 4.313531e-01;
    //////////////////////////


    real dydt[2*latent_stages + 3*inf_stages + 6];

    int num_compartments = 2*latent_stages + 3*inf_stages + 6;

    int target_ind  = 1;
    int latent_CF_start = target_ind + 1;
    int latent_CC_start = latent_CF_start + latent_stages;
    int inf_CF_start = latent_CC_start + latent_stages;
    int inf_CC_start = inf_CF_start + inf_stages;
    int inf_r_dead_ind = inf_CC_start + inf_stages;
    int inf_d_start = inf_r_dead_ind+1;
    int inf_d_dead_ind = inf_d_start + inf_stages;
    int v_ind = inf_d_dead_ind + 1;
    int CF_cum_ind = v_ind + 1;
    int CC_cum_ind = num_compartments;
    
    
    // count all infected cells
    real all_inf_cells = 0;
    for (cell_ind in inf_CF_start:(inf_CF_start+2*inf_stages-1)){
      all_inf_cells += y[cell_ind];
    }
    for (cell_ind in inf_d_start:(inf_d_start+inf_stages-1)){
      all_inf_cells += y[cell_ind];
    }
    
    
    //target cells
    dydt[target_ind] = -theta[1]*y[target_ind]*(all_inf_cells) - theta[2]*y[target_ind]*y[v_ind];
    
    //eclipse CF cells
    dydt[latent_CF_start] = theta[2]*y[target_ind]*y[v_ind] - latent_stages * gamma_param *y[latent_CF_start];
    for (latent in (latent_CF_start+1):(latent_CF_start+latent_stages-1)){
      dydt[latent] = latent_stages * gamma_param * (y[latent-1] - y[latent]);
    }

    //eclipse CC cells
    dydt[latent_CC_start] = theta[1]*y[target_ind]*(all_inf_cells) - latent_stages * gamma_param*y[latent_CC_start];
    for (latent in (latent_CC_start+1):(latent_CC_start+latent_stages-1)){
      dydt[latent] = latent_stages * gamma_param * (y[latent-1] - y[latent]);
    }
    
    //infected (recipient, CF) cells
    dydt[inf_CF_start] = latent_stages * gamma_param * y[latent_CF_start+latent_stages-1] - inf_stages * delta_param *y[inf_CF_start];
    for (inf_cell in (inf_CF_start+1):(inf_CF_start+inf_stages-1)){
      dydt[inf_cell] = inf_stages * delta_param * (y[inf_cell-1] - y[inf_cell]);
    }

    //infected (recipient, CC) cells
    dydt[inf_CC_start] = latent_stages * gamma_param * y[latent_CC_start+latent_stages-1] - inf_stages * delta_param *y[inf_CC_start];
    for (inf_cell in (inf_CC_start+1):(inf_CC_start+inf_stages-1)){
      dydt[inf_cell] = inf_stages * delta_param * (y[inf_cell-1] - y[inf_cell]);
    }
    
    //dead infected cells
    dydt[inf_r_dead_ind] = inf_stages * delta_param * (y[inf_CF_start+inf_stages-1] + y[inf_CC_start+inf_stages-1]);
    
    //infected (donor) cells
    dydt[inf_d_start] = - inf_stages * delta_param *y[inf_d_start];
    for (inf_cell in (inf_d_start+1):(inf_d_start+inf_stages-1)){
      dydt[inf_cell] = inf_stages * delta_param * (y[inf_cell-1] - y[inf_cell]);
    }

    //dead infected cells
    dydt[inf_d_dead_ind] = inf_stages * delta_param * y[inf_d_start+inf_stages-1];
    
    //virus
    dydt[v_ind] = exp(log_p_param)*all_inf_cells - c_param*y[v_ind];                      //CAUTION!!! CHECK IF USING LOG_P!!
    

    //cumulative CF cells
    dydt[CF_cum_ind] = latent_stages * gamma_param * y[latent_CF_start+latent_stages-1];

    //cumulative CC cells
    dydt[CC_cum_ind] = latent_stages * gamma_param * y[latent_CC_start+latent_stages-1];
    
    
    return dydt;} 
    
}




data {
  int<lower=0> N;
  real time_data[N];
  int<lower=0> num_infected_data[N];
  real T_0;
  real E_CF_0;
  real E_CC_0;
  real I_CF_0;
  real I_CC_0;
  real I_r_dead_0;
  real I_d_0;
  real I_d_dead_0;
  real V_0;
  real t0;
  
  int<lower=0> num_components;
  
  int target_ind;
  int latent_CF_start;
  int latent_CC_start;
  int inf_CF_start;
  int inf_CC_start;
  int inf_r_dead_ind;
  int inf_d_start;
  int inf_d_dead_ind;
  int virus_ind;
  int CF_cum_ind;
  int CC_cum_ind;
  
  int sheet_size;
  real phi_param;
}




transformed data {
  real x_r[0];
  int x_i[0];
}




parameters {
  real<lower=0, upper=5> alpha_param;
  real<lower=0, upper=1e-5> beta_param;
}




transformed parameters {
  real model_output[N,num_components];
  real Y0[num_components];
  real pred_num_infected[N];

  real theta[2];

  #initial conditions
  Y0[target_ind] = T_0;
  Y0[latent_CF_start] = E_CF_0;
  for (i in (latent_CF_start+1):(latent_CC_start-1)){
    Y0[i]=0;
  }
  Y0[latent_CC_start] = E_CC_0;
  for (i in (latent_CC_start+1):(inf_CF_start-1)){
    Y0[i]=0;
  }
  Y0[inf_CF_start] = I_CF_0;
    for (i in (inf_CF_start+1):(inf_CC_start-1)){
    Y0[i]=0;
  }
    Y0[inf_CC_start] = I_CC_0;
    for (i in (inf_CC_start+1):(inf_r_dead_ind-1)){
    Y0[i]=0;
  }
  Y0[inf_r_dead_ind] = I_r_dead_0;
  Y0[inf_d_start] = I_d_0;
    for (i in (inf_d_start+1):(inf_d_dead_ind-1)){
    Y0[i]=0;
  }
  Y0[inf_d_dead_ind] = I_d_dead_0;
  Y0[virus_ind] = V_0;
  Y0[CF_cum_ind] = 0;
  Y0[CC_cum_ind] = 0;


  theta[1] = alpha_param;
  theta[2] = beta_param;
  
  
  // run model
  model_output = integrate_ode_bdf(staged_TEIV_model, Y0, t0, time_data, theta, x_r, x_i);
  
  
  // extract number of fluorescent cells at observation times
  for (i in 1:N){
    pred_num_infected[i] = sheet_size * (model_output[i,CC_cum_ind] + model_output[i,CF_cum_ind]);
  }
}






model { 
  // priors
  alpha_param ~ uniform(0,5);
  beta_param ~ uniform(0,1e-5);
  
  // negative binomial likelihood
  num_infected_data ~ neg_binomial_2(pred_num_infected, phi_param);
}


