staged_TIV_catch_infection_modes <- function(t, Y, params){
  
  
  alpha_prm = params[1];
  beta_prm = params[2];
  gamma_prm = params[3];
  delta_prm = params[4];
  p_prm = params[5];
  c_prm = params[6];
  
  latent_stages = params[7];
  inf_stages = params[8];
  
  
  
  #indices
  num_compartments = length(Y);
  
  target_ind  = 1;
  latent_CF_start = target_ind + 1;
  latent_CC_start = latent_CF_start + latent_stages;
  inf_CF_start = latent_CC_start + latent_stages;
  inf_CC_start = inf_CF_start + inf_stages;
  inf_r_dead_ind = inf_CC_start + inf_stages;
  inf_d_start = inf_r_dead_ind+1;
  inf_d_dead_ind = inf_d_start + inf_stages;
  v_ind = inf_d_dead_ind + 1;
  
  tot_CF_inf_ind = num_compartments - 1;
  tot_CC_inf_ind = num_compartments;
  
  
  
  #======================function=======================
  vec_out = rep(0,num_compartments);
  
  #target cells
  vec_out[target_ind] = - alpha_prm * Y[target_ind] * (sum(Y[inf_CF_start:(inf_CF_start+inf_stages-1)]) + 
                                                         sum(Y[inf_CC_start:(inf_CC_start+inf_stages-1)]) + 
                                                         sum(Y[inf_d_start:(inf_d_start+inf_stages-1)])) -
    beta_prm * Y[target_ind] * Y[v_ind];
  
  
  
  #start of CF latent cycle
  vec_out[latent_CF_start] = beta_prm * Y[target_ind] * Y[v_ind] -
    latent_stages * gamma_prm * Y[latent_CF_start];         
  
  #rest of CF latent cycle 
  if (latent_stages>1){
    vec_out[(latent_CF_start+1):(latent_CF_start+latent_stages-1)] = 
      latent_stages * gamma_prm * (Y[latent_CF_start:(latent_CF_start+latent_stages-2)] - 
                                     Y[(latent_CF_start+1):(latent_CF_start+latent_stages-1)]); 
  }
  
  
  #start of CC latent cycle
  vec_out[latent_CC_start] = alpha_prm * Y[target_ind] * (sum(Y[inf_CF_start:(inf_CF_start+inf_stages-1)]) +
                                                            sum(Y[inf_CC_start:(inf_CC_start+inf_stages-1)]) +
                                                            sum(Y[inf_d_start:(inf_d_start+inf_stages-1)])) -
    latent_stages * gamma_prm * Y[latent_CC_start];              
  
  #rest of CC latent cycle 
  if (latent_stages>1){
    vec_out[(latent_CC_start+1):(latent_CC_start+latent_stages-1)] = 
      latent_stages * gamma_prm * (Y[latent_CC_start:(latent_CC_start+latent_stages-2)] - 
                                     Y[(latent_CC_start+1):(latent_CC_start+latent_stages-1)]); 
  }
  
  
  #start of infected cycle (recipient, cell free)
  vec_out[inf_CF_start] = latent_stages * gamma_prm * Y[latent_CF_start+latent_stages-1]-
    inf_stages * delta_prm * Y[inf_CF_start];
  
  #rest of infected cycle (recipient, cell free)
  if (inf_stages>1){
    vec_out[(inf_CF_start+1):(inf_CF_start+inf_stages-1)] = 
      inf_stages * delta_prm * (Y[inf_CF_start:(inf_CF_start+inf_stages-2)] -
                                  Y[(inf_CF_start+1):(inf_CF_start+inf_stages-1)]);  
  }
  
  
  #start of infected cycle (recipient, cell-cell)
  vec_out[inf_CC_start] = latent_stages * gamma_prm * Y[latent_CC_start+latent_stages-1]-
    inf_stages * delta_prm * Y[inf_CC_start];
  
  #rest of infected cycle (recipient, cell-cell)
  if (inf_stages>1){
    vec_out[(inf_CC_start+1):(inf_CC_start+inf_stages-1)] = 
      inf_stages * delta_prm * (Y[inf_CC_start:(inf_CC_start+inf_stages-2)] - 
                                  Y[(inf_CC_start+1):(inf_CC_start+inf_stages-1)]);  
  }
  
  
  #dead infected cells (recipient)
  vec_out[inf_r_dead_ind] = inf_stages * delta_prm * (Y[inf_CF_start+inf_stages-1] + Y[inf_CC_start + inf_stages-1]); 
  
  
  
  #start of second infected cycle (donor)
  vec_out[inf_d_start] = -inf_stages * delta_prm * Y[inf_d_start]; 
  
  #rest of second infected cycle (donor)
  if (inf_stages>1){
    vec_out[(inf_d_start+1):(inf_d_start+inf_stages-1)] = 
      inf_stages * delta_prm * (Y[inf_d_start:(inf_d_start+inf_stages-2)] - 
                                  Y[(inf_d_start+1):(inf_d_start+inf_stages-1)]);
  }
  
  
  
  #dead donor cells
  vec_out[inf_d_dead_ind] = inf_stages * delta_prm * Y[inf_d_start+inf_stages-1];
  
  
  
  #virus
  vec_out[v_ind] = p_prm * (sum(Y[inf_CF_start:(inf_CF_start+inf_stages-1)]) + 
                              sum(Y[inf_CC_start:(inf_CC_start+inf_stages-1)]) + 
                              sum(Y[inf_d_start:(inf_d_start+inf_stages-1)])) - c_prm * Y[v_ind];
  
  
  
  
  #TOTAL INFECTIONS (no death)
  vec_out[tot_CF_inf_ind] = latent_stages * gamma_prm * Y[latent_CF_start+latent_stages-1];
  vec_out[tot_CC_inf_ind] = latent_stages * gamma_prm * Y[latent_CC_start+latent_stages-1];
  
  return(list(vec_out))}