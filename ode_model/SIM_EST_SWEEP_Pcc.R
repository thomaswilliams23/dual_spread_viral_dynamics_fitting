# ===========  Load R libraries (install them before loading) ====================
library(rstan)
library(readxl)
library(tidyverse)
library(bridgesampling)
library(deSolve)


# ===========  (Meta) parameter setup ====================
#dispersion parameter
phi_vals = c(1e1, 1e2, 1e3, 1e4, 1e5) 

#folder template
results_base_template = 'RESULTS/fit_cumul_I'
dir.create(results_base_template)
results_template_2 = '/Pcc_'
folder_names_1 = c('01', '05', '09');
folder_names_2 = c('phi_1e1', 'phi_1e2', 'phi_1e3', 'phi_1e4', 'phi_1e5') 

#number of replicates
n_replicates = 10

#file containing the fit
stan_model_name = 'stan_helpers/delay_model_cumul_alpha_beta.stan'

#options to determine which parameters are fitted (FIXED)
params_to_fit = c(1,2)
fixed_params_are_zero = c()

#number of sammples and chains
n_iterates = 2000
n_burnin = 200
n_chains = 4

#assumed size of sheet (for observation model)
sheet_size = 2e5

#observation times
obs_times = c(3,6,9,12,15,18,21,24,27,30)

#initially infected proportion
moi = 0.01


# ===========  Model setup ====================

#parameter values
alpha_array = c(2.28849627e-01, 1.08801559e+00, 1.87383799e+00);
beta_array = c(1.37951410e-06, 7.19841731e-07, 1.36340545e-07);

true_gamma = 3.366934e-01;
true_delta = 8.256588e-02;
true_log_p = 1.409457e+01;
true_c = 4.313531e-01;

param_names = c('alpha', 'beta', 'gamma', 'delta', 'log_p', 'c')

latent_stages = 3
inf_stages = 1
stages = c(latent_stages, inf_stages)


#ICs
T_0 = 1-moi
E_CF_0 = 0
E_CC_0 = 0
I_CF_0 = 0
I_CC_0 = 0
I_r_dead_0 = 0
I_d_0 = moi
I_d_dead_0 = 0
V_0 = 0


#indices
num_compartments = 2*latent_stages + 3*inf_stages + 6

target_ind  = 1
latent_CF_start = target_ind + 1
latent_CC_start = latent_CF_start + latent_stages
inf_CF_start = latent_CC_start + latent_stages
inf_CC_start = inf_CF_start + inf_stages
inf_r_dead_ind = inf_CC_start + inf_stages
inf_d_start = inf_r_dead_ind+1
inf_d_dead_ind = inf_d_start + inf_stages
v_ind = inf_d_dead_ind + 1
CF_cum_ind = v_ind + 1
CC_cum_ind = num_compartments

#IC structure
init_conditions = rep(0, num_compartments)

init_conditions[target_ind] = T_0
init_conditions[latent_CF_start] = E_CF_0
init_conditions[latent_CC_start] = E_CC_0
init_conditions[inf_CF_start] = I_CF_0
init_conditions[inf_CC_start] = I_CC_0
init_conditions[inf_r_dead_ind] = I_r_dead_0
init_conditions[inf_d_start]= I_d_0
init_conditions[inf_d_dead_ind] = I_d_dead_0
init_conditions[v_ind] = V_0



# =========== Outer parameter loop ===================
for (Pcc_ind in 1:length(alpha_array)){
  
  #set up true parameters
  true_alpha = alpha_array[Pcc_ind];
  true_beta = beta_array[Pcc_ind];
  true_params = c(true_alpha, true_beta, true_gamma, true_delta, true_log_p, true_c)
  
  params_orig_run = c(true_params, stages)
  params_orig_run[5] = exp(params_orig_run[5])
  
  
  #output folder
  results_template_this_Pcc = paste(results_base_template, results_template_2, folder_names_1[Pcc_ind], sep="")
  dir.create(results_template_this_Pcc)
  
  
  # ===========  Original model run ====================
  
  #load in model as function
  source("stan_helpers/delay_ode_model_func.R")
  
  orig_model_out <- ode(init_conditions, c(0,obs_times), staged_TIV_catch_infection_modes, params_orig_run)
  orig_model_out = orig_model_out[2:dim(orig_model_out)[1],2:dim(orig_model_out)[2]]
  
  true_inf_cumul_data_raw = orig_model_out[,CF_cum_ind] + orig_model_out[,CC_cum_ind]
  
  
  # ===========  LOOP ====================
  for (noise_ind in 1:(length(phi_vals))){
    
    #output folder
    results_this_noise_level = paste(results_template_this_Pcc, "/", folder_names_2[noise_ind], sep="")
    dir.create(results_this_noise_level)
    results_stem = paste(results_this_noise_level, '/replicate_', sep="")
    
    
    #loop over replicates
    for (replicate in 1:n_replicates){
      
      #output folder
      results_folder = paste(results_stem, as.character(replicate), sep="")
      
      
      # apply observational noise
      obs_cumul_inf_data = rep(0, length(true_inf_cumul_data_raw));
      for (k in 1:length(true_inf_cumul_data_raw)){
        nb_mu = round(sheet_size*true_inf_cumul_data_raw[k]);
        obs_cumul_inf_data[k] = rnbinom(1, mu=nb_mu, size=phi_vals[noise_ind]);
      }
      
      
      # ========== Package data for stan ==================
      
      data_combined = list(N=length(obs_times),
                           time_data=obs_times,
                           T_0=T_0,
                           E_CF_0=E_CF_0,
                           E_CC_0=E_CC_0,
                           I_CF_0=I_CF_0,
                           I_CC_0=I_CC_0,
                           I_r_dead_0=I_r_dead_0,
                           I_d_0=I_d_0,
                           I_d_dead_0=I_d_dead_0,
                           V_0=V_0,
                           t0=0,
                           num_components=num_compartments,
                           target_ind=1,
                           latent_CF_start=target_ind + 1,
                           latent_CC_start=latent_CF_start + latent_stages,
                           inf_CF_start=latent_CC_start + latent_stages,
                           inf_CC_start=inf_CF_start + inf_stages,
                           inf_r_dead_ind=inf_CC_start + inf_stages,
                           inf_d_start=inf_r_dead_ind+1,
                           inf_d_dead_ind=inf_d_start + inf_stages,
                           virus_ind=inf_d_dead_ind + 1,
                           CF_cum_ind=v_ind + 1,
                           CC_cum_ind=num_compartments,
                           sheet_size=sheet_size,
                           phi_param = phi_vals[noise_ind],
                           num_infected_data = obs_cumul_inf_data)
      
      
      # ========== Initialise parameter selection for stan ==================
      #draw randomly from the pdf of all parameters being fit (FIXED as alpha and beta only)
      stan_init = c()
      for (chain in 1:n_chains){
        
        num_params = 6 #hard coded
        
        init_this_chain = list()
        for (prm in 1:num_params){
          
          #draw from priors for all params we are fitting
          if (1 %in% params_to_fit){
            init_this_chain$alpha_param = runif(1, 0, 5)
          }
          if (2 %in% params_to_fit){
            init_this_chain$beta_param = runif(1, 0, 5e-6)
          }
          if (3 %in% params_to_fit){
            init_this_chain$gamma_param = runif(1, 0, 1)
          }
          if (4 %in% params_to_fit){
            init_this_chain$delta_param = exp(rnorm(1, log(true_delta), 0.1*log(true_delta)))
          }
          if (5 %in% params_to_fit){
            init_this_chain$log_p_param = runif(1, 0, 20)
          }
          if (6 %in% params_to_fit){
            init_this_chain$c_param = exp(rnorm(1, log(true_c), 0.1*log(true_c)))
          }
        }
        
        stan_init[[chain]] = init_this_chain
      }
      
      
      # ========== Run Stan ==================
      options(mc.cores = 4)
      
      fitting_results = stan(stan_model_name,
                             data=data_combined,
                             iter=n_iterates,
                             warmup=n_burnin,
                             chains=n_chains,
                             init=stan_init,
                             verbose=TRUE)
      
      
      
      # ========== Extract some stuff ==================
      dir.create(results_folder)
      
      params_to_extract = c("lp__")
      for (param in 1:length(params_to_fit)){
        if (params_to_fit[param]==1){
          params_to_extract[param+1] = 'alpha_param';
        } else if (params_to_fit[param]==2){
          params_to_extract[param+1] = 'beta_param';
        } else if (params_to_fit[param]==3){
          params_to_extract[param+1] = 'gamma_param';
        } else if (params_to_fit[param]==4){
          params_to_extract[param+1] = 'delta_param';
        } else if (params_to_fit[param]==5){
          params_to_extract[param+1] = 'log_p_param';
        } else if (params_to_fit[param]==6){
          params_to_extract[param+1] = 'c_param';
        }
      }
      
      
      #all samples
      posterior_samples_all = rstan::extract(fitting_results, 'pars'=params_to_extract, 'permuted'=FALSE, 'inc_warmup'=TRUE)
      
      
      
      # =============== Save the fitting results and generate outputs in CSV file  ===================
      for (chain in 1:n_chains){
        posterior_sample_array = posterior_samples_all[,chain,];
        fname = paste(results_folder, '/all_samples_chain_', as.character(chain), '.csv', sep="")
        write.csv(posterior_sample_array, file=fname)
      }
      
      write.csv(obs_cumul_inf_data, file=paste(results_folder, '/obs_cumul_inf_data.csv', sep=""))
      
      save.image(file = paste(results_folder, "R_fitting_object.RData", sep="/"))
      
    }
  }
}




