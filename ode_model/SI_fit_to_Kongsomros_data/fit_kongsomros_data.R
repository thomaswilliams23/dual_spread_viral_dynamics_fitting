rm(list=ls())  # clear memory
# ===========  Load R libraries (install them before loading) ====================
library(rstan)
library(readxl)
library(tidyverse)
library(bayesplot)
library(bridgesampling)

# ===========  num stages =============== 
latent_stages = 3
inf_stages = 1
num_components = 3 + latent_stages + 2*inf_stages

outind1 = 2 + latent_stages #start of range for output 
outind2 = outind1 + inf_stages #end of range for output

target_ind = 1;
latent_start = 2;
inf_r_start = latent_start + latent_stages;
inf_dead_ind = latent_start + latent_stages + inf_stages;
inf_d_start = latent_start + latent_stages + inf_stages + 1;
virus_ind = latent_start + latent_stages + inf_stages + 1 + inf_stages;


# ===========  load data ===============
kongsomros_data= read_excel("Kongsomros_data.xlsx",1)

# read out vectors
time_data = kongsomros_data$time_hours[1:4]
percent_infected_data = kongsomros_data$percent_infected_h1n1_coculture[1:4]

#set assumed number of samples (controls noise of likelihood)
num_samples = 100;

#convert percentage data to number of infected cells
num_infected_data = round((num_samples/100) * percent_infected_data)


# merge data and fixed parameter values in a list
data_combined = list(N = length(time_data),
                     time_data = time_data,
                     num_infected_data = num_infected_data,
                     T0 = 0.5,
                     E0 = 0,
                     I_r_viable_0 = 0,
                     I_r_dead_0 = 0,
                     I_d_0 = 0.5,
                     V0 = 0,
                     t0 = 0,
                     num_components = num_components,
                     outind1 = outind1,
                     outind2 = outind2,
                     target_ind = 1,
                     latent_start = 2,
                     inf_r_start = latent_start + latent_stages,
                     inf_dead_ind = latent_start + latent_stages + inf_stages,
                     inf_d_start = latent_start + latent_stages + inf_stages + 1,
                     virus_ind = latent_start + latent_stages + inf_stages + 1 + inf_stages,
                     num_samples = num_samples)


set.seed(202004)   # set random seed for reproducibility

# set the initial conditions to be a list
init1 = list(
  theta = c(0.2, 1e-3*2.167e-3, 0.167, 0.217, log(384), 0.3125))

init2 = list(
  theta = c(runif(1,0.1,0.5), runif(1,1e-3*1e-3,1e-3*5e-3), runif(1,0.05,0.5), runif(1,0.1,0.5), runif(1,log(100),log(1000)), runif(1,0.1,1)))


# set number of cores to use
options(mc.cores=2)

# model fitting
fitting_results = stan("fit_helper_file_stan.stan",
                  data = data_combined,
                  pars = c("theta"),
                  seed = 90127,  # set random seed for reproducibility
                  iter = 5000,
                  chains = 2,
                  init = list(init1,init2),
                  warmup = 1000,
                  control = list(adapt_delta = 0.99, max_treedepth = 20)) 

# extract posterior samples for selected parameters
posterior_samples_all = rstan::extract(fitting_results, pars = c("theta"), inc_warmup = TRUE, permuted = FALSE)
posterior_samples_merged_after_burnin = rstan::extract(fitting_results, pars = c("theta"))
  

# =============== Save outputs in CSV file  ===================
write.csv(posterior_samples_merged_after_burnin, file="posterior_samples.csv")
