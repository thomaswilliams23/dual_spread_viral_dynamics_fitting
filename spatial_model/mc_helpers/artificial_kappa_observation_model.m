
function obs_kappa_data = artificial_kappa_observation_model(dispersion_param, true_kappa_data)


%applies negative binomial noise to kappa data -- the easy way

artificial_num_cells = 200000;

p_params = dispersion_param./(dispersion_param + artificial_num_cells*true_kappa_data);

obs_kappa_data = nbinrnd(dispersion_param, p_params)/artificial_num_cells;










% %applies negative binomial noise to kappa data
% 
% artificial_num_cells = 200000;
% 
% 
% dispersion_params = 6*base_dispersion_param*true_kappa_data;
% 
% tot_num_neighbours_sampled = 6*artificial_num_cells*obs_fluor_prop_data;
% r_params = base_dispersion_param./(base_dispersion_param + artificial_num_cells*true_kappa_data);
% 
% 
% obs_kappa_data = nbinrnd(dispersion_params, r_params)/tot_num_neighbours_sampled;
% 
% 


