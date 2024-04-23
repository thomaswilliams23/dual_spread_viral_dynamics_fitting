

function obs_fluor_data = artificial_fluor_observation_model(dispersion_param, fluor_prop_data)

%applies negative binomial noise to fluorescence data

artificial_num_cells = 200000;

p_params = dispersion_param./(dispersion_param + artificial_num_cells*fluor_prop_data);

obs_fluor_data = nbinrnd(dispersion_param, p_params)/artificial_num_cells;