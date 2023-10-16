

%establish defaults
true_params_struct = make_default_param_struct;

%override any changed values
true_params_struct.alpha_param = true_params(1);
true_params_struct.beta_param = true_params(2);
true_params_struct.gamma_param = true_params(3);
true_params_struct.delta_param = true_params(4);
true_params_struct.p = exp(true_params(5));
true_params_struct.c = true_params(6);

true_params_struct.latent_stages = latent_stages;

true_params_struct.moi = moi;

true_params_struct.vis_grid = 0;

true_params_struct.num_species = 1;

true_params_struct.random_initial = 1;
true_params_struct.initialise_in_window = 0;

true_params_struct.show_histogram_of_inf_neighbourhoods = 0;

true_params_struct.final_time = max(obs_times);

true_params_struct.stop_when_all_infected = 0;

true_params_struct.dt = 0.01;