function output  = ABC_rejection_sampler_fluoro_only(init_params_as_struct, params_to_fit, ...
    moi, sheet_size, sim_times, indices_of_obs_times, obs_inf_data, acceptance_rate, ...
    num_particles, prior_dists)
%naive ABC rejection sampler as implemented in Kypraios et al.

%args: parameters for model simulation, reference statistics, epsilon value
% required number of particles, prior functions



%initialise sample vector, sim count and particle index
n_params_to_fit = length(params_to_fit);


%initialise parameter structures
params_as_struct = init_params_as_struct;
params_as_array = [params_as_struct.alpha_param, params_as_struct.beta_param, ...
    params_as_struct.gamma_param,params_as_struct.delta_param, ...
    log(params_as_struct.p), params_as_struct.c];


%initialise structures for run
num_samples_to_run = round(num_particles/acceptance_rate);

all_samples = zeros(n_params_to_fit, num_samples_to_run);
err_vec = zeros(1, num_samples_to_run);
accepted_Pcc = zeros(1, num_samples_to_run);
accepted_ppc = zeros(length(sim_times), num_samples_to_run);


%parallel loop
parfor sim_ind = 1:num_samples_to_run

    %draw proposal sample from priors
    prop_params_as_struct = params_as_struct;
    prop_params_as_array = params_as_array;
    for prm = params_to_fit
        prior_dist = prior_dists{prm};
        prop_params_as_array(prm) = random(prior_dist,1);
    end

    prop_params_as_struct.alpha_param = prop_params_as_array(1);
    prop_params_as_struct.beta_param = prop_params_as_array(2);
    prop_params_as_struct.gamma_param = prop_params_as_array(3);
    prop_params_as_struct.delta_param = prop_params_as_array(4);
    prop_params_as_struct.p = exp(prop_params_as_array(5));
    prop_params_as_struct.c = prop_params_as_array(6);


    %call model with proposed sample
    [model_sim] = single_run_as_func(sim_times, prop_params_as_struct);

    %extract fluorescence data
    model_pred_cumul_all = model_sim.prop_infected  * (1-moi)*sheet_size;
    model_pred_cumul_inf = model_pred_cumul_all(indices_of_obs_times);

    
    %save into output structures
    for prm_ind = 1:n_params_to_fit
        all_samples(prm_ind, sim_ind) = prop_params_as_array(params_to_fit(prm_ind));
    end

    %error compared to data
    err_vec(sim_ind) = data_to_sim_dist(obs_inf_data, model_pred_cumul_inf);

    %actual Pcc of proposal sim
    accepted_Pcc(sim_ind) = model_sim.net_CC(end)/(model_sim.net_CC(end) + model_sim.net_CF(end));

    %full fluorescence time series for proposal sim
    accepted_ppc(:,sim_ind) = model_pred_cumul_all;

end


%find the kth best errors
[~,best_inds] = mink(err_vec,num_particles);


%crop to best samples
all_samples = all_samples(:, best_inds);
err_vec = err_vec(best_inds);
accepted_Pcc = accepted_Pcc(best_inds);
accepted_ppc = accepted_ppc(:,best_inds);




%output formatting
output.sim_count = num_samples_to_run;
output.samples = all_samples;
output.err_vec = err_vec;
output.accepted_Pcc = accepted_Pcc;
output.accepted_ppc = accepted_ppc;


