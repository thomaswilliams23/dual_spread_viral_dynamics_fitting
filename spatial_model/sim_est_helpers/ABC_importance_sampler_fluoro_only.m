function output = ABC_importance_sampler_fluoro_only(init_params_as_struct, params_to_fit, ...
    moi, sheet_size, sim_times, indices_of_obs_times, obs_inf_data, eps_val, num_particles, ...
    prior_funcs, prev_samples, prev_weights)

%args: parameters for model simulation, reference statistics, epsilon
%value, number of particles, previous sample vector, previous weight
%vector



%define standard dev of noise kernel (comes from Beaumont et al. 2009)
%V=cov(prev_samples');
mu_w = sum(prev_weights'.*prev_samples');
V = num_particles/(num_particles-1) * ...
    (prev_weights'.*(prev_samples'-mu_w))' * (prev_samples'-mu_w);


%initialise output structures
valid_samples = zeros(length(params_to_fit),num_particles);
valid_weights = zeros(1,num_particles);
err_vec = zeros(1,num_particles);

accepted_Pcc = zeros(1,num_particles);
accepted_peaks = zeros(1,num_particles);
accepted_ppc = zeros(length(sim_times),num_particles);

n_particles_found = 0;

n_params_to_fit = length(params_to_fit);


%count simulations
sim_count = 0;


%while particle index less than required number of particles
while n_particles_found<num_particles

    %number of sims to pass to CPUs
    MIN_SIMS_TO_CALL = 10;
    sims_remaining = num_particles-n_particles_found;
    sims_this_batch = max(sims_remaining,MIN_SIMS_TO_CALL);
    sim_count = sim_count + sims_this_batch;

    %parfor data structures
    all_samples = zeros(n_params_to_fit,sims_this_batch);
    all_weights = zeros(1,sims_this_batch);
    samples_are_valid = zeros(1,sims_this_batch);
    err_this_batch = zeros(1,sims_this_batch);
    Pcc_this_batch = zeros(1,sims_this_batch);
    peak_this_batch = zeros(1,sims_this_batch);
    ppc_this_batch = zeros(length(sim_times),sims_this_batch);

    %valid sample region indicator
    in_valid_region = zeros(1,sims_this_batch);


    %parfor
    parfor sim_ind = 1:sims_this_batch

        %draw from old sample vector with replacement
        rnd_draw = rand;
        cumul_weight = 0;
        smpl_ind=0;
        while cumul_weight<rnd_draw
            smpl_ind = smpl_ind + 1;
            cumul_weight = cumul_weight + prev_weights(smpl_ind);
        end


        %perturb sample using gaussian noise
        prop_params_to_fit = mvnrnd(prev_samples(:,smpl_ind)',2*V);


        %if proposed sample in valid region (i.e. if priors>0)
        in_valid_region(sim_ind) = 1;
        prob_eps = 1e-10;
        for prm = params_to_fit
            prior_this_prm = prior_funcs{prm};
            if prior_this_prm(prop_params_to_fit(prm))<prob_eps
                in_valid_region(sim_ind) = 0;
            end
        end

        if in_valid_region(sim_ind)


            %call model with proposed sample
            prop_params_as_struct = init_params_as_struct;
            for prm1 = 1:n_params_to_fit
                if params_to_fit(prm1)==1
                    prop_params_as_struct.alpha_param = prop_params_to_fit(prm1);
                elseif params_to_fit(prm1)==2
                    prop_params_as_struct.beta_param = prop_params_to_fit(prm1);
                elseif params_to_fit(prm1)==3
                    prop_params_as_struct.delta_param = prop_params_to_fit(prm1);
                elseif params_to_fit(prm1)==4
                    prop_params_as_struct.gamma_param = prop_params_to_fit(prm1);
                elseif params_to_fit(prm1)==5
                    prop_params_as_struct.p = exp(prop_params_to_fit(prm1));
                elseif params_to_fit(prm1)==6
                    prop_params_as_struct.c = prop_params_to_fit(prm1);
                end
            end
            prop_params_as_array = [prop_params_as_struct.alpha_param, prop_params_as_struct.beta_param, ...
                prop_params_as_struct.gamma_param, prop_params_as_struct.delta_param, ...
                log(prop_params_as_struct.p), prop_params_as_struct.c];


            [model_sim] = single_run_as_func(sim_times, prop_params_as_struct);
            model_pred_cumul_all = model_sim.prop_infected * (1-moi)*sheet_size;
            model_pred_cumul_inf = model_pred_cumul_all(indices_of_obs_times);


            %compute error from data
            err_this_sample=data_to_sim_dist(obs_inf_data, model_pred_cumul_inf);

            %if error within epsilon of reference
            if (err_this_sample<eps_val)

                %draw was valid
                samples_are_valid(sim_ind)=1;

                %save error and full fluorescence time series (posterior
                %predictive check - PCC), as well as actual Pcc
                err_this_batch(sim_ind)=err_this_sample;
                ppc_this_batch(:,sim_ind) = model_pred_cumul_all;
                Pcc_this_batch(sim_ind) = model_sim.net_CC(end)/(model_sim.net_CC(end) + model_sim.net_CF(end));

                %compute peak time
                [~,max_ind] = max(model_sim.net_I);
                peak_this_batch(sim_ind) = sim_times(max_ind);

                %save samples and weights
                all_weights(sim_ind)=1;
                for prm_ind = 1:n_params_to_fit
                    all_samples(prm_ind, sim_ind) = prop_params_as_array(params_to_fit(prm_ind));

                    prior_this_prm = prior_funcs{prm_ind};
                    all_weights(sim_ind) = all_weights(sim_ind)*...
                        prior_this_prm(prop_params_as_array(params_to_fit(prm_ind)));
                end

                %weights denominator
                weight_denom = 0;
                for smpl = 1:num_particles
                    weight_denom = weight_denom + prev_weights(smpl)*...
                        mvnpdf(prop_params_as_array(params_to_fit),prev_samples(:,smpl)',2*V);
                end

                %update weight
                all_weights(sim_ind) = all_weights(sim_ind)/weight_denom;

            end

        end
    end

    %total valid number of particles
    tot_valid_particles = n_particles_found+sum(samples_are_valid);
    valid_samples_this_batch = all_samples(:,(samples_are_valid>0));
    valid_weights_this_batch = all_weights(:,(samples_are_valid>0));
    valid_err_this_batch = err_this_batch((samples_are_valid>0));
    valid_Pcc_this_batch = Pcc_this_batch((samples_are_valid>0));
    valid_peak_this_batch = peak_this_batch(:,(samples_are_valid>0));
    valid_ppc_this_batch = ppc_this_batch(:,(samples_are_valid>0));

    %unpack into arrays
    if tot_valid_particles<num_particles
        valid_samples(:,(n_particles_found+1):tot_valid_particles) = valid_samples_this_batch;
        valid_weights((n_particles_found+1):tot_valid_particles) = valid_weights_this_batch;
        err_vec((n_particles_found+1):tot_valid_particles) = valid_err_this_batch;
        accepted_Pcc((n_particles_found+1):tot_valid_particles) = valid_Pcc_this_batch;
        accepted_peaks((n_particles_found+1):tot_valid_particles) = valid_peak_this_batch;
        accepted_ppc(:,(n_particles_found+1):tot_valid_particles) = valid_ppc_this_batch;
    else
        valid_samples(:,(n_particles_found+1):num_particles) = ...
            valid_samples_this_batch(:,1:(num_particles-n_particles_found));
        valid_weights((n_particles_found+1):num_particles) = ...
            valid_weights_this_batch(:,1:(num_particles-n_particles_found));
        err_vec((n_particles_found+1):num_particles) = ...
            valid_err_this_batch(:,1:(num_particles-n_particles_found));
        accepted_Pcc((n_particles_found+1):num_particles) = ...
            valid_Pcc_this_batch(1:(num_particles-n_particles_found));
        accepted_peaks(:,(n_particles_found+1):num_particles) = ...
            valid_peak_this_batch(:,1:(num_particles-n_particles_found));
        accepted_ppc(:,(n_particles_found+1):num_particles) = ...
            valid_ppc_this_batch(:,1:(num_particles-n_particles_found));
    end

    %update num particles found
    n_particles_found = n_particles_found + sum(samples_are_valid);

    %progress update
    fprintf('\nFound %d particles at an acceptance rate of %.3f\n', n_particles_found, n_particles_found/sim_count);

end



output.samples = valid_samples;
output.weights = valid_weights/sum(valid_weights);
output.sim_count = sim_count;
output.err_vec = err_vec;
output.accepted_Pcc = accepted_Pcc;
output.accepted_peaks = accepted_peaks;
output.accepted_ppc = accepted_ppc;
