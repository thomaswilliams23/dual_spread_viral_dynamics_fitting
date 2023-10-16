


%% RUNS AN ARRAY OF SIMULATION ESTIMATIONS FOR THE SPATIAL MODEL


%  iterates over Pcc = 0.1, 0.5, 0.9 and sample sizes for the clustering 
%  metric of 50, 100, 500, 1000 and all cells



rng('shuffle')

parpool(8)

addpath sim_est_helpers
addpath mc_helpers


base_folder = 'SAMPLE_SIZE_SWEEP_RESULTS';
mkdir(base_folder);

save_samples_at_all_generations = 1;


%% outer parameters setup

%Pcc sweep values
Pcc_sweep = [0.1, 0.5, 0.9];
Pcc_names = {'0.1', '0.5', '0.9'};

%previously computed alpha and beta target values
alpha_sweep_vals = [1.88035441e-01, 1.10716665e+00, 8.08312311e+00];
beta_sweep_vals = [9.61921398e-07, 3.91381868e-07, 5.74722863e-08];

%sample size values (-1 => all cells)
sample_size_sweep = [-1, 1000, 500, 100, 50];
sample_size_names = {'all', '1000', '500', '100', '50'};

%number of iterations per parameter combination
num_iters = 4;



%loop over Pcc values
for Pcc_ind = 1:length(Pcc_sweep)

    %loop over sample size
    for sample_size_ind = 1:length(sample_size_sweep)

        %create parent output folder
        array_point_folder = strcat(base_folder, '/fits_Pcc_', Pcc_names{Pcc_ind}, ...
            '_sample_', sample_size_names{sample_size_ind}, '_cells');
        mkdir(array_point_folder);
        

        %loop over iterations
        for iter = 1:num_iters


            %construct output folder
            output_folder = strcat(array_point_folder, '/replicate_', num2str(iter));
            strcat(output_folder);


            %% (meta) parameters
            params_to_fit =  1:2; %corresponds to alpha and beta [THIS IS FIXED]
            fixed_params_are_zero = []; %these ones are fixed as zero when fitting [THIS IS FIXED]

            %define number of accepted particles per generation
            num_particles = 500;

            %number of samples to take to calculate cluster metric (-1 => all cells)
            num_samples_for_clustering = sample_size_sweep(sample_size_ind);
            

            %thresholds
            num_thresh = 5; %number of rounds of importance sampling
            fluor_err_thresh_quantile = 0.5;
            clust_err_thresh_quantile = 0.5;

            %number of species [THIS IS FIXED]
            num_species = 1;


            %set of time values to export information
            sim_times = 0:0.1:30;

            %observation times for refitting
            obs_times = 3*(1:10);  %<--- NEEDS TO BE A SUBSET OF SIM_TIMES

            %percent of cell sheet initially infected [THIS IS FIXED]
            moi = 0.01;

            %basic biological parameters
            biology_parameters;  %<---  SET TISSUE SIZE HERE
            sheet_size = cells_wide*cells_long;


            %parameter names
            param_names = {'\alpha', '\beta', '\gamma', '\delta', 'log p', 'c'};

            %number of latent stages [THIS IS FIXED]
            latent_stages = 3;


            %% find indices of obs times in larger output times vector

            indices_of_obs_times = zeros(1,length(obs_times));

            EPS_SMALL = 1e-6;
            for obs_time = obs_times
                i=1;
                while ~(abs(sim_times(i) - obs_time)<EPS_SMALL)
                    i=i+1;
                end
                indices_of_obs_times((obs_time == obs_times)) = i;
            end


            %% set true model parameters

            %put in params for Pcc
            alpha_true = alpha_sweep_vals(Pcc_ind);
            beta_true = beta_sweep_vals(Pcc_ind);

            %other model parameters (fixed)
            gamma_true = 3.366934e-01;
            delta_true = 8.256588e-02;
            log_p_true = 1.409457e+01;
            c_true = 4.313531e-01;

            param_means = [alpha_true, beta_true, gamma_true, delta_true, log_p_true, c_true];
            true_params = param_means;

            %convert true_params to a struct
            set_true_params;

            %need to set number of species and instruct the sim to save grid data
            true_params_struct.num_species = num_species;
            true_params_struct.save_clustering_data = 1;
            true_params_struct.clustering_calc_times = obs_times;
            true_params_struct.num_samples_for_clustering = num_samples_for_clustering;


            %% initial model run --- generate observational data

            %initialise observational data structures
            obs_clustering_data = 0*obs_times;
            obs_cumul_inf_data = 0*obs_times;

            Pcc_estimates = 0*obs_times;
            peak_estimates = 0*obs_times;

            errs_encountered = 0*obs_times;

            %loop over observation times points and run one simulation per point
            parfor time_point = 1:length(obs_times)

                %run the model
                [true_output_point] = single_run_as_func(sim_times, true_params_struct);

                %extract one time point for the fluorescent cell and
                %clustering measures
                obs_clustering_data(time_point) = true_output_point.purity_data(time_point);
                obs_cumul_inf_data(time_point) = (1-moi)*true_output_point.prop_infected(indices_of_obs_times(time_point));


                %record the actual Pcc and peak time obtained from this simulation
                Pcc_estimates(time_point) = true_output_point.net_CC(end)/(true_output_point.net_CC(end) + true_output_point.net_CF(end));

                [~,max_ind] = max(true_output_point.net_I);
                peak_estimates(time_point) = sim_times(max_ind);

                %any errors in computing clustering metrics
                errs_encountered(time_point) = true_output_point.tot_purity_errs;
            end

            %displays if errors encountered in computing the clustering metric
            if sum(errs_encountered)
                fprintf('\nCAUTION: encountered %d void clustering calculations in constructing observed data\n', sum(errs_encountered));
            end


            %% set the prior functions

            make_prior_funcs;


            %% run PMC-MCMC - step 1 (rejection sampling)

            disp('RUNNING REJECTION SAMPLING')

            %accepts this proportion of pure prior samples
            prior_accept_rate = 0.3;


            %draw some initial parameters
            params_init_struct = make_init_params(true_params_struct, ...
                params_to_fit, prior_dists, fixed_params_are_zero);


            %initialise vectors to capture posterior data for each epsilson and number
            %of simulations called in each case
            num_sims = zeros(1,num_thresh);


            %do an initial run using rejection sampling to get first set of samples
            init_run = ABC_rejection_sampler_with_clustering_stats(params_init_struct, ...
                params_to_fit, moi, sheet_size, sim_times, indices_of_obs_times, ...
                obs_cumul_inf_data, obs_clustering_data, prior_accept_rate, ...
                num_particles, prior_dists);
            curr_samples = init_run.samples;
            tot_sims = init_run.sim_count;
            fluor_err_vec = init_run.fluor_err_vec;
            clust_err_vec = init_run.clust_err_vec;

            %initialise weights (all ones)
            curr_weights = (1/num_particles)* ones(1,num_particles);

            %update target error
            eps_fluor_err_target = quantile(fluor_err_vec, fluor_err_thresh_quantile);
            eps_clust_err_target = quantile(clust_err_vec, clust_err_thresh_quantile);


            %keep track of all the samples at each generation (if requested)
            if save_samples_at_all_generations
                initial_samples = curr_samples;
            end


            %print errors if any
            if (init_run.num_void_calc)
                fprintf('\nCAUTION: encountered %d void clustering calculations during sampling\n', init_run.num_void_calc);
            end

            disp('COMPLETED REJECTION SAMPLING')


            %% run PMC-MCMC - step 2

            % initialise a structure to save out all samples from each 
            % generation of sampling (if option is on)
            if save_samples_at_all_generations
                all_samples = zeros(num_thresh,num_particles,length(params_to_fit)+1);
            end

            %loop over generations
            for gen_num = 1:num_thresh

                fprintf('\nRUNNING ROUND %d OF IMPORTANCE SAMPLING\n', gen_num);


                %call ABC importance sampler function, update samples, weights, total
                %sim count
                this_run = ABC_importance_sampler_with_clustering_stats(params_init_struct, ...
                    params_to_fit, moi, sheet_size, sim_times, indices_of_obs_times, ...
                    obs_cumul_inf_data, obs_clustering_data, eps_fluor_err_target, ...
                    eps_clust_err_target, num_particles,prior_funcs, curr_samples, curr_weights);

                %save out data
                curr_samples = this_run.samples;
                curr_weights = this_run.weights;

                %tighten error bounds
                fluor_err_vec = this_run.fluor_err_vec;
                eps_fluor_err_target = quantile(fluor_err_vec,fluor_err_thresh_quantile);

                clust_err_vec = this_run.clust_err_vec;
                eps_clust_err_target = quantile(clust_err_vec,clust_err_thresh_quantile);


                %save this generation of samples if requested
                if save_samples_at_all_generations
                    all_samples(gen_num, :, 1) = curr_weights';
                    all_samples(gen_num, :, 2:end) = curr_samples';
                end


                %only save posterior predictive check (PPC) data at the last simulation
                if (gen_num == num_thresh)
                    final_ppc = this_run.accepted_ppc;
                    final_Pcc_samples = this_run.accepted_Pcc;
                    final_peak_samples = this_run.accepted_peaks;
                end


                %print errors if any
                if (this_run.num_void_calc)
                    fprintf('\nCAUTION: encountered %d void clustering calculations during sampling\n', this_run.num_void_calc);
                end

            end

            disp('COMPLETED IMPORTANCE SAMPLING')





            %% save stuff
            mkdir(output_folder)
            save(strcat(output_folder, '/curr_samples'), 'curr_samples');
            save(strcat(output_folder, '/curr_weights'), 'curr_weights');
            save(strcat(output_folder, '/obs_cumul_inf_data'), 'obs_cumul_inf_data');
            save(strcat(output_folder, '/Pcc_estimates'), 'Pcc_estimates');
            save(strcat(output_folder, '/peak_estimates'), 'peak_estimates');
            save(strcat(output_folder, '/final_Pcc_samples'), 'final_Pcc_samples');
            save(strcat(output_folder, '/final_peak_samples'), 'final_peak_samples');
            save(strcat(output_folder, '/final_ppc'), 'final_ppc');
            save(strcat(output_folder, '/err_vec'), 'fluor_err_vec');
            save(strcat(output_folder, '/pur_dist_vec'), 'clust_err_vec');
            if save_samples_at_all_generations
                save(strcat(output_folder, '/all_samples'), 'all_samples');
            end



        end
    end
end