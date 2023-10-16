function [Pcc_samples, peak_samples] = compute_Pcc_and_peak_prior(params_to_fit, ...
    fixed_params_are_zero, samples_per_param, final_time, moi, true_params)


%Draw random samples of the parameters to fit and compute Pcc in each case.
%Returns a vector of Pcc values from which a histogram could be
%constructed.

%ALSO ASSUMES PRIORS AREN'T GOING TO CHANGE

prev_path = cd;

cd ../
addpath sim_est_helpers



disp("CAUTION (in compute_Pcc_and_peak_prior): CHECK PRIORS IN SOURCE CODE")
disp("Also - no support yet for fixed zero parameters")

obs_times = 0:0.1:final_time;
latent_stages = 3;

set_true_params;
true_params_struct = true_params_struct;

Pcc_samples = zeros(1,samples_per_param^length(params_to_fit));
peak_samples = zeros(1,samples_per_param^length(params_to_fit));


parfor sample = 1:samples_per_param^length(params_to_fit)

    params_temp = true_params_struct;


    for prm=1:length(params_to_fit)
        if params_to_fit(prm)==1
            params_temp.alpha_param = 5*rand;
        end
        if params_to_fit(prm)==2
            params_temp.beta_param = 5e-6*rand;
        end
        if params_to_fit(prm)==3
            params_temp.gamma_param = 1*rand;
        end
        if params_to_fit(prm)==4
            params_temp.delta_param = exp(normrnd(log(true_params(4)),0.1*log(true_params(4))));
        end
        if params_to_fit(prm)==5
            params_temp.p = exp(20*rand);
        end
        if params_to_fit(prm)==6
            params_temp.c = exp(normrnd(log(true_params(6)),0.1*log(true_params(6))));
        end
    end



    [model_sim] = single_run_as_func(obs_times, params_temp);
    
    Pcc_samples(sample) = model_sim.net_CC(end)/(model_sim.net_CC(end) + model_sim.net_CF(end));

    [~,peak_ind] = max(model_sim.net_I);
    peak_samples(sample) = obs_times(peak_ind);


end


cd(prev_path);