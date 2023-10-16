function params_init_struct = make_init_params(true_params_struct, params_to_fit, prior_dists, fixed_params_are_zero)

%unpack priors
alpha_prior_dist = prior_dists{1};
beta_prior_dist = prior_dists{2};
gamma_prior_dist = prior_dists{3};
delta_prior_dist = prior_dists{4};
log_p_prior_dist = prior_dists{5};
c_prior_dist = prior_dists{6};

%random init params
alpha_init = random(alpha_prior_dist,1);
beta_init = random(beta_prior_dist,1);
gamma_init = random(gamma_prior_dist,1);
delta_init = random(delta_prior_dist, 1);
log_p_init = random(log_p_prior_dist, 1);
c_init = random(c_prior_dist,1);

%if not fitting, set param
for param = 1:length(true_params_struct)
    if ~(ismember(1,params_to_fit))
        if ismember(1,fixed_params_are_zero)
            alpha_init = 0;
        else
            alpha_init = true_params_struct.alpha_param;
        end
    end
    if ~(ismember(2,params_to_fit))
        if ismember(2,fixed_params_are_zero)
            beta_init = 0;
        else
            beta_init = true_params_struct.beta_param;
        end
    end
    if ~(ismember(3,params_to_fit))
        if ismember(3,fixed_params_are_zero)
            gamma_init = 0;
        else
            gamma_init = true_params_struct.gamma_param;
        end
    end
    if ~(ismember(4,params_to_fit))
        if ismember(4,fixed_params_are_zero)
            delta_init = 0;
        else
            delta_init = true_params_struct.delta_param;
        end
    end
    if ~(ismember(5,params_to_fit))
        if ismember(5,fixed_params_are_zero)
            log_p_init = -Inf;
        else
            log_p_init = log(true_params_struct.p);
        end
    end
    if ~(ismember(6,params_to_fit))
        if ismember(6,fixed_params_are_zero)
            c_init = 0;
        else
            c_init = true_params_struct.c;
        end
    end
end


%initialise
params_init_struct = true_params_struct;

params_init_struct.alpha_param = alpha_init;
params_init_struct.beta_param = beta_init;
params_init_struct.gamma_param = gamma_init;
params_init_struct.delta_param = delta_init;
params_init_struct.p = exp(log_p_init);
params_init_struct.c = c_init;


end