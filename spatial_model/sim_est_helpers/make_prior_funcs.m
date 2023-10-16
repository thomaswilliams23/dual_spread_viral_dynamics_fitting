% define priors 
% [NOTE: this code is desgined so that it can easily be extended to
% estimating all parameters, here, only the alpha and beta priors are
% actually ever used]


max_alpha = 10;
alpha_prior_dist = makedist('Uniform', 'lower',0, 'upper',max_alpha);
alpha_prior_pdf = @(x) unifpdf(x,0,max_alpha);

beta_prior_dist = makedist('Uniform', 'lower', 0, 'upper', 5e-6);
beta_prior_pdf = @(x) unifpdf(x,0,5e-6);


%UNUSED IN THIS WORK
max_gamma = 1;
gamma_prior_dist = makedist('Uniform', 'lower', 0, 'upper', max_gamma);
gamma_prior_pdf = @(x) unifpdf(x,0,max_gamma);

delta_prior_dist = makedist('Lognormal', 'mu', log(delta_true), 'sigma', abs(log(delta_true))*0.1);
delta_prior_pdf = @(x) normpdf(log(x),log(delta_true),abs(log(delta_true))*0.1);

max_log_p = 20;
log_p_prior_dist = makedist('Uniform', 'lower', 0, 'upper',max_log_p);
log_p_prior_pdf = @(x) unifpdf(x,0,max_log_p);

c_prior_dist = makedist('Lognormal', 'mu', log(c_true), 'sigma', abs(log(c_true))*0.1);
c_prior_pdf = @(x) normpdf(log(x),log(c_true),abs(log(c_true))*0.1);


%gather priors together
prior_funcs = {alpha_prior_pdf, beta_prior_pdf, gamma_prior_pdf, delta_prior_pdf, ...
    log_p_prior_pdf, c_prior_pdf};      

prior_dists = {alpha_prior_dist, beta_prior_dist, gamma_prior_dist, delta_prior_dist, ...
    log_p_prior_dist, c_prior_dist};
