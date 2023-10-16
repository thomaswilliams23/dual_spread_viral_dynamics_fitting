%define priors (and make some parameters fixed)
max_alpha = 10;%2;
alpha_prior_dist = makedist('Uniform', 'lower',0, 'upper',max_alpha);
alpha_prior_pdf = @(x) unifpdf(x,0,max_alpha);

beta_prior_dist = makedist('Uniform', 'lower', 0, 'upper', 5e-6);
beta_prior_pdf = @(x) unifpdf(x,0,5e-6);

max_gamma = 1;
gamma_prior_dist = makedist('Uniform', 'lower', 0, 'upper', max_gamma);
gamma_prior_pdf = @(x) unifpdf(x,0,max_gamma);

delta_prior_dist = delta_pdf;
delta_prior_pdf = @(x) normpdf(log(x),log(delta_mean),abs(log(delta_mean))*sig_prop);

max_log_p = 20;
log_p_prior_dist = makedist('Uniform', 'lower', 0, 'upper',max_log_p);
log_p_prior_pdf = @(x) unifpdf(x,0,max_log_p);

c_prior_dist = c_pdf;
c_prior_pdf = @(x) normpdf(log(x),log(c_mean),abs(log(c_mean))*sig_prop);


%relevant priors
prior_funcs = {alpha_prior_pdf, beta_prior_pdf, gamma_prior_pdf, delta_prior_pdf, ...
    log_p_prior_pdf, c_prior_pdf};      

prior_dists = {alpha_prior_dist, beta_prior_dist, gamma_prior_dist, delta_prior_dist, ...
    log_p_prior_dist, c_prior_dist};