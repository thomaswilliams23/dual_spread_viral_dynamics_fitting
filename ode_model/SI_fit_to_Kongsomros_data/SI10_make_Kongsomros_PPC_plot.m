

%Makes SI Figure 10




%set font of plots
set(0, 'DefaultTextInterpreter', 'latex')

%% setup

addpath ../matlab_helpers

%time vector
time_vec = 0:0.1:25;

%load in data
posterior_samples = readmatrix("posterior_samples.csv");
posterior_samples = posterior_samples(2:end,2:end);

%find 2.5 and 97.5 percentiles
prctile_025 = prctile(posterior_samples, 2.5, 1);
prctile_975 = prctile(posterior_samples, 97.5, 1);


%shade 95% CI region
p1 = patch('XData', [time_vec, fliplr(time_vec)], 'YData', [prctile_975, fliplr(prctile_025)],...
    'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.5, 'EdgeColor', 'none');



%% show the sample used for our fits

%parameters
moi = 0.5;
final_time = 25;
t_span = [0, final_time];


%stages
latent_stages = 3;
inf_stages = 1;

%default params
alpha_def = 0.95027068764583;
beta_def =  1.37914512086743e-06;
gamma_def = 3.366934e-01;
delta_def = 8.256588e-02;
p_def = exp(1.409457e+01);
c_def = 4.313531e-01;


%ICs
T_0 = 1-moi;
E_CF_0 = 0;
E_CC_0 = 0;
I_CF_0 = 0;
I_CC_0 = 0;
I_r_dead_0 = 0;
I_d_0 = moi;
I_d_dead_0 = 0;
V_0 = 0;


%indices setup
num_compartments = 2*latent_stages + 3*inf_stages + 6;

target_ind  = 1;
latent_CF_start = target_ind + 1;
latent_CC_start = latent_CF_start + latent_stages;
inf_CF_start = latent_CC_start + latent_stages;
inf_CC_start = inf_CF_start + inf_stages;
inf_r_dead_ind = inf_CC_start + inf_stages;
inf_d_start = inf_r_dead_ind+1;
inf_d_dead_ind = inf_d_start + inf_stages;
v_ind = inf_d_dead_ind + 1;
CF_cum_ind = v_ind + 1;
CC_cum_ind = num_compartments;


%ICs
init_conditions = zeros(num_compartments, 1);

init_conditions(target_ind) = T_0;
init_conditions(latent_CF_start) = E_CF_0;
init_conditions(latent_CC_start) = E_CC_0;
init_conditions(inf_CF_start) = I_CF_0;
init_conditions(inf_CC_start) = I_CC_0;
init_conditions(inf_r_dead_ind) = I_r_dead_0;
init_conditions(inf_d_start)= I_d_0;
init_conditions(inf_d_dead_ind) = I_d_dead_0;
init_conditions(v_ind) = V_0;

params = [alpha_def, beta_def, gamma_def, delta_def, p_def, c_def, ...
    latent_stages, inf_stages];


%solve
anon_ode = @(t, Y) dual_spread_ODE_model(t, Y, params);
[t_out, tiv_out] = ode45(anon_ode, t_span, init_conditions);


%add to plot
all_infected_recipient = sum(tiv_out(:,inf_CF_start:(inf_CF_start+inf_stages-1)),2)...
                        +sum(tiv_out(:,inf_CC_start:(inf_CC_start+inf_stages-1)),2)...
                        +tiv_out(:,inf_r_dead_ind);    


hold on
p2 = plot(t_out, 100*all_infected_recipient/init_conditions(1), 'LineWidth', 1.5, ...
    'Color', 'k', 'LineStyle','--');



%% plot data from Kongsomros et al.
data = 100*[0.022, 0.12,  0.933, 0.99,  1];       %H1N1 coculture data
t_data = [1, 3, 6, 12, 24];                       %time points from data

p3 = plot(t_data, data, 'rx', 'MarkerSize', 10);



%% formatting

ylim([0, 110])

xlabel('Time (h)')
ylabel('Percent fluorescent recipient cells')

legend([p1, p2, p3], {'95\% CI', 'Default parameter set', 'Kongsomros \textit{et al.} data'}, ...
    'Location', 'southeast', 'Interpreter','latex', 'FontSize', 9)

