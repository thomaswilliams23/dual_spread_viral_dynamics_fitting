

%Makes SI Figures 3abc from the manuscript



%set font of plots
set(0,'defaulttextinterpreter','latex')

%add folders
addpath helper_funcs



%% simulation data specs

%Pcc value to plot: choose one of the below
true_Pcc = 0.1;
% true_Pcc = 0.5;
% true_Pcc = 0.9;

%Pcc and parameter values for all data
all_Pcc = [0.1, 0.5, 0.9];
Pcc_names = {'01', '05', '09'};
all_alpha = [2.28849627e-01, 1.08801559e+00, 1.87383799e+00];
all_beta =  [1.37951410e-06, 7.19841731e-07, 1.36340545e-07];

%extract information relevant to the Pcc value to be plotted
Pcc_ind = find(all_Pcc == true_Pcc);
true_alpha = all_alpha(Pcc_ind);
true_beta = all_beta(Pcc_ind);


%folder for data
target_folder = strcat('../SIM_DATA_from_ms/SIM_EST_SWEEP/Pcc_', Pcc_names{Pcc_ind},...
    '/phi_1e2/replicate_1');

%colours
colour1 = [102,194,165]/256; 
colour2 = [141,160,203]/256; 
colour3 = [231,138,195]/256; 
colour4 = [166,216,84]/256; 


%% open data from each chain
alpha_beta_1 = readmatrix(strcat(target_folder, '/all_samples_chain_1.csv'));
alpha_beta_2 = readmatrix(strcat(target_folder, '/all_samples_chain_2.csv'));
alpha_beta_3 = readmatrix(strcat(target_folder, '/all_samples_chain_3.csv'));
alpha_beta_4 = readmatrix(strcat(target_folder, '/all_samples_chain_4.csv'));

alpha_beta_1 = alpha_beta_1(2:end, 3:4);
alpha_beta_2 = alpha_beta_2(2:end, 3:4);
alpha_beta_3 = alpha_beta_3(2:end, 3:4);
alpha_beta_4 = alpha_beta_4(2:end, 3:4);


%% plot scatters

%formatting
scatter_size = 100;
num_samples_to_show = 500; %%NOTE: to avoid crashing MATLAB, kept at 500 (instead of full 2000)

figure
hold on

%first round of scatters get handles for legend
s1 = scatter(alpha_beta_1(1,1), alpha_beta_1(1,2), 100, colour1, '.');
s2 = scatter(alpha_beta_2(1,1), alpha_beta_2(1,2), 100, colour2, '.');
s3 = scatter(alpha_beta_3(1,1), alpha_beta_3(1,2), 100, colour3, '.');
s4 = scatter(alpha_beta_4(1,1), alpha_beta_4(1,2), 100, colour4, '.');

%plot all other scatters one at a time so they intersperse 
for i = 2:500
    scatter(alpha_beta_1(i,1), alpha_beta_1(i,2), 100, colour1, '.');
    scatter(alpha_beta_2(i,1), alpha_beta_2(i,2), 100, colour2, '.');
    scatter(alpha_beta_3(i,1), alpha_beta_3(i,2), 100, colour3, '.');
    scatter(alpha_beta_4(i,1), alpha_beta_4(i,2), 100, colour4, '.');
end


%% plot true value
true_val = plot(true_alpha, true_beta, 'ro', 'MarkerSize', 15);


%% r contour

%other model parameters
true_gamma = 3.366934e-01;
true_delta = 8.256588e-02;
true_log_p = 1.409457e+01;
true_c = 4.313531e-01;
true_params = [true_alpha, true_beta, true_gamma, true_delta, exp(true_log_p), true_c];

%number of latent stages
stages = [3,1];

%set up fuction for the r curve (beta as a function of alpha on contour)
beta_curve = @(x) r_contour(x, true_params, stages);

%generate many alpha values and plot r contour on these
alpha_sweep = linspace(0, 5, 1000);
plot(alpha_sweep, beta_curve(alpha_sweep), 'r--', 'LineWidth',1);



%% Pcc contour

%range of values for alpha and beta used in lookup table
alpha_vals = 2.5*(0:0.05:1);
beta_vals = 2e-6*(0:0.05:1);

%initialise
alpha_on_Pcc_contour = zeros(1,length(beta_vals));

%loop over beta values and find corresponding alpha value on Pcc contour
%%% NOTE: omit first point as points near the origin are undefined for
%%% Pcc since infection does not infect all cells
for beta_ind = 2:length(beta_vals)
   alpha_on_Pcc_contour(beta_ind) = infer_alpha_on_Pcc_contour('../SIM_DATA_from_ms/lookup_tables/Pcc_lookup', ...
       alpha_vals, beta_vals, beta_vals(beta_ind), true_Pcc); 
end

%plot
plot([0,alpha_on_Pcc_contour(2:end)], [0,beta_vals(2:end)], 'r--', 'LineWidth', 1);


%% format
xlim([0, 2.3])
ylim([0, 2e-6])

xlabel('$\alpha$')
ylabel('$\beta$')

legend([s1, s2, s3, s4, true_val], ...
    {'Chain 1', 'Chain 2', 'Chain 3', 'Chain 4', 'True Value'},...
    'FontSize', 8, 'Interpreter', 'latex');

