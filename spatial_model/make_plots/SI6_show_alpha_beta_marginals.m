

%Makes SI Figure 6



%set of font of plots
set(0, 'DefaultTextInterpreter', 'latex')


%% setup

%data structure
base_folder = '../SIM_DATA_from_ms/SIM_EST_SWEEP_with_clustering';
template_1 = '/fits_Pcc_';
folder_names = {'0.1', '0.5', '0.9'};
template_2 = '_sample_all_cells';
num_replicates = 4;


%data specs
alpha_max = 10;
beta_max = 1.5e-6;

%axis limits
alpha_plot_max = 2.5;
beta_plot_max = 2.5e7;

%true alpha beta parameters for Pcc=0.1, 0.5, 0.9, tpeak = 18
true_alpha = [1.88035441e-01, 1.10716665e+00, 8.08312311e+00];
true_beta = [9.61921398e-07, 3.91381868e-07, 5.74722863e-08];



%% loop

%loop Pcc
for fldr_ind = 1:length(folder_names)

    %load in data
    fldr_path = strcat(base_folder, template_1, folder_names{fldr_ind}, template_2);
    load(strcat(fldr_path, '/replicate_1/curr_samples'))
    load(strcat(fldr_path, '/replicate_1/curr_weights'))

    %work out weighted marginal density
    [alpha_dens, alpha_x] = ksdensity(curr_samples(1,:), 'Weights', curr_weights/sum(curr_weights),...
        'Function', 'pdf');
    [beta_dens, beta_x] = ksdensity(curr_samples(2,:), 'Weights', curr_weights/sum(curr_weights),...
        'Function', 'pdf');



    %alpha density
    figure
    hold on
    p1 = patch('XData', [0, alpha_max, alpha_max, 0], 'YData', [0, 0, 1/alpha_max, 1/alpha_max], ...
        'edgecolor', 'none', 'facecolor', [0.7, 0.7, 0.7], 'facealpha', 0.3);
    p2 = patch('XData', [alpha_x, fliplr(alpha_x)], 'YData', [0*alpha_x, fliplr(alpha_dens)], ...
        'edgecolor', 'none', 'facecolor', [0.4940 0.1840 0.5560], 'facealpha', 0.3);
    p3 = plot(true_alpha(fldr_ind)*[1,1], [0,alpha_plot_max], 'k--', 'LineWidth', 1.5);
    xlim([0, alpha_max])
    ylim([0, alpha_plot_max])
    xlabel('$\alpha$')
    ylabel('Probability density')
    if fldr_ind==3
        legend([p1, p2, p3], {'Prior', 'Posterior', 'True Value'}, 'Interpreter', 'latex', 'Location', 'northwest')
    end


    %beta density
    figure
    hold on
    p1 = patch('XData', [0, beta_max, beta_max, 0], 'YData', [0, 0, 1/beta_max, 1/beta_max], ...
        'edgecolor', 'none', 'facecolor', [0.7, 0.7, 0.7], 'facealpha', 0.3);
    p2 = patch('XData', [beta_x, fliplr(beta_x)], 'YData', [0*beta_x, fliplr(beta_dens)], ...
        'edgecolor', 'none', 'facecolor', [0.6350 0.0780 0.1840], 'facealpha', 0.3);
    p3 = plot(true_beta(fldr_ind)*[1,1], [0,beta_plot_max], 'k--', 'LineWidth', 1.5);
    xlim([0,beta_max])
    ylim([0,beta_plot_max])
    xlabel('$\beta$')
    ylabel('Probability density')
    if fldr_ind==3
        legend([p1, p2, p3], {'Prior', 'Posterior', 'True Value'}, 'Interpreter', 'latex', 'Location', 'northeast')
    end


end

