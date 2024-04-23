

%Makes SI Figure 5



%set of font of plots
set(0, 'DefaultTextInterpreter', 'latex')


%% setup

%data structure
base_folder = '../SIM_DATA_from_ms/SIM_EST_SWEEP';
template_1 = '/Pcc_';
folder_names = {'01', '05', '09'};
template_2 = '/phi_1e2/replicate_1';

%data specs
alpha_max = 2.5;
beta_max = 2e-6;

%axis limits
alpha_plot_max = 1.5;
beta_plot_max = 14e5;

%true alpha beta parameters for Pcc=0.1, 0.5, 0.9, r=0.52
true_alpha = [2.28849627e-01, 1.08801559e+00, 1.87383799e+00];
true_beta = [1.37951410e-06, 7.19841731e-07, 1.36340545e-07];



%% loop

%loop Pcc
for fldr_ind = 1:length(folder_names)

    %load in data
    target_folder = strcat(base_folder, template_1, folder_names{fldr_ind}, template_2);

    alpha_beta_1 = readmatrix(strcat(target_folder, '/all_samples_chain_1.csv'));
    alpha_beta_2 = readmatrix(strcat(target_folder, '/all_samples_chain_2.csv'));
    alpha_beta_3 = readmatrix(strcat(target_folder, '/all_samples_chain_3.csv'));
    alpha_beta_4 = readmatrix(strcat(target_folder, '/all_samples_chain_4.csv'));

    all_alpha_samples = [alpha_beta_1(2:end,3); alpha_beta_2(2:end,3); ...
                         alpha_beta_3(2:end,3); alpha_beta_4(2:end,3)];
    all_beta_samples = [alpha_beta_1(2:end,4); alpha_beta_2(2:end,4); ...
                        alpha_beta_3(2:end,4); alpha_beta_4(2:end,4)];


    %work out weighted marginal density
    [alpha_dens, alpha_x] = ksdensity(all_alpha_samples, 'Function', 'pdf');
    [beta_dens, beta_x] = ksdensity(all_beta_samples, 'Function', 'pdf');



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
    xlim([0, beta_max])
    ylim([0, beta_plot_max])
    xlabel('$\beta$')
    ylabel('Probability density')
    if fldr_ind==3
        legend([p1, p2, p3], {'Prior', 'Posterior', 'True Value'}, 'Interpreter', 'latex', 'Location', 'northeast')
    end

end


