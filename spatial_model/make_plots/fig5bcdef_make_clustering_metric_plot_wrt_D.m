

%Makes Figures 5bcdef from the manuscript



%set font of plots
set(0, 'DefaultTextInterpreter', 'latex')


%% setup

%folder structure
data_source = '../SIM_DATA_from_ms/fig5bcdef_data_for_kappa';
template_1 = '/Pcc_';
template_2 = '_virus_diff_';

Pcc_vals = [0.1, 0.5, 0.9];
Pcc_names = {'0.1', '0.5', '0.9'};
diff_vals = [0.1, 1, 10, 100, Inf];
diff_names = {'01', '1', '10', '100', 'Inf'};

num_reps = 8;


%colour scheme
kappa_colours = {'#b3cde3', '#8c96c6', '#88419d'};
fluoro_plot_colours = [178,226,226; 102,194,164; 35,139,69]/256;
I_plot_colours = [254,232,200; 253,187,132; 227,74,51]/256;

%legend
legend_text = {'$P_{\mathrm{CC}}\approx0.1$', '$P_{\mathrm{CC}}\approx0.5$', '$P_{\mathrm{CC}}\approx0.9$'};


%time vector
dt = 0.1;
final_time = 40;
time_vec = dt*(0:(round(final_time/dt)-1));
num_time_points = length(time_vec);


%% loop 
%loop diffusion
for diff_ind = 1:length(diff_vals)


    figure
    hold on

    %legends
    kappa_legend = zeros(3,1);

    %loop Pcc
    for Pcc_ind = 1:length(Pcc_vals)
    
        %load data
        load(strcat(data_source, template_1, Pcc_names{Pcc_ind}, template_2, diff_names{diff_ind}, '/all_sim_data'));

        %loop over reps and plot
        for rep = 1:num_reps

            this_sim = all_sim_data{rep};
            kappa_legend(Pcc_ind) = ...
                plot(time_vec, this_sim.kappa_data, 'LineWidth', 0.5, 'Color', kappa_colours{Pcc_ind});

        end

    end


    %% work out max vertical distance

    %load data
    load(strcat(data_source, template_1, Pcc_names{1}, template_2, diff_names{diff_ind}, '/all_sim_data'));
    all_sim_data_upper = all_sim_data;

    load(strcat(data_source, template_1, Pcc_names{3}, template_2, diff_names{diff_ind}, '/all_sim_data'));
    all_sim_data_lower = all_sim_data;


    %form an average trajectory
    mean_upper_traj = zeros(1,num_time_points);
    mean_lower_traj = zeros(1,num_time_points);

    for rep = 1:num_reps
        mean_upper_traj = mean_upper_traj + all_sim_data_upper{rep}.kappa_data/num_reps;
        mean_lower_traj = mean_lower_traj + all_sim_data_lower{rep}.kappa_data/num_reps;
    end


    %now find vertical distance
    vert_distance_over_time = abs(mean_upper_traj-mean_lower_traj);

    %and report max
    [~,max_dist_ind] = max(vert_distance_over_time);


    %now plot
    plot(max_dist_ind*dt*[1,1], [mean_lower_traj(max_dist_ind), mean_upper_traj(max_dist_ind)], ...
        'LineStyle', '--', 'LineWidth', 1.3, 'Color', [0.8500 0.3250 0.0980])

    %label max vertical distance for last plot
    if diff_ind ==5
        text(20, 0.6, 'Max. Vert. Distance', ...
            'Interpreter','latex','Color', [0.8500 0.3250 0.0980], 'FontSize',8)
    end


    %% formatting

    %plot legend for last plot
    if diff_ind ==5
       l1 = legend(kappa_legend, legend_text, 'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 9);
    end


    %formatting
    ylim([0, 1.2])
    xlabel('Time (h)')
    ylabel('$\kappa(t)$')
    

end



