

%Makes SI Figures 7abc



%set font of plots
set(0, 'DefaultTextInterpreter', 'latex')


%% setup

%data structure
data_source = '../SIM_DATA_from_ms/fig5bcdef_data_for_kappa';
template_1 = '/Pcc_';
template_2 = '_virus_diff_';

Pcc_vals = [0.1, 0.5, 0.9];
Pcc_names = {'0.1', '0.5', '0.9'};
diff_vals = [0.1, 1, 10, 100, Inf];
diff_names = {'01', '1', '10', '100', 'Inf'};

diff_legend_names = {'$10^{-1}$', '$10^{0}$', '$10^{1}$', '$10^{2}$', '$\infty$'};
diff_colours = {'#fa9fb5', '#f768a1', '#c51b8a', '#7a0177', '#000000'}; 

num_reps = 8;


%time vector
dt = 0.1;
final_time = 40;
time_vec = dt*(0:(round(final_time/dt)-1));


%% loop
%loop Pcc
for Pcc_ind = 1:length(Pcc_vals)

    %legends
    fig_legend = zeros(length(diff_vals),1);

    %set up figure
    figure
    hold on

    %loop iterations
    for rep = 1:num_reps
    
        %loop over reps
        for diff_ind = length(diff_vals):-1:1
                
            %load data
            load(strcat(data_source, template_1, Pcc_names{Pcc_ind}, ...
                template_2, diff_names{diff_ind}, '/all_sim_data'));
            this_sim = all_sim_data{rep};

            %add to plot
            fig_legend(diff_ind) = ...
                plot(time_vec, this_sim.kappa_data, 'LineWidth', 0.5, 'Color', diff_colours{diff_ind});

        end
    end


    %fix formatting
    if Pcc_ind ==3
        l1 = legend(fig_legend, diff_legend_names, 'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 9);
        ylim([0, 1.2])
        xlabel('Time (h)')
        ylabel('$\kappa(t)$')
    else
        ylim([0, 1.2])
        xlabel('Time (h)')
        ylabel('$\kappa(t)$')
    end

end


