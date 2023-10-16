

%Makes Figure 3g from the manuscript



%set font of plots
set(0,'defaulttextinterpreter','latex')

%add folders
addpath helper_funcs


%% simulation data specs

%folder formatting
data_source = '../SIM_DATA_from_ms/fig3g_example_sims/';
folder_stem = 'CC_weight_factor_';
folder_index_num = [0.1, 0.5, 0.9];
num_Pcc_vals = length(folder_index_num);

%number of replicates
num_reps = 8;

%time vector
dt = 0.5;
final_time = 50;
num_time_points = round(final_time/dt);

%percent initially infected
moi = 0.01;

%grid size
cells_wide = 50;
cells_long = 50;

%plot colours
plot_colours = {'#b3cde3', '#8c96c6', '#88419d'};

%x limit
max_time_to_plot = 40;


%% first need to process data to extract clustering statistics
extract_clustering_stats_from_data;


%% make plot

plt_handles = zeros(1,num_Pcc_vals);
figure
hold on

%loop over Pcc value
for fldr = 1:num_Pcc_vals

    %loop over replicate
    for rep = 1:num_reps

        %find end of simulation (simulation terminated when all cells infected)
        [~,length_of_time_series] = max(all_clustering_data(2:end,rep,fldr));

        %first replicate gets a handle for the legend
        num_neighbouring_cells = 6;
        if rep==1
            %plot up to end of simulation
            plt_handles(fldr) = plot(dt*(1:length_of_time_series), ...
                (1/num_neighbouring_cells)*all_clustering_data(1:length_of_time_series,rep,fldr),...
                'Color', plot_colours{fldr});
        else
            %plot up to end of simulation
            plot(dt*(1:length_of_time_series), ...
                (1/num_neighbouring_cells)*all_clustering_data(1:length_of_time_series,rep,fldr),...
                'Color', plot_colours{fldr});
        end
        %after end of simulation, kappa(t) saturates at 1, add this to plot
        plot(dt*(length_of_time_series:round(max_time_to_plot/dt)),1+0*dt*(length_of_time_series:round(max_time_to_plot/dt)),...
            'Color', plot_colours{fldr})
    end
end


%plot formatting
ylim([0, 1.2])
ylabel('$\kappa(t)$')

xlim([0,max_time_to_plot])
xlabel('Time (h)')

box on


%legend
l = legend(plt_handles, {'$P_{\mathrm{CC}} \approx 0.1$', '$P_{\mathrm{CC}} \approx 0.5$', '$P_{\mathrm{CC}} \approx 0.9$'}, ...
    'Location', 'SouthEast', 'Interpreter', 'latex');

