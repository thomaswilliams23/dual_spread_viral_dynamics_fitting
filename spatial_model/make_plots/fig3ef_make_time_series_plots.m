

% Makes Figures 3def from the manuscript



%set font of plots
set(0,'defaulttextinterpreter','latex')


%% simulation data specs

%folder formatting
folder_stem = '../SIM_DATA_from_ms/fig3ef_example_sims/CC_weight_factor_';
folder_names = [0.1, 0.5, 0.9];
legend_names = {'$P_{\mathrm{CC}} \approx 0.1$', '$P_{\mathrm{CC}} \approx 0.5$', '$P_{\mathrm{CC}} \approx 0.9$'};

%number of replicates
num_reps = 8;

%vector of time values
dt = 0.1;
num_time_points = 500;
time_vec = dt*(1:num_time_points);


%colours
fluoro_plot_colours = [178,226,226; 102,194,164; 35,139,69]/256;
I_plot_colours = [254,232,200; 253,187,132; 227,74,51]/256;



%% net I
net_I_handles = zeros(1,3);

figure
%loop over Pcc values
for i = 1:length(folder_names)

    %load data
    load(strcat(folder_stem, num2str(folder_names(i)), '/net_I_data'));

    %plot (first one gets a handle for the legend)
    net_I_handles(i) = plot(time_vec, net_I_data(:,1), 'LineWidth', 0.5, 'Color', I_plot_colours(i,:));

    hold on
    plot(time_vec, net_I_data(:,2:end), 'LineWidth', 0.5, 'Color', I_plot_colours(i,:));
end

%plot formatting
peak_time = 17;
ymax = 0.8;
plot(peak_time*[1,1], [0,ymax], 'k:', 'LineWidth', 1.5)
ylim([0,ymax])
xlim([0,40])
text(18.5, 0.72, '$t_{\mathrm{peak}}$', 'FontSize', 10)

%add legend
l1 = legend(net_I_handles, legend_names, 'Location', 'SouthEast', 'Interpreter', 'latex');

%add labels
xlabel('Time (h)')
ylabel('Prop. Infected')



%% fluorescent cells
fluoro_handles = zeros(1,3);

figure
%loop over Pcc values
for i = 1:length(folder_names)

    %load data
    load(strcat(folder_stem, num2str(folder_names(i)), '/prop_infected_data'));

    %plot (first one gets a handle in the legend)
    fluoro_handles(i) = plot(time_vec, prop_infected_data(:,1), 'LineWidth', 0.5, 'Color', fluoro_plot_colours(i,:));

    hold on
    plot(time_vec, prop_infected_data(:,2:end), 'LineWidth', 0.5, 'Color', fluoro_plot_colours(i,:));
end

%add legend
l2 = legend(fluoro_handles, legend_names, 'Location', 'SouthEast', 'Interpreter', 'latex');


%axis labels and limits
xlabel('Time (h)')
ylabel('Prop. Fluorescent Susceptible')

ylim([0, 1.1])
xlim([0,40])




