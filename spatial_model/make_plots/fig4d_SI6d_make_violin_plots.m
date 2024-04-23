

%Makes Figure 4d from the manuscript (and also SI Figure 2d)



%set font of plots
set(0,'defaulttextinterpreter','latex')

%add folders
addpath helper_funcs




%% simulation data specs

%whether or not clustering data was used for fit (if 0, generates SI Figs 2d)
cluster_data = 1;

%folder for data
if cluster_data
    folder_template_1 = "../SIM_DATA_from_ms/SIM_EST_SWEEP_with_clustering/fits_Pcc_";
    folder_template_2 = "_sample_all_cells";
else
    folder_template_1 = "../SIM_DATA_from_ms/SIM_EST_SWEEP_fluoro_only/fits_Pcc_";
    folder_template_2 = "_cumul_I_only";
end

%folder for prior samples
prior_directory = 'helper_funcs/Pcc_and_peak_prior_densities';


%Pcc and parameter values for all data
folder_indices = {'0.1', '0.5', '0.9'};
true_Pcc_vals = [0.1, 0.5, 0.9];
all_alpha = [1.88035441e-01, 1.10716665e+00, 8.08312311e+00];
all_beta =  [9.61921398e-07, 3.91381868e-07, 5.74722863e-08];

%labels
section_labels = {{'True $P_{\mathrm{CC}}\approx$ 0.1, ','($\alpha=0.188$, $\beta=9.62\times 10^{-7}$)'},...
    {'True $P_{\mathrm{CC}}\approx$ 0.5, ','($\alpha=1.11$, $\beta=3.91\times 10^{-7}$)'},...
    {'True $P_{\mathrm{CC}}\approx$ 0.9, ','($\alpha=8.08$, $\beta=5.75\times 10^{-8}$)'}};


%formatting
plot_font_size = 9;
fig_width = 20;
fig_ratio = 4;
Pcc_fig_height = 1.5;


%data specs
num_samples_per_sim = 500;
num_iters_to_plot = 4;
tot_iters = 4;

%number of objects to plot
tot_num_obj = length(folder_indices)*(num_iters_to_plot+1)+1;

%highlights a particular section of the plot (in this case bolds the label)
highlight_first_replicate = 1;


%colours
mean_colour = 'k';
med_colour = 'k:';

prior_plot_colour = [0, 0, 0];
prior_interior_plot_colour = [0.4, 0.4, 0.4];

Pcc_plot_colour = [0, 0.4470, 0.7410];
peak_plot_colour = [49/256, 153/256, 2/256];

Pcc_interior_plot_colour = Pcc_plot_colour*1.2;
Pcc_interior_plot_colour(Pcc_interior_plot_colour>1)=1;

peak_interior_plot_colour = peak_plot_colour*1.2;
peak_interior_plot_colour(peak_interior_plot_colour>1)=1;

%labels
axis_labels = {'Prior'};
for i = 1:length(folder_indices)
    axis_labels{(i-1)*(num_iters_to_plot+1)+2} = 'Post. Wt. Means';

    if highlight_first_replicate
        axis_labels{(i-1)*(num_iters_to_plot+1)+3} ='$\textbf{Replicate 1}$';
    else
        axis_labels{(i-1)*(num_iters_to_plot+1)+3} = 'Replicate 1';
    end

    for j = 2:num_iters_to_plot
        axis_labels{(i-1)*(num_iters_to_plot+1)+2+j} = strcat(['Replicate ', num2str(j)]);
    end
end


%% load priors
load(strcat(prior_directory, '/peak_prior_samples.mat'));
load(strcat(prior_directory, '/Pcc_prior_samples.mat'));


% plot prior violin for Pcc
figure(1)
violin(Pcc_prior_samples', 'x', [0,1], 'facecolor', prior_interior_plot_colour, ...
    'edgecolor', prior_plot_colour, 'mc', mean_colour, 'medc', med_colour);
legend off
hold on

% plot prior violin for tpeak
figure(2)
violin(peak_prior_samples', 'x', [0,1], 'facecolor', prior_interior_plot_colour, ...
    'edgecolor', prior_plot_colour, 'mc', mean_colour, 'medc', med_colour);
legend off
hold on



%% loop over Pcc values

%initialise
actual_Pcc = zeros(1,length(folder_indices));
actual_peak = zeros(1,length(folder_indices));

%loop
for folder_ind = 1:length(folder_indices)

    %folder name
    source_folder = strcat(folder_template_1, folder_indices{folder_ind}, folder_template_2);


    %% initialise
    all_Pcc_estimates = zeros(num_samples_per_sim, num_iters_to_plot);
    all_peak_estimates = zeros(num_samples_per_sim, num_iters_to_plot);
    all_weights = zeros(num_samples_per_sim, num_iters_to_plot);

    actual_Pcc_estimates = zeros(1,tot_iters);
    actual_peak_estimates = zeros(1,tot_iters);


    %% calculate mean actual Pcc and peak from the simulations
    for iter = 1:tot_iters
        load(strcat(source_folder, '/replicate_', num2str(iter), '/Pcc_estimates'));
        load(strcat(source_folder, '/replicate_', num2str(iter), '/peak_estimates'));

        actual_Pcc_estimates(iter) = mean(Pcc_estimates);
        actual_peak_estimates(iter) = mean(peak_estimates);
    end

    actual_Pcc(folder_ind) = mean(actual_Pcc_estimates);
    actual_peak(folder_ind) = mean(actual_peak_estimates);


    %% loop replicates
    for iter = 1:num_iters_to_plot

        %load data
        load(strcat(source_folder, '/replicate_', num2str(iter), '/final_Pcc_samples'));
        load(strcat(source_folder, '/replicate_', num2str(iter), '/final_peak_samples'));
        load(strcat(source_folder, '/replicate_', num2str(iter), '/curr_weights'));

        %save into big array
        all_Pcc_estimates(:,iter)=final_Pcc_samples;
        all_peak_estimates(:,iter)=final_peak_samples;
        all_weights(:,iter)=curr_weights;
    end



    %% load weighted means
    load(strcat(source_folder, '/Pcc_WMs.mat'));
    load(strcat(source_folder, '/peak_WMs.mat'));


    %% plot posterior violins and weighted mean boxplots

    %positions of plots
    vlns_pos = (folder_ind-1)*(num_iters_to_plot+1) + 1 + (1:num_iters_to_plot);
    box_pos = (folder_ind-1)*(num_iters_to_plot+1) + 1;


    %Pcc figure
    figure(1)

    %plot weighted violins of all Pcc samples
    weighted_violin(all_Pcc_estimates, 'weights', all_weights,  'x', vlns_pos, 'facecolor', Pcc_interior_plot_colour,...
        'edgecolor', Pcc_plot_colour, 'mc', mean_colour, 'medc', []);

    %plot boxplot of weighted means for this Pcc
    boxplot(Pcc_WMs, 'positions', box_pos, 'colors', Pcc_plot_colour, 'widths', 0.7);

    %adjust colours
    h=findobj(gca,'tag','Outliers');
    set(h, 'MarkerEdgeColor', Pcc_plot_colour)

    %formatting
    hold on
    legend off

    %plot divider
    if folder_ind>1
        plot((vlns_pos(1)-1.5)*[1,1], [-0.3, Pcc_fig_height], 'k:');
    else
        plot((vlns_pos(1)-1.5)*[1,1], [-0.3, Pcc_fig_height], 'k', 'LineWidth',1.5);
    end

    %annotation
    text(vlns_pos(1)-1.25, 0.86*Pcc_fig_height, section_labels{folder_ind}, 'FontSize', plot_font_size);




    %tpeak figure
    figure(2)

    %plot weighted violins of all tpeak samples
    weighted_violin(all_peak_estimates, 'weights', all_weights,  'x', vlns_pos, 'facecolor', peak_interior_plot_colour, ...
        'edgecolor', peak_plot_colour, 'mc', mean_colour, 'medc', []);

    %plot boxplot of weighted means for this Pcc
    boxplot(peak_WMs, 'positions', box_pos, 'colors', peak_plot_colour, 'widths', 0.7);

    %adjust colours
    h=findobj(gca,'tag','Outliers');
    set(h, 'MarkerEdgeColor', peak_plot_colour)

    %formatting
    hold on
    legend off

    %plot divider
    if folder_ind>1
        plot((vlns_pos(1)-1.5)*[1,1], [0, max(peak_prior_samples)], 'k:');
    else
        plot((vlns_pos(1)-1.5)*[1,1], [0, max(peak_prior_samples)], 'k', 'LineWidth',1.5);
    end
end



%% formatting

%plot true Pcc value
figure(1)
for Pcc_ind = 1:length(actual_Pcc)
    plot([(Pcc_ind-1)*(num_iters_to_plot+1)+0.5, Pcc_ind*(num_iters_to_plot+1)+0.5], ...
        actual_Pcc(Pcc_ind)*[1, 1], 'k');
end

%formatting
hAxes = gca;
hAxes.TickLabelInterpreter = 'latex';

xticks([])
yticks(0:0.2:1)

ylabel('$\hat{P}_{\mathrm{CC}}$')
xlim([-0.5, tot_num_obj-0.5])
ylim([-0.2, Pcc_fig_height])


%plot true tpeak value
figure(2)
for peak_ind = 1:length(actual_peak)
    plot([(peak_ind-1)*(num_iters_to_plot+1)+0.5, peak_ind*(num_iters_to_plot+1)+0.5], ...
        actual_peak(peak_ind)*[1, 1], 'k');
end

%formatting
hAxes = gca;
hAxes.TickLabelInterpreter = 'latex';

xticks(0:(tot_num_obj+1))
xticklabels(axis_labels)
xtickangle(45)

ylabel('$\hat{t}_{\mathrm{peak}}$ (h)')
xlim([-0.5, tot_num_obj-0.5])
ylim([0, max(peak_prior_samples)])




%% put in one plot

% put both plots into one figure
figlist=get(groot,'Children');

newfig=figure;
tcl=tiledlayout(newfig,'flow', 'TileSpacing','Tight','Padding','Tight');

for i = 1:numel(figlist)
    figure(figlist(i));
    ax=gca;
    ax.PlotBoxAspectRatioMode='manual';
    ax.PlotBoxAspectRatio=[fig_ratio,1,1];
    ax.FontSize = round(plot_font_size);
    ax.Parent=tcl;
    ax.Layout.Tile=numel(figlist)+1-i;
end

close 1
close 2



%fix aspect ratio
fig = gcf;
fig_pos = fig.Position;
fig_pos(4) = round(fig_pos(3)/1.7);
fig.Position = fig_pos;




