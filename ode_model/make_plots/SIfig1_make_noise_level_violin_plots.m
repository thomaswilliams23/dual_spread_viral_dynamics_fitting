

%Makes SI Figure 1



%set font of plots
set(0,'defaulttextinterpreter','latex')

%add folders
addpath helper_funcs



%% simulation data specs

%folder template for data
folder_stem = "../SIM_DATA_from_ms/SIM_EST_SWEEP/";

%folder for prior samples
prior_directory = 'helper_funcs/Pcc_and_r_prior_densities';

%true Pcc value of data to plot
true_Pcc = 0.5;
fname = 'Pcc_05/';

%all noise levels available to plot
all_folders = ["phi_1e1", "phi_1e2", "phi_1e3", "phi_1e4", "phi_1e5"];
all_labels = {"$\phi=10^1$", "$\phi=10^2$", "$\phi=10^3$", "$\phi=10^4$", "$\phi=10^5$"};

%noise levels to show in plot (indices)
folders_to_use = [1, 2, 3, 4];

%truncate folder names and labels to ones used in the plot
folder_ids = all_folders(folders_to_use);
division_labels = all_labels(folders_to_use);

%formatting
plot_font_size = 9;
fig_width = 18; 
fig_ratio = 3.5;

%data specs
num_samples_per_sim = 4000;
num_iters_to_plot = 4;
tot_iters = 10;

%number of objects to plot
tot_num_obj = length(folder_ids)*(num_iters_to_plot+1)+1;

%highlights a particular section of the plot
highlight_noise_level = 2; %<-- 0 is none, counts from the left

%true value of r
true_r = 0.52;


%colours
mean_colour = 'k:';
med_colour = 'k';

prior_plot_colour = [0, 0, 0];
prior_interior_plot_colour = [0.4, 0.4, 0.4];

Pcc_plot_colour = [0, 0.4470, 0.7410]; 
r_plot_colour = [0.8500, 0.3250, 0.0980];

Pcc_interior_plot_colour = Pcc_plot_colour*1.2;
Pcc_interior_plot_colour(Pcc_interior_plot_colour>1)=1;

r_interior_plot_colour = r_plot_colour*1.2;
r_interior_plot_colour(r_interior_plot_colour>1)=1;


%labels
plot_labels = {'Prior'};
for i = 1:length(folder_ids)
    plot_labels{(i-1)*(num_iters_to_plot+1)+2} = 'Post. Meds';
    for j = 1:num_iters_to_plot
        plot_labels{(i-1)*(num_iters_to_plot+1)+2+j} = strcat(['Replicate ', num2str(j)]);
    end
end


%% load priors
load(strcat(prior_directory, '/r_prior_samples.mat'));
load(strcat(prior_directory, '/Pcc_prior_samples.mat'));


% plot prior violin for Pcc
figure(1)
violin(Pcc_prior_samples', 'x', [0,1], 'facecolor', prior_interior_plot_colour, ...
    'edgecolor', prior_plot_colour, 'mc', mean_colour, 'medc', med_colour);
legend off
hold on

% plot prior violin for r
figure(2)
violin(r_prior_samples', 'x', [0,1], 'facecolor', prior_interior_plot_colour, ...
    'edgecolor', prior_plot_colour, 'mc', mean_colour, 'medc', med_colour);
legend off
hold on


%% put in highlight
if highlight_noise_level
    figure(1)
    patch([((highlight_noise_level-1)*(num_iters_to_plot+1) + 0.5),...
        ((highlight_noise_level)*(num_iters_to_plot+1) + 0.5),...
        ((highlight_noise_level)*(num_iters_to_plot+1) + 0.5),...
        ((highlight_noise_level-1)*(num_iters_to_plot+1) + 0.5)],...
        [-0.2, -0.2, 1.3, 1.3], ...
        [0.9, 0.9, 0.9], 'EdgeColor', 'None');
    hold on

    figure(2)
    patch([((highlight_noise_level-1)*(num_iters_to_plot+1) + 0.5),...
        ((highlight_noise_level)*(num_iters_to_plot+1) + 0.5),...
        ((highlight_noise_level)*(num_iters_to_plot+1) + 0.5),...
        ((highlight_noise_level-1)*(num_iters_to_plot+1) + 0.5)],...
        [0, 0, max(r_prior_samples), max(r_prior_samples)], ...
        [0.9, 0.9, 0.9], 'EdgeColor', 'None');
    hold on
end



%% loop over noise levels
for folder_ind = 1:length(folder_ids)

    folder_id = folder_ids(folder_ind);
    source_folder = strcat(folder_stem, fname, folder_id);
    
    
    %% initialise
    all_Pcc_estimates = zeros(num_samples_per_sim, num_iters_to_plot);
    all_r_estimates = zeros(num_samples_per_sim, num_iters_to_plot);
    
    
    %% loop replicates
    for iter = 1:num_iters_to_plot

        %load data
        load(strcat(source_folder, '/replicate_', num2str(iter), '/Pcc_samples'));
        load(strcat(source_folder, '/replicate_', num2str(iter), '/r_samples'));
    
        %save into big array
        all_Pcc_estimates(:,iter)=all_sims_Pcc_ppc;
        all_r_estimates(:,iter)=all_sims_r_ppc;
    
    end
    
    
    
    %% load medians
    load(strcat(source_folder, '/Pcc_medians.mat'));
    load(strcat(source_folder, '/r_medians.mat'));
    
    
    %% plot posterior violins and median boxplots

    %positions of plots
    vlns_pos = (folder_ind-1)*(num_iters_to_plot+1) + 1 + (1:num_iters_to_plot);
    box_pos_l = (folder_ind-1)*(num_iters_to_plot+1) + 1;
    
    %Pcc figure
    figure(1)

    %plot violins of all Pcc samples
    violin(all_Pcc_estimates, 'x', vlns_pos, 'facecolor', Pcc_interior_plot_colour, ...
        'edgecolor', Pcc_plot_colour, 'mc', mean_colour, 'medc', med_colour);

    %plot boxplot of all medians this noise level
    boxplot(Pcc_medians, 'positions', box_pos_l, 'colors', Pcc_plot_colour, 'widths', 0.7);

    %adjust colours
    h=findobj(gca,'tag','Outliers');
    set(h, 'MarkerEdgeColor', Pcc_plot_colour)

    %formatting
    hold on
    legend off

    %plot divider
    if folder_ind>1
        plot((vlns_pos(1)-1.5)*[1,1], [-0.3, 1.3], 'k:');
    else
        plot((vlns_pos(1)-1.5)*[1,1], [-0.3, 1.3], 'k', 'LineWidth',1.5);
    end

    %annotation
    text(vlns_pos(1)-1.2, 0.92*1.3, division_labels{folder_ind}, 'FontSize', plot_font_size);
    
    
    %r plot
    figure(2)

    %plot violins of all r samples
    violin(all_r_estimates, 'x', vlns_pos, 'facecolor', r_interior_plot_colour, ...
        'edgecolor', r_plot_colour, 'mc', mean_colour, 'medc', med_colour);

    %plot boxplot of all medians for this noise level
    boxplot(r_medians, 'positions', box_pos_l, 'colors', r_plot_colour, 'widths', 0.7);

    %adjust colours
    h=findobj(gca,'tag','Outliers');
    set(h, 'MarkerEdgeColor', r_plot_colour)

    %formatting
    hold on
    legend off

    %plot divider
    if folder_ind>1
        plot((vlns_pos(1)-1.5)*[1,1], [min(r_prior_samples), max(r_prior_samples)], 'k:');
    else
        plot((vlns_pos(1)-1.5)*[1,1], [min(r_prior_samples), max(r_prior_samples)], 'k', 'LineWidth',1.5);
    end

    %annotation
    text(vlns_pos(1)-1.2, 0.92*max(r_prior_samples), division_labels{folder_ind}, 'FontSize', plot_font_size)

    

end


%% formatting

%plot true Pcc value
figure(1)
plot([0.5, tot_num_obj+2], true_Pcc*[1, 1], 'k');

%formatting
hAxes = gca;
hAxes.TickLabelInterpreter = 'latex';

xticks(0:(tot_num_obj+1))
xticklabels(plot_labels)
xtickangle(45)

yticks(0:0.2:1)

ylabel('$\hat{P}_{\mathrm{CC}}$')
xlim([-0.5, tot_num_obj-0.5])
ylim([-0.2, 1.3])


%plot true r value
figure(2)
plot([0.5, tot_num_obj+2], true_r*[1,1], 'k');

%formatting
hAxes = gca;
hAxes.TickLabelInterpreter = 'latex';

xticks(0:(tot_num_obj+1))
xticklabels(plot_labels)
xtickangle(45)

ylabel('$\hat{r}$ (h$^{-1}$)')
xlim([-0.5, tot_num_obj-0.5])
ylim([0, max(r_prior_samples)])


