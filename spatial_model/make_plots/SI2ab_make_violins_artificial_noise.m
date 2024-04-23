

%Makes SI Figures 2ab



%set font of plots
set(0, 'DefaultTextInterpreter', 'latex')


%% setup

%paths
addpath ../sim_est_helpers/peak_and_Pcc_lookup/
addpath ../sim_est_helpers/

%data structure (includes reference case with no artificial noise)
folder_stem = "../SIM_DATA_from_ms/SIM_EST_SWEEP_artificial_noise";
ref_data_folder = '../SIM_DATA_from_ms/SIM_EST_SWEEP_with_clustering/fits_Pcc_0.5_sample_all_cells';

all_folders = ["phi_1e1", "phi_1e2", "phi_1e3", "phi_1e4", "REF"];
all_labels = {"$\phi=10^1$", "$\phi=10^2$", "$\phi=10^3$", "$\phi=10^4$", "$\phi\rightarrow\infty$"};

folders_to_use = [1, 2, 3, 4, 5];

folder_ids = all_folders(folders_to_use);
division_labels = all_labels(folders_to_use);


%data specs
num_samples_per_sim = 500;
num_iters_to_plot = 4;
tot_num_obj = length(folder_ids)*(num_iters_to_plot+1)+1;
highlight_noise_level = 5; %0 is none, counts from the left



%colour scheme
mean_colour = 'k:';
med_colour = 'k';

prior_plot_colour = [0, 0, 0];
prior_interior_plot_colour = [0.4, 0.4, 0.4];

Pcc_plot_colour = [0, 0.4470, 0.7410]; 
peak_plot_colour = [49/256, 153/256, 2/256];

Pcc_interior_plot_colour = Pcc_plot_colour*1.2;
Pcc_interior_plot_colour(Pcc_interior_plot_colour>1)=1;

peak_interior_plot_colour = peak_plot_colour*1.2;
peak_interior_plot_colour(peak_interior_plot_colour>1)=1;



%xtick labels
plot_labels = {'Prior'};
for i = 1:length(folder_ids)
    plot_labels{(i-1)*(num_iters_to_plot+1)+2} = 'Post. Wt. Means';
    for j = 1:num_iters_to_plot
        plot_labels{(i-1)*(num_iters_to_plot+1)+2+j} = strcat(['Replicate ', num2str(j)]);
    end
end


%% load priors
prior_source_folder = 'helper_funcs/Pcc_and_peak_prior_densities';
load(strcat(prior_source_folder, '/peak_prior_samples.mat'));
load(strcat(prior_source_folder, '/Pcc_prior_samples.mat'));


% plot prior violins

%Pcc
figure(1)
violin(Pcc_prior_samples', 'x', [0,1], 'facecolor', prior_interior_plot_colour, ...
    'edgecolor', prior_plot_colour, 'mc', mean_colour, 'medc', med_colour);
legend off
hold on

%t_peak
figure(2)
violin(peak_prior_samples', 'x', [0,1], 'facecolor', prior_interior_plot_colour, ...
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
        [0, 0, max(peak_prior_samples), max(peak_prior_samples)], ...
        [0.9, 0.9, 0.9], 'EdgeColor', 'None');
    hold on
end



%% loop
for folder_ind = 1:length(folder_ids)

    folder_id = folder_ids(folder_ind);

    %catch reference case (no additional noise)
    if strcmp(folder_ids(folder_ind), 'REF')
        source_folder = ref_data_folder;
    else
        source_folder = strcat(folder_stem, '/', folder_id);
    end
    
    
    %initialise
    all_Pcc_estimates = zeros(num_samples_per_sim, num_iters_to_plot);
    all_peak_estimates = zeros(num_samples_per_sim, num_iters_to_plot);
    all_weights = zeros(num_samples_per_sim, num_iters_to_plot);
    
    
    %% extract Pcc and t_peak estimates
    for iter = 1:num_iters_to_plot
        load(strcat(source_folder, '/replicate_', num2str(iter), '/final_Pcc_samples'));
        load(strcat(source_folder, '/replicate_', num2str(iter), '/final_peak_samples'));
        load(strcat(source_folder, '/replicate_', num2str(iter), '/curr_weights'));
    
        all_Pcc_estimates(:,iter)=final_Pcc_samples;
        all_peak_estimates(:,iter)=final_peak_samples;
        all_weights(:,iter)=curr_weights;
    
    end
        
    %load weighted means
    load(strcat(source_folder, '/Pcc_WMs.mat'));
    load(strcat(source_folder, '/peak_WMs.mat'));

    
    %% compute actual Pcc and peak
    for iter = 1:num_iters_to_plot
        load(strcat(source_folder, '/replicate_', num2str(iter), '/Pcc_estimates'));
        load(strcat(source_folder, '/replicate_', num2str(iter), '/peak_estimates'));

        actual_Pcc_estimates(iter) = mean(Pcc_estimates);
        actual_peak_estimates(iter) = mean(peak_estimates);
    end

    actual_Pcc = mean(actual_Pcc_estimates);
    actual_peak = mean(actual_peak_estimates);


    %% plot posterior violins

    %plot positions
    vlns_pos = (folder_ind-1)*(num_iters_to_plot+1) + 1 + (1:num_iters_to_plot);
    box_pos_l = (folder_ind-1)*(num_iters_to_plot+1) + 1;
    

    %Pcc
    figure(1)

    %plot weighted violin
    weighted_violin(all_Pcc_estimates, 'weights', all_weights, 'x', vlns_pos, ...
        'facecolor', Pcc_interior_plot_colour, 'edgecolor', Pcc_plot_colour, ...
        'mc', mean_colour, 'medc', med_colour);

    %plot box plot
    boxplot(Pcc_WMs, 'positions', box_pos_l, 'colors', Pcc_plot_colour, 'widths', 0.7);
    h=findobj(gca,'tag','Outliers');
    set(h, 'MarkerEdgeColor', Pcc_plot_colour)

    %plot true Pcc
    plot([vlns_pos(1)-1.5, vlns_pos(end)+0.5], actual_Pcc*[1, 1], 'k');
    
    
    hold on
    legend off


    %plot divider
    if folder_ind>1
        plot((vlns_pos(1)-1.5)*[1,1], [-0.3, 1.3], 'k:');
    else
        plot((vlns_pos(1)-1.5)*[1,1], [-0.3, 1.3], 'k', 'LineWidth',1.5);
    end

    %plot label
    text(vlns_pos(1)-1.2, 0.92*1.3, division_labels{folder_ind}, 'FontSize', 11);
    
    

    %t_peak
    figure(2)

    %plot weighted violin
    weighted_violin(all_peak_estimates, 'weights', all_weights, 'x', vlns_pos, ...
        'facecolor', peak_interior_plot_colour, 'edgecolor', peak_plot_colour, ...
        'mc', mean_colour, 'medc', med_colour);

    %plot box plot
    boxplot(peak_WMs, 'positions', box_pos_l, 'colors', peak_plot_colour, 'widths', 0.7);
    h=findobj(gca,'tag','Outliers');
    set(h, 'MarkerEdgeColor', peak_plot_colour)

    %plot true t_peak
    plot([vlns_pos(1)-1.5, vlns_pos(end)+0.5], actual_peak*[1, 1], 'k');


    hold on
    legend off



    %plot divider
    if folder_ind>1
        plot((vlns_pos(1)-1.5)*[1,1], [0, max(peak_prior_samples)], 'k:');
    else
        plot((vlns_pos(1)-1.5)*[1,1], [0, max(peak_prior_samples)], 'k', 'LineWidth',1.5);
    end

    %plot label
    text(vlns_pos(1)-1.2, 0.92*max(peak_prior_samples), division_labels{folder_ind}, 'FontSize', 11)

    

end


%% formatting


%Pcc plot

figure(1)

hAxes = gca;
hAxes.TickLabelInterpreter = 'latex';

xticks(0:(tot_num_obj+1))
xticklabels(plot_labels)
xtickangle(45)

yticks(0:0.2:1)

ylabel('$\hat{P}_{\mathrm{CC}}$')
xlim([-0.5, tot_num_obj-0.5])
ylim([-0.2, 1.3])


%t_peak plot
figure(2)

hAxes = gca;
hAxes.TickLabelInterpreter = 'latex';

xticks(0:(tot_num_obj+1))
xticklabels(plot_labels)
xtickangle(45)

ylabel('$\hat{t}_{\mathrm{peak}}$ (h)')
xlim([-0.5, tot_num_obj-0.5])
ylim([0, max(peak_prior_samples)])

