

%Makes Figure 5 from the manuscript (and also SI Figure 3)



%set font of plots
set(0,'defaulttextinterpreter','latex')

%add folders
addpath helper_funcs


%% simulation data specs

%make a plot for each of these Pcc values
Pcc_vals_to_plot = [0.1, 0.5, 0.9];

%folder setup
sample_sizes = [50, 100, 500, 1000, NaN];
base_folder = '../SIM_DATA_from_ms/SIM_EST_SWEEP_with_clustering/';

template_1 = 'fits_Pcc_';
template_2 = '_sample_';
template_3 = '_cells';


%dimensions
num_samples_per_sim = 500;
num_iters_to_plot = 4;
tot_iters = 4;


%formatting
plot_font_size = 9;
Pcc_fig_height = 1.5;


%labels
folder_labels = {'50', '100', '500', '1000', 'all'};
division_labels = {'50 Cells', '100 Cells', '500 Cells', '1000 Cells', 'All Cells (2500)'};
plot_labels = {'Prior'};
for i = 1:length(folder_labels)
    plot_labels{(i-1)*(num_iters_to_plot+1)+2} = 'Post. Wt. Means';
    for j = 1:num_iters_to_plot
        plot_labels{(i-1)*(num_iters_to_plot+1)+2+j} = strcat(['Replicate ', num2str(j)]);
    end
end


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


%% loop over Pcc values
for Pcc_ind = 1:length(Pcc_vals_to_plot)
    Pcc_to_plot = Pcc_vals_to_plot(Pcc_ind);


    %% load and plot priors
    prior_source_folder = 'helper_funcs/Pcc_and_peak_prior_densities';
    load(strcat(prior_source_folder, '/peak_prior_samples.mat'));
    load(strcat(prior_source_folder, '/Pcc_prior_samples.mat'));
    
    
    %plot prior violin for Pcc
    figure((Pcc_ind-1)*2+1)
    violin(Pcc_prior_samples', 'x', [0,1], 'facecolor', prior_interior_plot_colour, ...
        'edgecolor', prior_plot_colour, 'mc', mean_colour, 'medc', med_colour);
    legend off
    hold on
    
    %plot prior violin for tpeak
    figure((Pcc_ind-1)*2+2)
    violin(peak_prior_samples', 'x', [0,1], 'facecolor', prior_interior_plot_colour, ...
        'edgecolor', prior_plot_colour, 'mc', mean_colour, 'medc', med_colour);
    legend off
    hold on
    
    
    %% plot violins for estimates

    %total number of objects to plot
    tot_num_obj = length(folder_labels)*(num_iters_to_plot+1)+1;
    
    %initialise
    actual_Pcc = 0;
    actual_peak = 0;
    
    %loop over sample sizes
    for folder_ind = 1:length(sample_sizes)
    
        %find relevant folder
        source_folder = strcat(base_folder, '/', template_1, num2str(Pcc_to_plot), ...
            template_2, folder_labels{folder_ind}, template_3);
        
        
        %% initialise
        all_Pcc_estimates = zeros(num_samples_per_sim, num_iters_to_plot);
        all_peak_estimates = zeros(num_samples_per_sim, num_iters_to_plot);
        all_weights = zeros(num_samples_per_sim, num_iters_to_plot);
    
        actual_Pcc_estimates = zeros(1,tot_iters);
        actual_peak_estimates = zeros(1,tot_iters);
    
    
        %% load weighted means
        load(strcat(source_folder, '/Pcc_WMs'));
        load(strcat(source_folder, '/peak_WMs'));
    
    
        %% compute mean actual Pcc and peak from the simulations (across all sample sizes)
        for iter = 1:tot_iters
            load(strcat(source_folder, '/replicate_', num2str(iter), '/Pcc_estimates'));
            load(strcat(source_folder, '/replicate_', num2str(iter), '/peak_estimates'));
    
            actual_Pcc_estimates(iter) = mean(Pcc_estimates);
            actual_peak_estimates(iter) = mean(peak_estimates);
        end
    
        actual_Pcc = actual_Pcc + mean(actual_Pcc_estimates)/length(sample_sizes);
        actual_peak = actual_peak + mean(actual_peak_estimates)/length(sample_sizes);
    
    
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
        
        
        
        %% plot posterior violins and weighted mean boxplots

        %positions of plots
        vlns_pos = (folder_ind-1)*(num_iters_to_plot+1) + 1 + (1:num_iters_to_plot);
        box_pos_l = (folder_ind-1)*(num_iters_to_plot+1) + 1;
    
        %Pcc figure
        figure((Pcc_ind-1)*2+1)

        %plot weighted violins of all Pcc samples
        weighted_violin(all_Pcc_estimates, 'weights', all_weights, 'x', vlns_pos, ...
            'facecolor', Pcc_interior_plot_colour, 'edgecolor', Pcc_plot_colour, 'mc', mean_colour, 'medc', []);
        
        %plot boxplot of weighted means for this sample size 
        boxplot(Pcc_WMs, 'positions', box_pos_l, 'colors', Pcc_plot_colour, 'widths', 0.7);

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
        

        %tpeak figure
        figure((Pcc_ind-1)*2+2)

        %plot weighted violins of all tpeak samples
        weighted_violin(all_peak_estimates, 'weights', all_weights, 'x', vlns_pos, 'facecolor', peak_interior_plot_colour, ...
            'edgecolor', peak_plot_colour, 'mc', mean_colour, 'medc', []);

        %plot boxplot of weighted means for this sample size
        boxplot(peak_WMs, 'positions', box_pos_l, 'colors', peak_plot_colour, 'widths', 0.7);

        %adjust colours
        h=findobj(gca,'tag','Outliers');
        set(h, 'MarkerEdgeColor', peak_plot_colour)

        %formatting
        hold on
        legend off
    
        %plot divider
        if folder_ind>1
            plot((vlns_pos(1)-1.5)*[1,1], [min(peak_prior_samples), max(peak_prior_samples)], 'k:');
        else
            plot((vlns_pos(1)-1.5)*[1,1], [min(peak_prior_samples), max(peak_prior_samples)], 'k', 'LineWidth',1.5);
        end

        %annotation
        text(vlns_pos(1)-1.2, 0.92*max(peak_prior_samples), division_labels{folder_ind}, 'FontSize', plot_font_size)

    end
    
    
    %% formatting
    
    %plot true Pcc value
    figure((Pcc_ind-1)*2+1)
    plot([0.5, tot_num_obj+2], actual_Pcc*[1, 1], 'k');
    
    %formatting
    hAxes = gca;
    hAxes.TickLabelInterpreter = 'latex';
    
    xticks(0:(tot_num_obj+1))
    xticklabels(plot_labels)
    xtickangle(60)

    yticks(0:0.2:1)

    xlim([-0.5, tot_num_obj-0.5])
    ylim([-0.2, 1.3])

    
    %plot true tpeak value
    figure((Pcc_ind-1)*2+2)
    plot([0.5, tot_num_obj+2], actual_peak*[1,1], 'k');
    
    %formatting
    hAxes = gca;
    hAxes.TickLabelInterpreter = 'latex';
    
    xticks(0:(tot_num_obj+1))
    xticklabels(plot_labels)
    xtickangle(45)

    xlim([-0.5, tot_num_obj-0.5])
    ylim([min(peak_prior_samples), max(peak_prior_samples)])


end





