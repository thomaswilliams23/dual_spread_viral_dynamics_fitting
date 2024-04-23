

%Makes Figure 6 from the manuscript (and also SI Figure 8)



%set font of plots
set(0, 'DefaultTextInterpreter', 'latex')


%% setup

%paths
addpath ../sim_est_helpers/peak_and_Pcc_lookup/
addpath ../sim_est_helpers/


%data structure
Pcc_vals_to_plot = [0.1, 0.5, 0.9];
fnames = {'Pcc_01', 'Pcc_05', 'Pcc_09'};
diff_vals = [0.1, 1, 10, 100, Inf];

base_folder = '../SIM_DATA_from_ms/SIM_EST_SWEEP_diffusion_data/';
template_1 = 'fits_Pcc_';
template_2 = '_true_diff_';


%this is the reference case where data is generated with infinite diffusion
ref_fit_1 = '../SIM_DATA_from_ms/SIM_EST_SWEEP_with_clustering/fits_Pcc_';
ref_fit_2 = '_sample_all_cells';


%dimensions
num_samples_per_sim = 500;
num_iters_to_plot = 4;
tot_iters = 4;


%labels
folder_labels = {'01', '1', '10', '100', 'Inf'};
division_labels = {'$D=0.1$', '$D=1$', '$D=10$', '$D=100$', '$D\rightarrow\infty$'};
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
peak_plot_colour = [49/256, 153/256, 2/256]; %[0.8500, 0.3250, 0.0980];

Pcc_interior_plot_colour = Pcc_plot_colour*1.2;
Pcc_interior_plot_colour(Pcc_interior_plot_colour>1)=1;

peak_interior_plot_colour = peak_plot_colour*1.2;
peak_interior_plot_colour(peak_interior_plot_colour>1)=1;



%% loop
%loop Pcc
for Pcc_ind = 1:length(Pcc_vals_to_plot)

    Pcc_to_plot = Pcc_vals_to_plot(Pcc_ind);
    fname = fnames{Pcc_ind};


    %% load and plot priors
    prior_source_folder = 'helper_funcs/Pcc_and_peak_prior_densities';
    load(strcat(prior_source_folder, '/peak_prior_samples.mat'));
    load(strcat(prior_source_folder, '/Pcc_prior_samples.mat'));
    
    
    % plot prior violin
    figure((Pcc_ind-1)*2+1)
    violin(Pcc_prior_samples', 'x', [0,1], 'facecolor', prior_interior_plot_colour, ...
        'edgecolor', prior_plot_colour, 'mc', mean_colour, 'medc', med_colour);
    legend off
    hold on
    
    figure((Pcc_ind-1)*2+2)
    violin(peak_prior_samples', 'x', [0,1], 'facecolor', prior_interior_plot_colour, ...
        'edgecolor', prior_plot_colour, 'mc', mean_colour, 'medc', med_colour);
    legend off
    hold on
    

    
    %% main plot code
    tot_num_obj = length(folder_labels)*(num_iters_to_plot+1)+1;
    for folder_ind = 1:length(diff_vals)
    
        %get the source folder for this case
        if ~(diff_vals(folder_ind)==Inf)
            source_folder = strcat(base_folder, '/', template_1, num2str(Pcc_to_plot), ...
                template_2, folder_labels{folder_ind});
        else
            source_folder = strcat(ref_fit_1, num2str(Pcc_to_plot), ref_fit_2);
        end
        
        
        %% initialise
        all_Pcc_estimates = zeros(num_samples_per_sim, num_iters_to_plot);
        all_peak_estimates = zeros(num_samples_per_sim, num_iters_to_plot);
        all_weights = zeros(num_samples_per_sim, num_iters_to_plot);
    
        actual_Pcc_estimates = zeros(1,tot_iters);
        actual_peak_estimates = zeros(1,tot_iters);
    
    
        %% load weighted means
        load(strcat(source_folder, '/Pcc_WMs'));
        load(strcat(source_folder, '/peak_WMs'));
    
    
        %% compute actual Pcc and peak
        for iter = 1:tot_iters
            load(strcat(source_folder, '/replicate_', num2str(iter), '/Pcc_estimates'));
            load(strcat(source_folder, '/replicate_', num2str(iter), '/peak_estimates'));
    
            actual_Pcc_estimates(iter) = mean(Pcc_estimates);
            actual_peak_estimates(iter) = mean(peak_estimates);
        end
    
        actual_Pcc = mean(actual_Pcc_estimates);
        actual_peak = mean(actual_peak_estimates);
    
    
        %% extract data from replicates
        for iter = 1:num_iters_to_plot
            load(strcat(source_folder, '/replicate_', num2str(iter), '/final_Pcc_samples'));
            load(strcat(source_folder, '/replicate_', num2str(iter), '/final_peak_samples'));
            load(strcat(source_folder, '/replicate_', num2str(iter), '/curr_weights'));
        
            all_Pcc_estimates(:,iter)=final_Pcc_samples;
            all_peak_estimates(:,iter)=final_peak_samples;
            all_weights(:,iter)=curr_weights;
        end
        
        
        
        %% plot posterior violins and median boxplots
    
        %plot positions
        vlns_pos = (folder_ind-1)*(num_iters_to_plot+1) + 1 + (1:num_iters_to_plot);
        box_pos_l = (folder_ind-1)*(num_iters_to_plot+1) + 1;



        %Pcc
    
        %get relevant figure
        figure((Pcc_ind-1)*2+1)

        %plot weighted violin
        weighted_violin(all_Pcc_estimates, 'weights', all_weights, 'x', vlns_pos, ...
            'facecolor', Pcc_interior_plot_colour, 'edgecolor', Pcc_plot_colour, 'mc', mean_colour, 'medc', []);
        
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

        %optionally add text labels
        if Pcc_ind==1
            text(vlns_pos(1)-1.2, 0.92*1.3, division_labels{folder_ind}, 'FontSize', 11);
        end
        


        %t_peak

        %get relevant figure
        figure((Pcc_ind-1)*2+2)
        
        %plot weighted violin
        weighted_violin(all_peak_estimates, 'weights', all_weights, 'x', vlns_pos, 'facecolor', peak_interior_plot_colour, ...
            'edgecolor', peak_plot_colour, 'mc', mean_colour, 'medc', []);
        
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
            plot((vlns_pos(1)-1.5)*[1,1], [min(peak_prior_samples), max(peak_prior_samples)], 'k:');
        else
            plot((vlns_pos(1)-1.5)*[1,1], [min(peak_prior_samples), max(peak_prior_samples)], 'k', 'LineWidth',1.5);
        end

        %optionally put in text labels
        if Pcc_ind==1
            text(vlns_pos(1)-1.2, 0.92*max(peak_prior_samples), division_labels{folder_ind}, 'FontSize', 11)
        end
        
    
    end
    
    
    %% formatting
    

    %Pcc figure

    figure((Pcc_ind-1)*2+1)
    
    hAxes = gca;
    hAxes.TickLabelInterpreter = 'latex';
    
    xticks(0:(tot_num_obj+1))
    xticklabels(plot_labels)
    xtickangle(60)
    

    if ~(Pcc_ind ==3)
        xticks([]);
    end
    
    yticks(0:0.2:1)

    xlim([-0.5, tot_num_obj-0.5])
    ylim([-0.2, 1.3])
    
    


    %t_peak figure
    
    figure((Pcc_ind-1)*2+2)
    
    hAxes = gca;
    hAxes.TickLabelInterpreter = 'latex';
    
    xticks(0:(tot_num_obj+1))
    xticklabels(plot_labels)
    xtickangle(45)
    

    if ~(Pcc_ind ==3)
        xticks([]);
    end


    xlim([-0.5, tot_num_obj-0.5])
    ylim([min(peak_prior_samples), max(peak_prior_samples)])

end



