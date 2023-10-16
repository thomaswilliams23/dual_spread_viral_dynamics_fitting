

%Makes Figure 2d from the manuscript



%set font of plots
set(0,'defaulttextinterpreter','latex')

%add folders
addpath helper_funcs



%% simulation data specs

%folder for data
folder_template = "../SIM_DATA_from_ms/SIM_EST_SWEEP/Pcc_";
noise_level = "phi_1e2";

%folder for prior samples
prior_directory = 'helper_funcs/Pcc_and_r_prior_densities';

%Pcc values for all data
folder_indices = {'01', '05', '09'};
true_Pcc_vals = [0.1, 0.5, 0.9];


%label
section_labels = {{'True $P_{\mathrm{CC}}\approx$ 0.1, ','($\alpha=0.229$, $\beta=1.38\times 10^{-6}$)'},...
    {'True $P_{\mathrm{CC}}\approx$ 0.5, ','($\alpha=1.09$, $\beta=7.20\times 10^{-7}$)'},...
    {'True $P_{\mathrm{CC}}\approx$ 0.9, ','($\alpha=1.87$, $\beta=1.36\times 10^{-7}$)'}};

%formatting
plot_font_size = 9;
fig_width = 20;
fig_ratio = 4;
Pcc_fig_height = 1.5;

%data specs
num_samples_per_sim = 4000;
num_iters_to_plot = 4;
tot_iters = 10;

%number of objects to plot
tot_num_obj = length(folder_indices)*(num_iters_to_plot+1)+1;

%highlights a particular section of the plot (in this case bolds the label)
highlight_first_replicate = 1;

%true value of r
true_r = 0.52;


%colours
mean_colour = 'k:';
med_colour = 'k';

prior_plot_colour = [0, 0, 0];
prior_interior_plot_colour = [0.4, 0.4, 0.4];

Pcc_plot_colour = [0, 0.4470, 0.7410];
r_plot_colour = [0.8500, 0.3250, 0.0980]; % [0.9290, 0.6940, 0.1250];

Pcc_interior_plot_colour = Pcc_plot_colour*1.2;
Pcc_interior_plot_colour(Pcc_interior_plot_colour>1)=1;

r_interior_plot_colour = r_plot_colour*1.2;
r_interior_plot_colour(r_interior_plot_colour>1)=1;


%axis labels
axis_labels = {'Prior'};
for i = 1:length(folder_indices)
    axis_labels{(i-1)*(num_iters_to_plot+1)+2} = 'Post. Meds';

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



%% loop over Pcc values
for folder_ind = 1:length(folder_indices)

    %folder name
    source_folder = strcat(folder_template, folder_indices{folder_ind}, '/', noise_level);


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
    box_pos = (folder_ind-1)*(num_iters_to_plot+1) + 1;

    %Pcc figure
    figure(1)

    %plot violins of all Pcc samples
    violin(all_Pcc_estimates, 'x', vlns_pos, 'facecolor', Pcc_interior_plot_colour,...
        'edgecolor', Pcc_plot_colour, 'mc', mean_colour, 'medc', med_colour);

    %plot boxplot of all medians for this Pcc
    boxplot(Pcc_medians, 'positions', box_pos, 'colors', Pcc_plot_colour, 'widths', 0.7);

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


    %r plot
    figure(2)

    %plot violins of all r samples
    violin(all_r_estimates, 'x', vlns_pos, 'facecolor', r_interior_plot_colour, ...
        'edgecolor', r_plot_colour, 'mc', mean_colour, 'medc', med_colour);

    %plot boxplot of all medians for this Pcc
    boxplot(r_medians, 'positions', box_pos, 'colors', r_plot_colour, 'widths', 0.7);

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

end


%% formatting

%plot true Pcc value
figure(1)
for Pcc_ind = 1:length(true_Pcc_vals)
    plot([(Pcc_ind-1)*(num_iters_to_plot+1)+0.5, Pcc_ind*(num_iters_to_plot+1)+0.5], ...
        true_Pcc_vals(Pcc_ind)*[1, 1], 'k');
end

%formatting
hAxes = gca;
hAxes.TickLabelInterpreter = 'latex';

xticks([])
yticks(0:0.2:1)

ylabel('$\hat{P}_{\mathrm{CC}}$')
xlim([-0.5, tot_num_obj-0.5])
ylim([-0.2, Pcc_fig_height])


%plot true r value
figure(2)
plot([0.5, tot_num_obj+2.5], true_r*[1,1], 'k');

%formatting
hAxes = gca;
hAxes.TickLabelInterpreter = 'latex';

xticks(0:(tot_num_obj+1))
xticklabels(axis_labels)
xtickangle(45)

yticks(0:0.2:1)

ylabel('$\hat{r}$ (h$^{-1}$)')
xlim([-0.5, tot_num_obj-0.5])
ylim([0, max(r_prior_samples)])



%% put in one plot

%put both plots into one figure
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

