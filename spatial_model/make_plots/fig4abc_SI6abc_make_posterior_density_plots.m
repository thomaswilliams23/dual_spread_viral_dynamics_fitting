

%Makes Figures 4abc from the manuscript



%set font of plots
set(0,'defaulttextinterpreter','latex')

%add folders
addpath helper_funcs



%% simulation data specs

%Pcc value to plot: choose one of the below
true_Pcc = 0.1;
%true_Pcc = 0.5;
%true_Pcc = 0.9;

%whether or not clustering data was used for fit (if 0, generates SI Figs 2abc)
cluster_data = 1;

%folder for data
if cluster_data
    target_folder = strcat('../SIM_DATA_from_ms/SIM_EST_SWEEP_with_clustering/fits_Pcc_', num2str(true_Pcc), ...
        '_sample_all_cells/replicate_1');
else
    target_folder = strcat('../SIM_DATA_from_ms/SIM_EST_SWEEP_fluoro_only/fits_Pcc_', num2str(true_Pcc), ...
        '_fluoro_only/replicate_1');
end


%Pcc and parameter values for all data
all_Pcc = [0.1, 0.5, 0.9];
Pcc_names = {'01', '05', '09'};
all_alpha = [1.88035441e-01, 1.10716665e+00, 8.08312311e+00];
all_beta =  [9.61921398e-07, 3.91381868e-07, 5.74722863e-08];


%regions to show in insets - if using the fit with clustering data
has_insets = cluster_data;
x_inset_lims = [0,1; 0.2, 2.2; 5.9, 10];
y_inset_lims = [0.6e-6, 1.25e-6; 0.2e-6, 0.65e-6; 0, 0.15e-6];


%extract information relevant to the Pcc value to be plotted
Pcc_ind = find(all_Pcc == true_Pcc);
true_alpha = all_alpha(Pcc_ind);
true_beta = all_beta(Pcc_ind);
x_inset_lims = x_inset_lims(Pcc_ind,:);
y_inset_lims = y_inset_lims(Pcc_ind,:);


%% load in data
load(strcat(target_folder, '/curr_samples'));
alpha_beta_samples = curr_samples;

load(strcat(target_folder, '/curr_weights'));
alpha_beta_weights = curr_weights;


%% plot density
%do this by binning accepted samples and tallying total weight per bin

%number of bins on each axis
num_bins_per_dim = 120; %50;

%set the edges of the bins
alpha_inc = 10/(num_bins_per_dim-2);
alpha_edges = alpha_inc*(0:num_bins_per_dim);

beta_inc = 1.5e-6/(num_bins_per_dim-2);
beta_edges = beta_inc*(0:num_bins_per_dim);

%initialise
hist_counts = zeros(num_bins_per_dim);


%loop over accepted samples
for smpl = 1:length(alpha_beta_weights)

    %find alpha_ind
    alpha_ind = 1;
    while alpha_beta_samples(1,smpl)>alpha_edges(alpha_ind+1)
        alpha_ind = alpha_ind + 1;
    end

    %find beta_ind
    beta_ind = 1;
    while alpha_beta_samples(2,smpl)>beta_edges(beta_ind+1)
        beta_ind = beta_ind + 1;
    end

    %add weight to that element in array
    hist_counts(beta_ind, alpha_ind) = hist_counts(beta_ind, alpha_ind) + alpha_beta_weights(smpl);

end

%normalise
patch_size = (alpha_edges(2)-alpha_edges(1))*(beta_edges(2)-beta_edges(1));
hist_counts = hist_counts/(sum(sum(hist_counts))*patch_size);


%now plot
PLOT_THRESH = 1e-4;
figure
hold on
[~,contourObj] = contourf(alpha_edges(1:(end-1)), beta_edges(1:(end-1)), hist_counts, 'EdgeAlpha',0);
contourObj.LevelList = [0, PLOT_THRESH, contourObj.LevelList(2:end)];
colormap parula


%consistent colour limits
if cluster_data
    %for cluster fit
    clim([0, 1.2e7]);

else
    %for fluoro fit
    clim([0, 2e6]);
end

%hide colourbar
colorbar('visible', 'off');



%make base contour transparent (from WILL GRANT on MATLAB answers):
% This is the secret that 'keeps' the transparency.
eventFcn = @(srcObj, e) updateTransparency(srcObj);
addlistener(contourObj, 'MarkedClean', eventFcn);
%




%% plot true parameter values
true_val = plot(true_alpha, true_beta, 'ro', 'MarkerSize', 15);


%% add contours

%% first load in data

%alpha and beta range for the overall contour
load('../SIM_DATA_from_ms/lookup_tables/alpha_range.mat') %contains proportions of max used for lookup table (e.g. 0, 0.1, 0.2, etc)
load('../SIM_DATA_from_ms/lookup_tables/beta_range.mat')

alpha_max = 10;
beta_max = 5e-6;

alpha_range_outer = alpha_max*alpha_range;
beta_range_outer = beta_max*beta_range;


%gets alpha and beta range for inset
lookup_folder = strcat('../SIM_DATA_from_ms/lookup_tables/insets/Pcc_', Pcc_names{Pcc_ind});
load(strcat(lookup_folder, '/alpha_inc'));
load(strcat(lookup_folder, '/beta_inc'));
load(strcat(lookup_folder, '/alpha_offset'));
load(strcat(lookup_folder, '/beta_offset'));

alpha_range_inner = alpha_offset + alpha_inc * (0:0.1:1);
beta_range_inner = beta_offset + beta_inc * (0:0.1:1);


%% plot the true peak contour

%compute true peak time from lookup table
true_peak_time = query_peak_lookup_table('../SIM_DATA_from_ms/lookup_tables/tpeak_lookup', ...
    alpha_range_outer, beta_range_outer, true_alpha, true_beta);


%%%%% BUILD A CONTOUR FOR THE WHOLE PLOT

%initialise
beta_on_peak_contour_outer = zeros(1,length(alpha_range_outer));

%compute corresponding beta value for each alpha value in the range
for beta_ind = 1:length(alpha_range_outer)
    beta_on_peak_contour_outer(beta_ind) = infer_beta_on_peak_time_contour('../SIM_DATA_from_ms/lookup_tables/tpeak_lookup', ...
        alpha_range_outer, beta_range_outer, alpha_range_outer(beta_ind), true_peak_time);
end



%%%%% BUILD A CONTOUR FOR THE INSET REGION

%initialise
beta_on_peak_contour_inner = zeros(1,length(alpha_range_inner));

%compute corresponding beta value for each alpha value in the range
% NOTE!! these lookup tables were subject to Gaussian smoothing with 
% sigma=1 due to noise at small region of parameter space
for beta_ind = 1:length(alpha_range_inner)
    beta_on_peak_contour_inner(beta_ind) = infer_beta_on_peak_time_contour(strcat(lookup_folder, '/peak_inset_lookup_SMOOTHED'), ...
        alpha_range_inner, beta_range_inner, alpha_range_inner(beta_ind), true_peak_time);
end


%%%%% now join up the two contour sections

%find left and right edges of the inset region
inset_min = min(alpha_range_inner);
inset_max = max(alpha_range_inner);

outer_cutoff_L = find(alpha_range_outer<inset_min, 1, 'last');
outer_cutoff_R = find(alpha_range_outer>inset_max, 1, 'first');

%assemble overall contour
all_alpha_vals_peak = [alpha_range_outer(1:outer_cutoff_L), ...
    alpha_range_inner, alpha_range_outer(outer_cutoff_R:end)];
all_beta_vals_peak = [beta_on_peak_contour_outer(1:outer_cutoff_L), ...
    beta_on_peak_contour_inner, beta_on_peak_contour_outer(outer_cutoff_R:end)];

%and plot
plot(all_alpha_vals_peak, all_beta_vals_peak, 'r--', 'LineWidth',1);



%% Pcc contour

%%%%% BUILD A CONTOUR FOR THE WHOLE PLOT

%initialise
alpha_on_Pcc_contour_outer = zeros(1,length(beta_range_outer));

%compute corresponding alpha value for each beta value in the range
% NOTE: we omit the first few points in the range as Pcc is undefined near
% the origin since simulations don't infect all cells
for beta_ind = 4:length(beta_range_outer)
    alpha_on_Pcc_contour_outer(beta_ind) = infer_alpha_on_Pcc_contour('../SIM_DATA_from_ms/lookup_tables/Pcc_lookup', ...
        alpha_range_outer, beta_range_outer, beta_range_outer(beta_ind), true_Pcc);
end


%%%%% BUILD A CONTOUR FOR THE INSET REGION

%initialise
alpha_on_Pcc_contour_inner = zeros(1,length(beta_range_inner));

%compute corresponding alpha value for each beta value in the range
for beta_ind = 1:length(beta_range_inner)
    alpha_on_Pcc_contour_inner(beta_ind) = infer_alpha_on_Pcc_contour(strcat(lookup_folder, '/Pcc_inset_lookup_SMOOTHED'), ...
        alpha_range_inner, beta_range_inner, beta_range_inner(beta_ind), true_Pcc);
end


%%%%% now join up two contour regions
if ~Pcc_ind ==3

    %find left and right edges of the inset region
    inset_min = min(beta_range_inner);
    inset_max = max(beta_range_inner);
    
    outer_cutoff_L = find(beta_range_outer<inset_min, 1, 'last');
    outer_cutoff_R = find(beta_range_outer>inset_max, 1, 'first');
    
    %assemble
    all_alpha_vals_Pcc = [alpha_on_Pcc_contour_outer(1:outer_cutoff_L), ...
        alpha_on_Pcc_contour_inner, alpha_on_Pcc_contour_outer(outer_cutoff_R:end)];
    all_beta_vals_Pcc = [beta_range_outer(1:outer_cutoff_L), ...
        beta_range_inner, beta_range_outer(outer_cutoff_R:end)];

else %%% SPECIAL CASE: for the Pcc=0.9 case, we include more beta points 
     %%% from the outer contour to smooth the join between the two regions

    %find left and right edges of the inset region
    inset_min = beta_range_inner(4);
    inset_max = max(beta_range_inner);
    
    outer_cutoff_L = find(beta_range_outer<inset_min, 1, 'last');
    outer_cutoff_R = find(beta_range_outer>inset_max, 1, 'first');
    
    %assemble
    all_alpha_vals_Pcc = [alpha_on_Pcc_contour_outer(1:outer_cutoff_L), ...
        alpha_on_Pcc_contour_inner(4:end), alpha_on_Pcc_contour_outer(outer_cutoff_R:end)];
    all_beta_vals_Pcc = [beta_range_outer(1:outer_cutoff_L), ...
        beta_range_inner(4:end), beta_range_outer(outer_cutoff_R:end)];

end


%and plot
plot(all_alpha_vals_Pcc, all_beta_vals_Pcc, 'r--', 'LineWidth',1);


%% format
x_max = 10;
y_max = 1.5e-6;

xlim([0, x_max])
ylim([0, y_max])

xlabel('$\alpha$')
ylabel('$\beta$')


%% make inset plots if needed
if has_insets

    %add a box around inset region
    p = patch([x_inset_lims, fliplr(x_inset_lims)], ...
        [y_inset_lims(1), y_inset_lims, y_inset_lims(2)],...
        'k', 'FaceAlpha', 0, 'EdgeColor', [0,0,0]);



    %% now make a new figure and plot everything again
    figure
    hold on
    [~,contourObj] = contourf(alpha_edges(1:(end-1)), beta_edges(1:(end-1)), hist_counts, 'EdgeAlpha',0);
    contourObj.LevelList = [0, PLOT_THRESH, contourObj.LevelList(2:end)];
    colormap parula
    if cluster_data
        clim([0, 1.2e7]);
    else
        clim([0, 2e6]);
    end
    colorbar('visible', 'off');
    
    %make base contour transparent (from WILL GRANT on MATLAB answers):
    % This is the secret that 'keeps' the transparency.
    eventFcn = @(srcObj, e) updateTransparency(srcObj);
    addlistener(contourObj, 'MarkedClean', eventFcn);
    %

    true_val = plot(true_alpha, true_beta, 'ro', 'MarkerSize', 15);
    plot(all_alpha_vals_Pcc, all_beta_vals_Pcc, 'r--', 'LineWidth',1);
    plot(all_alpha_vals_peak, all_beta_vals_peak, 'r--', 'LineWidth',1);
    p = patch([x_inset_lims, fliplr(x_inset_lims)], ...
        [y_inset_lims(1), y_inset_lims, y_inset_lims(2)],...
        'k', 'FaceAlpha', 0, 'EdgeColor', [0,0,0]);


    %% compute aspect ratio of inset
    fig_ratio = 1.1;
    norm_x_width = (x_inset_lims(2)-x_inset_lims(1))/(x_max*fig_ratio);
    norm_y_width = (y_inset_lims(2)-y_inset_lims(1))/y_max;
    inset_aspect_ratio = norm_x_width/norm_y_width;

    %% reformat
    %now change the x and y limits and make the inset figure
    xlim(x_inset_lims);
    ylim(y_inset_lims);

    xticks([]);
    yticks([]);

    xlabel([]);
    ylabel([]);

    title([]);

    box off
    grid off

    ax = gca;
    asp_ratio_vec = ax.PlotBoxAspectRatio;
    asp_ratio_vec(2) = 1/inset_aspect_ratio;
    ax.PlotBoxAspectRatio = asp_ratio_vec;
end




%% from WILL GRANT on MATLAB answers:
% Elsewhere in script, a separate file, or another method of your class.
function updateTransparency(contourObj)
contourFillObjs = contourObj.FacePrims;
for i = 1%:length(contourFillObjs)
    % Have to set this. The default is 'truecolor' which ignores alpha.
    contourFillObjs(i).ColorType = 'truecoloralpha';
    % The 4th element is the 'alpha' value. First 3 are RGB. Note, the
    % values expected are in range 0-255.
    contourFillObjs(i).ColorData(4) = 0;
end
end
