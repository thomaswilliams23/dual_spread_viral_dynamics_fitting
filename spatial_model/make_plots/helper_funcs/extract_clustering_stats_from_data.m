
% this version assumes that no simulations are ended early


%initialise
all_clustering_data = zeros(num_time_points, num_reps, length(folder_index_num));


%loop over folders
for fldr = 1:length(folder_index_num)

    fprintf('\n-----------------------------------\n')
    fprintf('\n\nWorking on folder number %d', fldr)
    
    %string version of folder name
    folder_names{fldr} = num2str(folder_index_num(fldr));
    
    
    %open data
    load(strcat(data_source, folder_stem, folder_names{fldr}, '/all_grid_data.mat'));
      
  
    %% analysis

    clustering_this_fldr = zeros(num_time_points,num_reps);

    %loop over replicates
    for rep = 1:num_reps

        fprintf('\nProcessing rep %d:\n', rep);

        %loop over time points
        for time_point = 1:num_time_points

            fprintf('\n --> Analysing time point %d', time_point);
            
            %extract the grid state at time point
            grid_to_analyse = abs(squeeze(all_grid_data(:,:,time_point,rep)));

            %compute clustering metric on grid
            clustering_this_fldr(time_point,rep) = compute_clustering_stat(grid_to_analyse);
        end
    end

    %save out
    all_clustering_data(:,:,fldr) = clustering_this_fldr;
end

