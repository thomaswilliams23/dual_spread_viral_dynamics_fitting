function [mean_neighbour_clustering, ERR_FLAG] = ...
    compute_mean_neighbour_clustering(cell_grid, num_species, num_samples_for_clustering)


    %compute distribution of num. fluorescent neighbours per fluorescent cell
    neighbour_dist = compute_neighbour_clustering_distribution(cell_grid,num_species, num_samples_for_clustering);
    
    %then compute mean
    mean_neighbour_clustering = sum((0:(length(neighbour_dist)-1)).*neighbour_dist);

    %report error flag if passed no infected cells (i.e. sample turned up no infected cells)
    if ~sum(neighbour_dist)
        ERR_FLAG = 1;
    else
        ERR_FLAG = 0;
    end

end