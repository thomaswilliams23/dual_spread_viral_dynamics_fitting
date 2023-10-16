function number_same_neighbours_prob = ...
    compute_neighbour_clustering_distribution(cell_grid, num_species, num_samples_for_clustering)

%sizes
cells_wide = size(cell_grid,1);
cells_long = size(cell_grid,2);

total_cells = numel(cell_grid);



%work out the target number of cells to analyse: either ALL the infected
%cells, or limited to a specified number
COUNT_ALL_INF_CELLS_FLAG = -1;
if num_samples_for_clustering == COUNT_ALL_INF_CELLS_FLAG
    num_cells_to_analyse = total_cells;
else
    num_cells_to_analyse = num_samples_for_clustering;
end


%initialising
all_cell_inds_shuffled = randperm(total_cells, num_cells_to_analyse);

neighbour_clustering_all_cells = zeros(num_cells_to_analyse, 1);


%loop
num_inf_or_dead = 0;
is_inf_or_dead = zeros(num_cells_to_analyse, 1);
for sample_ind=1:num_cells_to_analyse

    %draw random cell_ind WITHOUT replacement
    cell_ind = all_cell_inds_shuffled(sample_ind);

    %compute cell_i and cell_j from cell_ind
    [cell_i, cell_j] = ind2sub(size(cell_grid), cell_ind);

    %get the species at these indices
    species_this_cell = get_infected_species_incl_dead(cell_grid(cell_i, cell_j), num_species);

    %go to the next round of the loop if cell not infected or dead
    if ~species_this_cell
        continue
    end

    %increment count of inf/dead cells
    num_inf_or_dead = num_inf_or_dead + 1;
    is_inf_or_dead(cell_ind) = 1;

    %otherwise continue
    if mod(cell_i,2) %odd row, left-aligned

        neighbour_species = ...
            [get_infected_species_incl_dead(cell_grid( mod(cell_i-2,cells_wide)+1, mod(cell_j-2,cells_long)+1 ),num_species), ...
            get_infected_species_incl_dead(cell_grid( mod(cell_i-2,cells_wide)+1, cell_j ),num_species), ...
            get_infected_species_incl_dead(cell_grid( cell_i, mod(cell_j-2,cells_long)+1 ),num_species), ...
            get_infected_species_incl_dead(cell_grid( cell_i, mod(cell_j,cells_long)+1 ),num_species), ...
            get_infected_species_incl_dead(cell_grid( mod(cell_i,cells_wide)+1, mod(cell_j-2,cells_long)+1 ),num_species), ...
            get_infected_species_incl_dead(cell_grid( mod(cell_i,cells_wide)+1, cell_j),num_species)];

    else %even row, right-aligned

        neighbour_species = ...
            [get_infected_species_incl_dead(cell_grid( mod(cell_i-2,cells_wide)+1, mod(cell_j,cells_long)+1 ),num_species), ...
            get_infected_species_incl_dead(cell_grid( mod(cell_i-2,cells_wide)+1, cell_j),num_species), ...
            get_infected_species_incl_dead(cell_grid( cell_i, mod(cell_j-2,cells_long)+1),num_species), ...
            get_infected_species_incl_dead(cell_grid( cell_i, mod(cell_j,cells_long)+1 ),num_species), ...
            get_infected_species_incl_dead(cell_grid( mod(cell_i,cells_wide)+1, mod(cell_j,cells_long)+1 ),num_species), ...
            get_infected_species_incl_dead(cell_grid( mod(cell_i,cells_wide)+1, cell_j ) ,num_species)];

    end

    tot_same_species_neighbours = sum((neighbour_species == species_this_cell));
    neighbour_clustering_all_cells(cell_ind) = tot_same_species_neighbours;

end

%crop neighbour_purity_all_cells just to the entries corresponding to
%fluorescent cells
neighbour_clustering_all_cells = neighbour_clustering_all_cells((is_inf_or_dead>0));


%for each possible number of same-species neighbours, report the proportion
%of the analysed cells that had that number of same-species neighbours
number_same_neighbours_prob = zeros(1,7);
if num_inf_or_dead>0
    for num_same_neighbours = 0:6
        number_same_neighbours_prob(num_same_neighbours+1) = ...
            sum((neighbour_clustering_all_cells==num_same_neighbours))/num_inf_or_dead;
    end
end


%NOTE: if analyses no infected cells, returns vector of zeros
