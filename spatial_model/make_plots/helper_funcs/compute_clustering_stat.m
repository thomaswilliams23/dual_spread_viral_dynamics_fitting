function clustering_stat = compute_clustering_stat(cell_grid)
%computes kappa(t) from a given state of the cell grid


%dimensions of cell grid
cells_wide = size(cell_grid,1);
cells_long = size(cell_grid,2);

%number of species and neighbours [FIXED]
num_species = 1;
num_neighbours = 6;

%initialise
clustering_metric_all_cells = zeros(cells_wide*cells_long, 1);
num_fluor = 0;


%loop over cells
for cell_i = 1:cells_wide
    for cell_j = 1:cells_long

        %retrieves species number of cell 
        %(here, since only one species, either 1 for fluorescent or 0 for not)
        this_cell_inf = get_infected_species_incl_dead(cell_grid(cell_i, cell_j), num_species);

        % compute neighbours if this cell fluorescent, otherwise skip
        if this_cell_inf
            num_fluor = num_fluor + 1;
        else
            continue
        end


        %loop over neighbours, check if they are fluorescent
        if mod(cell_i,2) %odd row, left-aligned

            fluor_neighbours = ...
                [get_infected_species_incl_dead(cell_grid( mod(cell_i-2,cells_wide)+1, mod(cell_j-2,cells_long)+1 ),num_species), ...
                get_infected_species_incl_dead(cell_grid( mod(cell_i-2,cells_wide)+1, cell_j ),num_species), ...
                get_infected_species_incl_dead(cell_grid( cell_i, mod(cell_j-2,cells_long)+1 ),num_species), ...
                get_infected_species_incl_dead(cell_grid( cell_i, mod(cell_j,cells_long)+1 ),num_species), ...
                get_infected_species_incl_dead(cell_grid( mod(cell_i,cells_wide)+1, mod(cell_j-2,cells_long)+1 ),num_species), ...
                get_infected_species_incl_dead(cell_grid( mod(cell_i,cells_wide)+1, cell_j),num_species)];

        else %even row, right-aligned

            fluor_neighbours = ...
                [get_infected_species_incl_dead(cell_grid( mod(cell_i-2,cells_wide)+1, mod(cell_j,cells_long)+1 ),num_species), ...
                get_infected_species_incl_dead(cell_grid( mod(cell_i-2,cells_wide)+1, cell_j),num_species), ...
                get_infected_species_incl_dead(cell_grid( cell_i, mod(cell_j-2,cells_long)+1),num_species), ...
                get_infected_species_incl_dead(cell_grid( cell_i, mod(cell_j,cells_long)+1 ),num_species), ...
                get_infected_species_incl_dead(cell_grid( mod(cell_i,cells_wide)+1, mod(cell_j,cells_long)+1 ),num_species), ...
                get_infected_species_incl_dead(cell_grid( mod(cell_i,cells_wide)+1, cell_j ) ,num_species)];

        end

        %tally total fluorescent neighbours
        tot_fluor_neighbours = sum(fluor_neighbours);

        %save this into an array
        clustering_metric_all_cells(cell_i + cells_wide*(cell_j-1)) = tot_fluor_neighbours;

    end
end

%loop over each possible number of neighbours, compute prob of having that
%number of fluorescent neighbours
num_fluor_neighbours_prob = zeros(1,num_neighbours+1);
for num_fluor_neighbours = 0:num_neighbours
    num_fluor_neighbours_prob(num_fluor_neighbours+1) = ...
        sum((clustering_metric_all_cells==num_fluor_neighbours))/num_fluor;
end

%now compute mean
clustering_stat = sum((0:num_neighbours).*num_fluor_neighbours_prob);

end