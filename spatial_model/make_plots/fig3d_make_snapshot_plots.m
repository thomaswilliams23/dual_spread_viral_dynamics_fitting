

%Makes Figure 3d from the manuscript



%add folders
addpath helper_funcs


%% simulation data specs 

%data source
data_source = '../SIM_DATA_from_ms/fig3d_example_sims/';

% specify a data structure and times and get snapshots of the cell grid at those times
sim_to_load = strcat(data_source, 'low_CC_sim.mat');
% sim_to_load = strcat(data_source, 'equal_CC_sim.mat');
% sim_to_load = strcat(data_source, 'high_CC_sim.mat');

%time points to extract (max is 50)
dt = 0.1;
times_to_extract = round([0.1, 4, 8, 12, 16]/dt);

%size
cells_wide = 50;
cells_long = 50;

%species
num_species = 25; %corresponds to one colour for each initially infected cell
species_order = randperm(num_species); %permute the colour order

%cell IDs
target_id = 1;
eclipse_ids = 2:(2+num_species-1);
infected_recip_ids = (eclipse_ids(end) + 1):(eclipse_ids(end) + 1 + num_species - 1);
infected_recip_dead_ids = (infected_recip_ids(end)+1):(infected_recip_ids(end)+1 + num_species-1);
infected_donor_ids = (infected_recip_dead_ids(end) + 1):(infected_recip_dead_ids(end) + 1 + num_species - 1);
infected_donor_dead_ids = (infected_donor_ids(end)+1):(infected_donor_ids(end)+1+num_species-1);
num_cell_types = infected_donor_dead_ids(end);


%% plot

%load simulation data
load(sim_to_load)

%loop over time points
for time = times_to_extract

    %extract cell grid
    cell_grid = abs(sim_out.grid_frames(:,:,time));

    %plot
    plot_grid;
end

