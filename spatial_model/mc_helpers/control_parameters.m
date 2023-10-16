

%% model options
moi = 0.01;

periodic_BCs = 1;

latent_stages = 3;

random_initial = 1;
initialise_in_window = 0;

show_histogram_of_inf_neighbourhoods = 0;

window_width = 200;

num_species = 1;


%% time options
dt = 0.01;
final_time = 24;

stop_when_all_infected = 0;
stop_at_infected_peak = 0;


%% plotting
vis_grid = 1;

save_clustering_data = 0;
clustering_calc_times = [];
num_samples_for_clustering = -1; %this is the flag for taking all samples

frame_interval = round(1/dt);


%% cell IDs

target_id = 1;
eclipse_ids = 2:(2+num_species-1);
infected_recip_ids = (eclipse_ids(end) + 1):(eclipse_ids(end) + 1 + num_species - 1);
infected_recip_dead_ids = (infected_recip_ids(end)+1):(infected_recip_ids(end)+1 + num_species-1);
infected_donor_ids = (infected_recip_dead_ids(end) + 1):(infected_recip_dead_ids(end) + 1 + num_species - 1);
infected_donor_dead_ids = (infected_donor_ids(end)+1):(infected_donor_ids(end)+1+num_species-1);
num_cell_types = infected_donor_dead_ids(end);


%ALL INFECTED IDs NEGATIVE FOR FASTER LOOKUP
infected_recip_ids = -infected_recip_ids;
infected_donor_ids = -infected_donor_ids;




