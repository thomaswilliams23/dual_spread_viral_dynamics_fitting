function default_prm = make_default_param_struct()

    %uses params in the parameter files to populate a param struct
    
    biology_parameters;
    control_parameters;
    
    
    %read in
    default_prm.moi = moi;

    default_prm.dt = dt;
    
    default_prm.total_cells = total_cells;
    
    default_prm.random_initial = random_initial;
    default_prm.initialise_in_window = initialise_in_window;

    default_prm.show_histogram_of_inf_neighbourhoods = show_histogram_of_inf_neighbourhoods;
    
    default_prm.window_width = window_width;
    
    default_prm.stop_when_all_infected = stop_when_all_infected;
    default_prm.stop_at_infected_peak = stop_at_infected_peak;
    
    default_prm.final_time = final_time;
    
    default_prm.latent_stages = latent_stages;
    
    default_prm.num_species = num_species;
    
    default_prm.T_tot = T_tot;
    
    default_prm.vis_grid = vis_grid;


    default_prm.save_clustering_data = save_clustering_data;
    default_prm.clustering_calc_times = clustering_calc_times;
    default_prm.num_samples_for_clustering = num_samples_for_clustering;
    
    
    default_prm.frame_interval = frame_interval;
    
    
    default_prm.delta_param = delta_param;

    default_prm.beta_param = beta_param;
    default_prm.alpha_param = alpha_param;
    
    default_prm.gamma_param = gamma_param;

    default_prm.p = p;

    default_prm.c = c;

    default_prm.na_inhib_prod = na_inhib_prod;
    
end