%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                     DUAL SPREAD VIRUS SYSTEM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Implements a CA virus model with two modes of viral spread
%
%
% Model structure:
%
%  - Cells: Target, Eclipse, Infected (recipient), Infected (recipient, 
%    dead), Infected (donor), Infected (donor, dead). 
%    Arranged in hexagonal grid
%
%  - Extracellular virus: ODE
%  - Intracellular virus: neighbour-based rule

function [sim_out] = single_run_as_func(out_times, params)


    addpath mc_helpers

    %load params
    biology_parameters;
    control_parameters;
    
    
    %argument parameters (overrides built-in)
    if exist('params','var')
        
        %basic
        moi = params.moi;
        dt=params.dt;
        total_cells = params.total_cells;
        
        %termination
        final_time = params.final_time;
        stop_when_all_infected = params.stop_when_all_infected;
        stop_at_infected_peak = params.stop_at_infected_peak;
        

        %visualisation
        vis_grid = params.vis_grid;
        frame_interval = params.frame_interval;
        latent_stages = params.latent_stages;
        num_species = params.num_species;


        %clustering
        save_clustering_data = params.save_clustering_data;
        clustering_calc_times = params.clustering_calc_times;
        num_samples_for_clustering = params.num_samples_for_clustering;
        
        
        %biological
        beta_param = params.beta_param;
        alpha_param = params.alpha_param;
        gamma_param = params.gamma_param;
        delta_param = params.delta_param;
        p = params.p;
        c = params.c;
    end
    
    
    %set up output arrays
    if exist('out_times', 'var')

        num_output_points = length(out_times);
        
        prop_infected = 0*out_times;

        net_T = 0*out_times;
        net_E = 0*out_times;
        net_I = 0*out_times;
        net_V = 0*out_times;
        
        net_CC = 0*out_times;
        net_CF = 0*out_times;

        if save_clustering_data
            clustering_data = 0*clustering_calc_times;
            curr_clustering_output_ind = 1;
            tot_clustering_errs = 0;
        end
        
        curr_output_ind = 1;
    end


    %set up cell IDs
    target_id = 1;
    eclipse_ids = 2:(2+num_species-1);
    infected_recip_ids = (eclipse_ids(end) + 1):(eclipse_ids(end) + 1 + num_species - 1);
    infected_recip_dead_ids = (infected_recip_ids(end)+1):(infected_recip_ids(end)+1 + num_species-1);
    infected_donor_ids = (infected_recip_dead_ids(end) + 1):(infected_recip_dead_ids(end) + 1 + num_species - 1);
    infected_donor_dead_ids = (infected_donor_ids(end)+1):(infected_donor_ids(end)+1+num_species-1);
    num_cell_types = infected_donor_dead_ids(end);    
    
    %infected IDs negative for faster lookup
    infected_recip_ids = -infected_recip_ids;
    infected_donor_ids = -infected_donor_ids;

    %set up species order
    species_order = randperm(num_species);
    

    %initialise
    initialise_cells_random_state;



    %count number of initially infected cells
    num_init_infected = 0;
    for species = 1:num_species
        num_init_infected = num_init_infected + ...
            sum(sum((cell_grid == infected_donor_ids(species)) +...
                    (cell_grid == infected_recip_ids(species))));
    end

    
    %initalise some more structures
    eclipse_wait_times = 0*cell_grid;
    v_ext=zeros(1, num_species);
    
    cumulative_CC_infections = 0;
    cumulative_CF_infections = 0;


    %make plot figure
    if vis_grid
        grid_fig = figure;
    end


    %main loop
    t=0;
    while t<=final_time

        %make next grid
        make_next_grid;

        %update v_ext
        update_v_ext;
        
        
        %update the grid
        cell_grid = next_grid;
        
        %update eclipse wait times
        eclipse_wait_times = eclipse_wait_times - dt;


        %visualise
        if vis_grid && ~mod(round(t/dt), frame_interval)
            pause(0.5)
            plot_grid;
        end



        %optional output at specified times
        if exist('out_times', 'var')
            if curr_output_ind <= length(out_times) &&...
                    abs(t - out_times(curr_output_ind))<0.5*dt
        
                %save proportion of target cells
                net_T(curr_output_ind) = sum(sum((cell_grid == target_id)))/total_cells;
                
                %also virus
                net_V(curr_output_ind) = sum(v_ext);
                
                
                for species = 1:num_species
                    %save the proportion of recipient infected cells (incl. dead)
                    prop_infected(curr_output_ind) = prop_infected(curr_output_ind) + ...
                        (sum(sum((cell_grid==infected_recip_ids(species)) + ...
                                 (cell_grid==infected_recip_dead_ids(species)))));
                        
                    %save the proportion of target and infected cells
                    net_E(curr_output_ind) = net_E(curr_output_ind) + sum(sum((cell_grid == eclipse_ids(species))))/total_cells;
                    net_I(curr_output_ind) = net_I(curr_output_ind) + sum(sum((cell_grid == infected_recip_ids(species)) + ...
                        (cell_grid == infected_donor_ids(species))))/total_cells;
                end

                %normalise
                prop_infected(curr_output_ind) = prop_infected(curr_output_ind)/(numel(cell_grid) - num_init_infected);
                
                %save cumulative CC and CF infections
                net_CC(curr_output_ind) = cumulative_CC_infections/(numel(cell_grid) - num_init_infected);
                net_CF(curr_output_ind) = cumulative_CF_infections/(numel(cell_grid) - num_init_infected);
                


                %update the output index
                curr_output_ind = curr_output_ind+1;
            end
        end


        %optional clustering data output
        if save_clustering_data
            if curr_clustering_output_ind <= length(clustering_calc_times) &&...
                    abs(t - clustering_calc_times(curr_clustering_output_ind))<0.5*dt

                %compute clustering metric
                [mean_clustering, ERR_FLAG] = compute_mean_neighbour_clustering(abs(cell_grid),num_species,num_samples_for_clustering);
                clustering_data(curr_clustering_output_ind) = mean_clustering;
                tot_clustering_errs = tot_clustering_errs + ERR_FLAG;
                    
                %update clustering metric output index
                curr_clustering_output_ind = curr_clustering_output_ind + 1;
            end
        end        
        

        t=t+dt;
    end



    %optional output at specified times (final time point)
    if exist('out_times', 'var')
        if curr_output_ind <= length(out_times) &&...
                abs(t - out_times(curr_output_ind))<0.5*dt
    
            %save proportion of target cells
            net_T(curr_output_ind) = sum(sum((cell_grid == target_id)))/total_cells;
            
            %also virus
            net_V(curr_output_ind) = sum(v_ext);
            
            
            for species = 1:num_species
                %save the proportion of recipient infected cells (incl. dead)
                prop_infected(curr_output_ind) = prop_infected(curr_output_ind) + ...
                    (sum(sum((cell_grid==infected_recip_ids(species)) + ...
                             (cell_grid==infected_recip_dead_ids(species)))));
                    
                %save the proportion of target and infected cells
                net_E(curr_output_ind) = net_E(curr_output_ind) + sum(sum((cell_grid == eclipse_ids(species))))/total_cells;
                net_I(curr_output_ind) = net_I(curr_output_ind) + sum(sum((cell_grid == infected_recip_ids(species)) + ...
                    (cell_grid == infected_donor_ids(species))))/total_cells;
            end

            prop_infected(curr_output_ind) = prop_infected(curr_output_ind)/(numel(cell_grid) - num_init_infected);
            
            
            %save cumulative CC and CF infections
            net_CC(curr_output_ind) = cumulative_CC_infections/(numel(cell_grid) - num_init_infected);
            net_CF(curr_output_ind) = cumulative_CF_infections/(numel(cell_grid) - num_init_infected);
        end
    end



    %optional clustering data output (final time point)
    if save_clustering_data
        if curr_clustering_output_ind <= length(clustering_calc_times) &&...
                abs(t - clustering_calc_times(curr_clustering_output_ind))<0.5*dt

            %compute clustering metric
            [mean_clustering, ERR_FLAG] = compute_mean_neighbour_clustering(abs(cell_grid),num_species,num_samples_for_clustering);
            clustering_data(curr_clustering_output_ind) = mean_clustering;
            tot_clustering_errs = tot_clustering_errs + ERR_FLAG;
        end
    end
  

    
    %output information in structure
    sim_out.prop_infected = prop_infected;
    sim_out.net_T = net_T;
    sim_out.net_E = net_E;
    sim_out.net_I = net_I;
    sim_out.net_V = net_V/max(net_V);
    sim_out.net_CC = net_CC;
    sim_out.net_CF = net_CF;
    sim_out.final_grid_state = abs(cell_grid);
    sim_out.num_output_points = num_output_points;

    if save_clustering_data
        sim_out.clustering_data = clustering_data;
        sim_out.tot_clustering_errs = tot_clustering_errs;
    end

end
    
    
    
    
    
    
    
    
    
    
    