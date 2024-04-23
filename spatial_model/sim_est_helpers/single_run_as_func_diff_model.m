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
%  - Extracellular virus: PDE
%  - Intracellular virus: neighbour-based rule


function [sim_out] = single_run_as_func_diff_model(out_times, params)


    addpath mc_helpers

    %load params
    biology_parameters_w_diff;
    control_parameters_w_diff;
    
    
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

        virus_diff = params.virus_diff;
    end
    
    
    %set up output array
    if exist('out_times', 'var')
        
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

    %set up mesh
    mesh_parameters;
    make_diff_matrix;
    

    %initialise
    if random_initial
        initialise_cells_random_state;
    else
        initialise_cells_custom;
    end


    %specify moi here!
    moi = sum(sum((cell_grid~=target_id)))/total_cells;


    %count number of initially infected cells
    num_init_infected = 0;
    for species = 1:num_species
        num_init_infected = num_init_infected + ...
            sum(sum((cell_grid == infected_donor_ids(species)) +...
                    (cell_grid == infected_recip_ids(species))));
    end

    

    eclipse_wait_times = 0*cell_grid;
    v_ext=zeros(cells_wide, cells_long, num_species);
    
    cumulative_CC_infections = 0;
    cumulative_CF_infections = 0;


    %grid plot
    if vis_grid
        grid_fig = figure;
    end


    %main loop
    t=0;
    stopping_criterion_reached=0;
    num_output_points = length(out_times);
    while t<=final_time && ~stopping_criterion_reached

        %make next grid
        make_next_grid_w_diff;

        %update v_ext
        update_v_ext_w_diff;
        
        
        %update the grid
        cell_grid = next_grid;
        
        %update eclipse wait times
        eclipse_wait_times = eclipse_wait_times - dt;

        

        %visualise
        if vis_grid && ~mod(round(t/dt), frame_interval)
            pause(0.5)
            figure(grid_fig)
            plot_grid;
        end


        %optional output at specified times
        if exist('out_times', 'var')
            if curr_output_ind <= length(out_times) &&...
                    abs(t - out_times(curr_output_ind))<0.5*dt
        
                %save proportion of target cells
                net_T(curr_output_ind) = sum(sum((cell_grid == target_id)))/total_cells;
                
                %also virus
                net_V(curr_output_ind) = cell_area*sum(sum(sum(v_ext)));
                
                
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
                
                
                curr_output_ind = curr_output_ind+1;
            end
        end


        %optional clustering data output
        if save_clustering_data
            if curr_clustering_output_ind <= length(clustering_calc_times) &&...
                    abs(t - clustering_calc_times(curr_clustering_output_ind))<0.5*dt

                [mean_clustering, ERR_FLAG] = compute_mean_neighbour_purity(abs(cell_grid),num_species,num_samples_for_clustering);
                clustering_data(curr_clustering_output_ind) = mean_clustering;
                tot_clustering_errs = tot_clustering_errs + ERR_FLAG;
                    
                curr_clustering_output_ind = curr_clustering_output_ind + 1;
            end
        end
        
        
        %check if stopping criterion reached
        if stop_when_all_infected
            if sum(sum((cell_grid == target_id)))==0
                stopping_criterion_reached = 1;
                num_output_points=curr_output_ind-1;
                
                if (num_output_points < length(out_times))
                    disp('WARNING: stopped sim (all cells infected)')
                end
                
            end
        end

        if stop_at_infected_peak

            MIN_SIM_TIME = 10;
            PROP_OF_PEAK = 0.5;

            if t>MIN_SIM_TIME && net_I(curr_output_ind-1)<PROP_OF_PEAK*max(net_I)
                stopping_criterion_reached = 1;
                num_output_points=curr_output_ind-1;
                
                if (num_output_points < length(out_times))
                    disp('WARNING: stopped sim (passed infected peak)')
                end
                
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
            net_V(curr_output_ind) = cell_area*sum(sum(sum(v_ext)));


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


            curr_output_ind = curr_output_ind+1;
        end
    end


    %optional clustering data output (final point)
    if save_clustering_data
        if curr_clustering_output_ind <= length(clustering_calc_times) &&...
                abs(t - clustering_calc_times(curr_clustering_output_ind))<0.5*dt

            [mean_clustering, ERR_FLAG] = compute_mean_neighbour_purity(abs(cell_grid),num_species,num_samples_for_clustering);
            clustering_data(curr_clustering_output_ind) = mean_clustering;
            tot_clustering_errs = tot_clustering_errs + ERR_FLAG;

            curr_clustering_output_ind = curr_clustering_output_ind + 1;
        end
    end




    %output information
    sim_out.prop_infected = prop_infected(1:num_output_points);
    sim_out.net_T = net_T(1:num_output_points);
    sim_out.net_E = net_E(1:num_output_points);
    sim_out.net_I = net_I(1:num_output_points);
    sim_out.net_V = net_V(1:num_output_points)/max(net_V);
    sim_out.net_CC = net_CC(1:num_output_points);
    sim_out.net_CF = net_CF(1:num_output_points);
    sim_out.final_grid_state = abs(cell_grid);
    sim_out.num_output_points = num_output_points;

    if save_clustering_data
        sim_out.clustering_data = clustering_data;
        sim_out.tot_clustering_errs = tot_clustering_errs;
    end
    
end
    
    
    
    
    
    
    
    
    
    
    