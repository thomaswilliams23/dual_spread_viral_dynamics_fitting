next_grid = cell_grid;

%loop over cells
for i=1:cells_wide
    for j=1:cells_long

        %check cell type:
        if cell_grid(i,j)==target_id  %cell is TARGET

            %compute probability of being infected on next time step

            tot_neighbours = 6;
            
            %NOTE: MAKE ALL INFECTED SPECIES IDS NEGATIVE FOR SPEED
            if mod(i,2) %odd row, left-aligned
                tot_inf_neighbours = (cell_grid( mod(i-2,cells_wide)+1, mod(j-2,cells_long)+1 ) <0) + ... %recip
                     (cell_grid( mod(i-2,cells_wide)+1, j ) <0) + ...
                     (cell_grid( i, mod(j-2,cells_long)+1 ) <0) + ...
                     (cell_grid( i, mod(j,cells_long)+1 ) <0) + ...
                     (cell_grid( mod(i,cells_wide)+1, mod(j-2,cells_long)+1 ) <0) + ...
                     (cell_grid( mod(i,cells_wide)+1, j) <0);

            else %even row, right-aligned

                tot_inf_neighbours = (cell_grid( mod(i-2,cells_wide)+1, mod(j,cells_long)+1 ) <0) + ... %recip
                     (cell_grid( mod(i-2,cells_wide)+1, j) <0) + ...
                     (cell_grid( i, mod(j-2,cells_long)+1) <0) + ...
                     (cell_grid( i, mod(j,cells_long)+1 ) <0) + ...
                     (cell_grid( mod(i,cells_wide)+1, mod(j,cells_long)+1 ) <0) + ...
                     (cell_grid( mod(i,cells_wide)+1, j ) <0);

            end

             
            prob_t_to_e = 1 - exp(-(beta_param*sum(v_ext) + alpha_param*(tot_inf_neighbours/tot_neighbours))*dt);
            
            
            %decide mode of transmission
            prob_CC_unscaled = 1-exp(-alpha_param*(tot_inf_neighbours/tot_neighbours)*dt);
            prob_CF_unscaled = 1-exp(-beta_param*sum(v_ext)*dt);
            
            prob_CC_inf = prob_t_to_e * (prob_CC_unscaled/(prob_CC_unscaled + prob_CF_unscaled));
            prob_CF_inf = prob_t_to_e * (prob_CF_unscaled/(prob_CC_unscaled + prob_CF_unscaled));
            
            
            
            rnd_draw = rand;
            
            %CELL-CELL INFECTION
            if (rnd_draw <= prob_CC_inf)

                inf_neighbours = zeros(1,num_species);

                if mod(i,2) %odd row, left-aligned
                    
                    bl_species = get_infected_species(abs(cell_grid( mod(i-2,cells_wide)+1, mod(j-2,cells_long)+1 )), num_species);
                    if bl_species>0
                        inf_neighbours(bl_species) = inf_neighbours(bl_species) + 1;
                    end
                    
                    l_species = get_infected_species(abs(cell_grid( mod(i-2,cells_wide)+1, j )), num_species);
                    if l_species>0
                        inf_neighbours(l_species) = inf_neighbours(l_species) + 1;
                    end
                    
                    b_species = get_infected_species(abs(cell_grid( i, mod(j-2,cells_long)+1 )), num_species);
                    if b_species>0
                        inf_neighbours(b_species) = inf_neighbours(b_species) + 1;
                    end
                    
                    t_species = get_infected_species(abs(cell_grid( i, mod(j,cells_long)+1 )), num_species);
                    if t_species>0
                        inf_neighbours(t_species) = inf_neighbours(t_species) + 1;
                    end
                    
                    tr_species = get_infected_species(abs(cell_grid( mod(i,cells_wide)+1, mod(j-2,cells_long)+1 )), num_species);
                    if tr_species>0
                        inf_neighbours(tr_species) = inf_neighbours(tr_species) + 1;
                    end
                    
                    r_species = get_infected_species(abs(cell_grid(  mod(i,cells_wide)+1, j )), num_species);
                    if r_species>0
                        inf_neighbours(r_species) = inf_neighbours(r_species) + 1;
                    end
                    
                else       %even row, right-aligned
                    
                    bl_species = get_infected_species(abs(cell_grid( mod(i-2,cells_wide)+1, mod(j,cells_long)+1 )), num_species);
                    if bl_species>0
                        inf_neighbours(bl_species) = inf_neighbours(bl_species) + 1;
                    end
                    
                    l_species = get_infected_species(abs(cell_grid( mod(i-2,cells_wide)+1, j)), num_species);
                    if l_species>0
                        inf_neighbours(l_species) = inf_neighbours(l_species) + 1;
                    end
                    
                    b_species = get_infected_species(abs(cell_grid( i, mod(j-2,cells_long)+1  )), num_species);
                    if b_species>0
                        inf_neighbours(b_species) = inf_neighbours(b_species) + 1;
                    end
                    
                    t_species = get_infected_species(abs(cell_grid( i, mod(j,cells_long)+1 )), num_species);
                    if t_species>0
                        inf_neighbours(t_species) = inf_neighbours(t_species) + 1;
                    end
                    
                    tr_species = get_infected_species(abs(cell_grid( mod(i,cells_wide)+1, mod(j,cells_long)+1 )), num_species);
                    if tr_species>0
                        inf_neighbours(tr_species) = inf_neighbours(tr_species) + 1;
                    end
                    
                    r_species = get_infected_species(abs(cell_grid( mod(i,cells_wide)+1, j )), num_species);
                    if r_species>0
                        inf_neighbours(r_species) = inf_neighbours(r_species) + 1;
                    end
                    
                end

                
                %quantify contributions from each species
                cumulative_unscaled_species_probs = zeros(1, num_species);
                cumulative_unscaled_prob = 0;
                for species = 1:num_species
                    cumulative_unscaled_prob = cumulative_unscaled_prob + ...
                        1 - exp(-(alpha_param*inf_neighbours(species)/tot_neighbours)*dt);
                    cumulative_unscaled_species_probs(species) = cumulative_unscaled_prob;
                end
                species_probs = prob_CC_inf * (cumulative_unscaled_species_probs/cumulative_unscaled_species_probs(end));


                for species = 1:num_species
                    if rnd_draw<=species_probs(species)
                        next_grid(i,j) = eclipse_ids(species);
                        wait_time = gamrnd(latent_stages, 1/(latent_stages*gamma_param));
                        eclipse_wait_times(i,j) = wait_time;
                        break
                    end
                end
                
                
                %keep track of CC infections
                if exist('out_times', 'var')
                    cumulative_CC_infections = cumulative_CC_infections + 1;
                end
                    
            
                
            %CELL-FREE INFECTION    
            elseif ((rnd_draw<=prob_t_to_e) && (rnd_draw>prob_CC_inf))
                    
                %quantify contributions from each species
                cumulative_unscaled_species_probs = zeros(1, num_species);
                cumulative_unscaled_prob = 0;
                for species = 1:num_species
                    cumulative_unscaled_prob = cumulative_unscaled_prob + ...
                        1 - exp(-beta_param*v_ext(species)*dt);
                    cumulative_unscaled_species_probs(species) = cumulative_unscaled_prob;
                end
                species_probs = prob_CC_inf + prob_CF_inf * (cumulative_unscaled_species_probs/cumulative_unscaled_species_probs(end));    
                    
                for species = 1:num_species
                    if rnd_draw<=species_probs(species)
                        next_grid(i,j) = eclipse_ids(species);
                        wait_time = gamrnd(latent_stages, 1/(latent_stages*gamma_param));
                        eclipse_wait_times(i,j) = wait_time;
                        break
                    end
                end
                
                %keep track of CF infections
                if exist('out_times', 'var')
                    cumulative_CF_infections = cumulative_CF_infections + 1;
                end
                
            end
            

            
            
            
            
        elseif ismember(cell_grid(i,j),eclipse_ids)    %cell is ECLIPSE
            
            if eclipse_wait_times(i,j) < 0
                [~,species] = ismember(cell_grid(i,j), eclipse_ids);
                next_grid(i,j) = infected_recip_ids(species);
            end
            
            

        elseif ismember(cell_grid(i,j),infected_recip_ids) %cell is INFECTED (recipient)

            %compute prob of cell dying in next time step
            prob_i_to_d = 1-exp(-delta_param*dt);

            %draw random number, decide next state of cell
            if rand<prob_i_to_d %becomes dead
                [~,species] = ismember(cell_grid(i,j), infected_recip_ids);
                next_grid(i,j) = infected_recip_dead_ids(species);
            end
        
        
        
        elseif ismember(cell_grid(i,j),infected_donor_ids) %cell is INFECTED (donor)

            %compute prob of cell dying in next time step
            prob_i_to_d = 1-exp(-delta_param*dt);

            %draw random number, decide next state of cell
            if rand<prob_i_to_d %becomes dead
                [~,species] = ismember(cell_grid(i,j), infected_donor_ids);
                next_grid(i,j) = infected_donor_dead_ids(species);
            end
            
        end  %dead cells do nothing
        
    end
end