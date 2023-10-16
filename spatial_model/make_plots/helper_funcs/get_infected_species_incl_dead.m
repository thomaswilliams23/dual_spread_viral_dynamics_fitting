function species = get_infected_species_incl_dead(id, true_num_species)
% finds the species number of a cell given its cell id 
% if only one species, returns 1 for fluorescent and 0 for not fluorescent

    %all indices
    eclipse_ids = 2:(2+true_num_species-1);
    infected_recip_ids = (eclipse_ids(end) + 1):(eclipse_ids(end) + 1 + true_num_species - 1);
    infected_recip_dead_ids = (infected_recip_ids(end)+1):(infected_recip_ids(end)+1 + true_num_species-1);
    infected_donor_ids = (infected_recip_dead_ids(end) + 1):(infected_recip_dead_ids(end) + 1 + true_num_species - 1);
    infected_donor_dead_ids = (infected_donor_ids(end)+1):(infected_donor_ids(end)+1+true_num_species-1);

    %check membership
    [~,s1] = ismember(id, infected_recip_ids);
    [~,s2] = ismember(id, infected_donor_ids);

    [~,s3] = ismember(id, infected_recip_dead_ids);
    [~,s4] = ismember(id, infected_donor_dead_ids);

    %compute species number
    species = s1+s2+s3+s4;
end
