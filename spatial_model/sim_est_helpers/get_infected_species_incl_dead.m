function species = get_infected_species_incl_dead(id, num_species)

    eclipse_ids = 2:(2+num_species-1);
    infected_recip_ids = (eclipse_ids(end) + 1):(eclipse_ids(end) + 1 + num_species - 1);
    infected_recip_dead_ids = (infected_recip_ids(end)+1):(infected_recip_ids(end)+1 + num_species-1);
    infected_donor_ids = (infected_recip_dead_ids(end) + 1):(infected_recip_dead_ids(end) + 1 + num_species - 1);
    infected_donor_dead_ids = (infected_donor_ids(end)+1):(infected_donor_ids(end)+1+num_species-1);

    [~,s1] = ismember(id, infected_recip_ids);
    [~,s2] = ismember(id, infected_donor_ids);

    [~,s3] = ismember(id, infected_recip_dead_ids);
    [~,s4] = ismember(id, infected_donor_dead_ids);


    species = s1+s2+s3+s4;

end
