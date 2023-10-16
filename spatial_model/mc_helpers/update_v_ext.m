%update the global extracellular virus

production_vector = 0*v_ext;

for i=1:num_species
    production_vector(i) = sum(sum( (cell_grid == infected_donor_ids(i)) + (cell_grid == infected_recip_ids(i))))/numel(cell_grid);
end
               
v_ext = v_ext + dt*((1-na_inhib_prod) * p * production_vector - c*v_ext);