%makes the initial cell grid, choosing initially infected cells at random
cell_grid = target_id * ones(cells_long, cells_wide);


for inf_cell = 1:round(moi*total_cells)
    x_c=randi(cells_long);
    y_c=randi(cells_wide);

    %check not already infected
    while ismember(cell_grid(x_c,y_c),infected_donor_ids)
        x_c=randi(cells_long);
        y_c=randi(cells_wide);
    end

    cell_grid(x_c,y_c)=infected_donor_ids(species_order(mod(inf_cell-1,num_species)+1));
end

