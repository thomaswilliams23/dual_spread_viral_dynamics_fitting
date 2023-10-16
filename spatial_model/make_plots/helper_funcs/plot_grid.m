
%colour scheme
target_colour = target_id;
species_colours = abs([eclipse_ids; infected_recip_ids; infected_donor_ids]);
dead_colour = [infected_recip_dead_ids, infected_donor_dead_ids];


cell_colours = parula(num_cell_types);

infected_cell_colours = jet(num_species);
infected_cell_colours = infected_cell_colours(species_order,:);


cell_colours(target_colour,:) = [0.85, 0.85, 0.85];    %target cells are grey
for i = 1:length(dead_colour)
    cell_colours(dead_colour(i),:) = [0,0,0];        %dead cells are black
end
for i = 1:size(species_colours,2) % infected cells all get the same colour per species
    for j = 1:size(species_colours,1)
        cell_colours(species_colours(j,i),:) = infected_cell_colours(i,:);
    end
end


%% plot
figure
hold on

%loop over cell types
for cell_ind = 1:num_cell_types

    %find cells with this cell type
    [x_vals, y_vals] = find(cell_grid==cell_ind);

    %offset for hexagonal structure
    y_vals = y_vals + 0.5*(~mod(x_vals,2));

    %if any cells with this cell type
    if ~isempty(x_vals)

        %plot circular cells in the correct loaction
        circles(x_vals, y_vals, 0.5*ones(size(x_vals)), 'FaceColor', cell_colours(cell_ind,:), ...
            'EdgeColor', 'none');
    end
end

%formatting
xlim([0,size(cell_grid,1)+0.5]);
ylim([0,size(cell_grid,2)+1]);

xticks([])
yticks([])

axis off

ax = gca;
ax.PlotBoxAspectRatio = [1,1,1];


