

%for reproducibility
rng(0)


%colour scheme
cell_colours = jet(num_species);
cell_colours = cell_colours(randperm(num_species),:);

%any remaining target cells black
target_colour = [0,0,0];



%% preprocess the grid

get_spec_num = @(id) get_infected_species_incl_dead(id, num_species);
cell_grid_species = arrayfun(get_spec_num, abs(cell_grid));


%% make the figure
figure

%loop over species (assumes everything is infected)
for spec = 1:num_species

    %find cells with this cell type
    [x_vals, y_vals] = find(cell_grid_species==spec);

    %offset for hexagonal structure
    y_vals = y_vals + 0.5*(mod(x_vals,2));

    %if any cells with this cell type
    if ~isempty(x_vals)

        %plot circular cells in the correct loaction
        circles(x_vals, y_vals, 0.5*ones(size(x_vals)), 'FaceColor', cell_colours(spec,:), ...
                'EdgeColor', 'none');
        hold on
    end

end

%% now plot any target cells

target_id = 0;

%plot target cells
[x_vals, y_vals] = find(cell_grid_species==target_id);
y_vals = y_vals + 0.5*(mod(x_vals,2));
if ~isempty(x_vals)
    circles(x_vals, y_vals, 0.5*ones(size(x_vals)), 'FaceColor', target_colour, ...
            'EdgeColor', 'none');
end



%% formatting

%take off the axis labels
xticks([])
yticks([])

xlim([0,size(cell_grid,1)+0.5])
ylim([0,size(cell_grid,2)+1])


axis off

x_size = size(cell_grid,1);
y_size = size(cell_grid,2);
ax = gca;
ax.PlotBoxAspectRatio = [1,y_size/x_size,1];

fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = [0.1, 0.1, 0.8, 0.8];

hold off