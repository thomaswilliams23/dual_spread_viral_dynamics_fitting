

%Makes Figure 5a from the manuscript



%set font of plots
set(0,'defaulttextinterpreter','latex')


%% setup

addpath helper_funcs/

%folder structure
source_folder = '../SIM_DATA_from_ms/fig5a_example_sims';
template_1 = '/Pcc_';
template_2 = '_virus_diff_';


%sweep variables
Pcc_sweep = [0.1, 0.5, 0.9];
Pcc_names = {'01', '05', '09'};
diff_names = {'01', '1', '10', '100', 'inf'};


%% loop
%loop Pcc
for Pcc_ind = 1:length(Pcc_sweep)

    %loop diffusion
    for diff_ind = 1:length(diff_names)

        %load data
        load(strcat(source_folder, template_1, num2str(Pcc_sweep(Pcc_ind)), ...
            template_2, diff_names{diff_ind}, '/all_sim_data.mat'));

        this_sim = all_sim_data{1};

        %name variables
        cell_grid = this_sim.final_grid_state;
        num_species = round(0.01*numel(this_sim.final_grid_state));


        %plot
        plot_final_grid;


        %export
        exportgraphics(gcf, strcat('Pcc_',Pcc_names{Pcc_ind}, '_diff_', diff_names{diff_ind}, '.png'), 'Resolution', 300)
        pause(1)
        close

    end
end