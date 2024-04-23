

%Makes Figure 5g from the manuscript



%set font of plots
set(0, 'DefaultTextInterpreter', 'latex')


%% setup

%folder structure
data_source = '../SIM_DATA_from_ms/fig5bcdef_data_for_kappa';
template_1 = '/Pcc_';
template_2 = '_virus_diff_';

upper_Pcc_name = '0.9';
lower_Pcc_name = '0.1';

diff_vals = [0.1, 10^(-0.5), 1, 10^0.5,  10,  10^(1.5), 100];
diff_names = {'01', '10^-05', '1', '10^05',  '10',  '10^15', '100'};

num_reps = 8;


%time vector
dt = 0.1;
final_time = 40;
num_time_points = round(final_time/dt);


%initialise
max_vert_distance = zeros(1, length(diff_vals));


%% loop to compute variation
for diff_ind = 1:length(diff_vals)


    %load data
    load(strcat(data_source, template_1, upper_Pcc_name, template_2, diff_names{diff_ind}, '/all_sim_data'));
    all_sim_data_upper = all_sim_data;

    load(strcat(data_source, template_1, lower_Pcc_name, template_2, diff_names{diff_ind}, '/all_sim_data'));
    all_sim_data_lower = all_sim_data;


    %form an average trajectory
    mean_upper_traj = zeros(1,num_time_points);
    mean_lower_traj = zeros(1,num_time_points);

    for rep = 1:num_reps
        mean_upper_traj = mean_upper_traj + all_sim_data_upper{rep}.kappa_data/num_reps;
        mean_lower_traj = mean_lower_traj + all_sim_data_lower{rep}.kappa_data/num_reps;
    end


    %now find vertical distance
    vert_distance_over_time = abs(mean_upper_traj-mean_lower_traj);

    %and report max
    max_vert_distance(diff_ind) = max(vert_distance_over_time);

end

%% add infinite case

%load data
load(strcat(data_source, template_1, upper_Pcc_name, template_2, 'Inf', '/all_sim_data'));
all_sim_data_upper = all_sim_data;

load(strcat(data_source, template_1, lower_Pcc_name, template_2, 'Inf', '/all_sim_data'));
all_sim_data_lower = all_sim_data;


%form an average trajectory
mean_upper_traj = zeros(1,num_time_points);
mean_lower_traj = zeros(1,num_time_points);

for rep = 1:num_reps
    mean_upper_traj = mean_upper_traj + all_sim_data_upper{rep}.kappa_data/num_reps;
    mean_lower_traj = mean_lower_traj + all_sim_data_lower{rep}.kappa_data/num_reps;
end


%now find vertical distance
vert_distance_over_time = abs(mean_upper_traj-mean_lower_traj);

%and report max
inf_vert_distance = max(vert_distance_over_time);



%% plot

figure
plt = semilogx(diff_vals, max_vert_distance, 'color', [0.8500 0.3250 0.0980]);
hold on
plot([0.1, 100], inf_vert_distance*[1,1], 'k--')
text(10^(-0.9), 0.55, '$D\rightarrow\infty$', 'Interpreter','latex', 'FontSize',9)

ylim([0, 0.7])

set(plt, 'linewidth', 1.5)
set(plt, 'marker', 'o')

xlabel('$D$');
ylabel('Max Vertical Distance')

box off

