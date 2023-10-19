



%%% RUNS ONE ITERATION OF THE ODE MODEL AND PLOTS TIME SERIES


%set up path
addpath matlab_helpers/



%% parameters

%proportion of sheet initially infected
moi = 0.01;

%time values
final_time = 50;
dt = 0.1;
t_span = 0:dt:final_time;


%number of delay stages
latent_stages = 3;
inf_stages = 1;

%default params
alpha_val = 1.87383799e+00;
beta_val =  1.36340545e-07;
gamma_val = 3.366934e-01;
delta_val = 8.256588e-02;
p_val = exp(1.409457e+01);
c_val = 4.313531e-01;

%ICs
T_0 = 1-moi;
E_CF_0 = 0;
E_CC_0 = 0;
I_CF_0 = 0;
I_CC_0 = 0;
I_r_dead_0 = 0;
I_d_0 = moi;
I_d_dead_0 = 0;
V_0 = 0;


%compartment indices
num_compartments = 2*latent_stages + 3*inf_stages + 6;

target_ind  = 1;
latent_CF_start = target_ind + 1;
latent_CC_start = latent_CF_start + latent_stages;
inf_CF_start = latent_CC_start + latent_stages;
inf_CC_start = inf_CF_start + inf_stages;
inf_r_dead_ind = inf_CC_start + inf_stages;
inf_d_start = inf_r_dead_ind+1;
inf_d_dead_ind = inf_d_start + inf_stages;
v_ind = inf_d_dead_ind + 1;
CF_cum_ind = v_ind + 1;
CC_cum_ind = num_compartments;


%package initial conditions
init_conditions = zeros(num_compartments, 1);

init_conditions(target_ind) = T_0;
init_conditions(latent_CF_start) = E_CF_0;
init_conditions(latent_CC_start) = E_CC_0;
init_conditions(inf_CF_start) = I_CF_0;
init_conditions(inf_CC_start) = I_CC_0;
init_conditions(inf_r_dead_ind) = I_r_dead_0;
init_conditions(inf_d_start)= I_d_0;
init_conditions(inf_d_dead_ind) = I_d_dead_0;
init_conditions(v_ind) = V_0;

%package parameters
params = [alpha_val, beta_val, gamma_val, delta_val, p_val, c_val, ...
    latent_stages, inf_stages];


%% solve
anon_ode = @(t, Y) dual_spread_ODE_model(t, Y, params);
[t_out, tiv_out] = ode45(anon_ode, t_span, init_conditions);


%% plots

%plot compartments
figure
hold on
plot(t_out, tiv_out(:,target_ind)) %target
plot(t_out, sum(tiv_out(:,latent_CF_start:(latent_CF_start+latent_stages-1)),2)) %latent CF
plot(t_out, sum(tiv_out(:,inf_CF_start:(inf_CF_start+inf_stages-1)),2)) %infected CF
plot(t_out, sum(tiv_out(:,latent_CC_start:(latent_CC_start+latent_stages-1)),2)) %latent CC
plot(t_out, sum(tiv_out(:,inf_CC_start:(inf_CC_start+inf_stages-1)),2)) %infected CC
plot(t_out, tiv_out(:,v_ind)/max(tiv_out(:,v_ind))) %virus

legend({'T', 'E_{CF}', 'I_{CF}', 'E_{CC}', 'I_{CC}', 'V'})

ylim([0,1.1])


%now plot prop fluorescent
fluor_prop = sum(tiv_out(:,inf_CF_start:(inf_CF_start+inf_stages-1)),2)...
                        +sum(tiv_out(:,inf_CC_start:(inf_CC_start+inf_stages-1)),2)...
                        +tiv_out(:,inf_r_dead_ind);

figure
plot(t_out, fluor_prop/init_conditions(1))

ylim([0,1.1])

