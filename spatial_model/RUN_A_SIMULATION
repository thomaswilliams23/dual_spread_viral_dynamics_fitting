%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                     DUAL SPREAD VIRUS SYSTEM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Implements a CA virus model with two modes of viral spread
% 
% Model structure:
%  - Cells: Target, Eclipse, Infected, Dead, arranged in hexagonal grid
%  - Extracellular virus: ODE


addpath mc_helpers
addpath sim_est_helpers
addpath make_plots/helper_funcs/


%load default params
prms = make_default_param_struct;


%%%%%%%%%%% CHANGE PARAMETERS HERE

%enable visualisation
prms.vis_grid = 1;

%prop of cells initially infected
prms.moi = 0.01;
prms.random_initial = 1;

%number different strains
prms.num_species = 1;

%model parameters
prms.alpha_param = 1.08801559e+00; 
prms.beta_param = 7.19841731e-07;
prms.gamma_param = 3.366934e-01;
prms.delta_param = 8.256588e-02;
prms.p = exp(1.409457e+01);
prms.c = 4.313531e-01;
prms.latent_stages = 3;

%termination conditions
prms.final_time = 50;
prms.stop_when_all_infected = 0;

%time and visualisation intervals
prms.dt = 0.01;
prms.frame_interval = 1;

%%%%%%%%%%%


%specify any data output times
vis_dt = 0.1;
out_times = 0:vis_dt:(prms.final_time-vis_dt);



%run simulation
[sim_out] = single_run_as_func(out_times, prms);


