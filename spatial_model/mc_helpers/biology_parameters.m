

%% size and scaling parameters
cells_wide = 50; 
cells_long = 50; 

total_cells = cells_wide * cells_long;

T_tot=4e8;     %this is the total number of cells in the respiratory tract, used to derive beta, etc
prop_of_tract = T_tot/total_cells;  %used to weight contact parameters for the size of this patch



%% model parameters

%infectivity
beta_param = 1.3e-6;
alpha_param = 1;

%eclipse rate
gamma_param = 0.337;

%non-specific infected cell death 
delta_param = 0.00826;

%virus production rate
p = 1.32e6;

%virus clearance
c = 0.431;

