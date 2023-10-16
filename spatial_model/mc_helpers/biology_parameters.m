%%  model description


% Pr(T_i -> E_i) = 1 - exp(-(alpha * total_cells * f(x_i) + beta * v_ext) * dt)
% t(E_i -> I_i) ~ gamma(N, gamma_param) (N=20 generally)
% Pr(I_i -> D_i) = 1 - exp(-delta_param * dt)
% 
% d/dt(v_ext) = p_ext * n(infected)/total_cells - c_ext * v_ext


% "infected", or "dead" means the characteristic surface of infected or dead cells
% n(.) means the number of cells with property '.'
% f(x) \in [0,1] is the approximate local infected cell density at location x
%           note: should make this roughly equivalent to n(I)


%% size and scaling parameters
%disp('CAUTION!!! check size of cell grid')
cells_wide = 50; %20; %200;%120;
cells_long = 50; %20; %200;%120;

total_cells = cells_wide * cells_long;

T_tot=4e8;     %this is the total number of cells in the respiratory tract, used to derive beta, etc
prop_of_tract = T_tot/total_cells;  %used to weight contact parameters for the size of this patch



%% model parameters

%infectivity
beta_param = 0.01*prop_of_tract*1.58e-8/24.0; %watered down from Velasco-Hernandez
alpha_param = 5*prop_of_tract*1.58e-8/24.0; %Velasco-Hernandez

%eclipse rate
gamma_param = 2.89;

%non-specific infected cell death 
delta_param = 1.04/24.0;             %ibid

%virus production rate
p = 5.36/24.0;                       %ibid

%virus clearance
c = 2.4/24.0;                %ibid, but stronger outside of cells



%%% drugs:
%oseltamivir (NA-inhibitor, restricts viral budding) 
na_inhib_prod = 0;                    %proportion of production inhibited (default 0)
