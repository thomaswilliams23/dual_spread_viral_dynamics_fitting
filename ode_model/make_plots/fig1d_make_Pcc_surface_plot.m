

%Makes Figure 1d from the manuscript



%set font of plots
set(0,'defaulttextinterpreter','latex')


%% load in data

%load lookup table data
source_data = '../SIM_DATA_from_ms/lookup_tables/';
load(strcat(source_data, 'Pcc_lookup'));

%maximum values for alpha and beta
alpha_max = 2.5;
beta_max = 2e-6;

%normalised range of lookup values
alpha_range = 0:0.05:1;
beta_range = 0:0.05:1;


%% Pcc surface

figure(1)

%plot contours
Pcc_levels = 0:0.1:1;
[cp2, ch2] = contourf(alpha_max*alpha_range, beta_max*beta_range, prop_CC_data', Pcc_levels);
clabel(cp2, ch2, 'Interpreter', 'tex')

%colour scheme
colormap sky
colorbar

%labels
xlabel('$\alpha$')
ylabel('$\beta$')

%formatting
ylim([0, 1.5e-6])

