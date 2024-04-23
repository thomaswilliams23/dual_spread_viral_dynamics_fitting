

%Makes Figure 8 from the manuscript



%set font of plots
set(0,'defaulttextinterpreter','latex')


%% load in data

%load lookup table data
source_data = '../SIM_DATA_from_ms/lookup_tables/';
load(strcat(source_data, 'alpha_range')); %normalised between 0 and 1
load(strcat(source_data, 'beta_range')); %normalised between 0 and 1
load(strcat(source_data, 'Pcc_lookup'));
load(strcat(source_data, 'tpeak_lookup'));

%maximum values for alpha and beta
alpha_max = 10;
beta_max = 5e-6;


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



%% tpeak surface

figure(2)

%plot contours
tpeak_levels = [6, 8, 10, 12, 14, 16, 18, 20, 30];
[cp, ch] = contourf(alpha_max*alpha_range, beta_max*beta_range, peak_times_mean', tpeak_levels);
clabel(cp, ch, 'Interpreter', 'tex')

%colour scheme
colormap summer
colorbar

%labels
xlabel('$\alpha$')
ylabel('$\beta$')

%formatting
ylim([0, 1.5e-6])


