function beta_val = infer_beta_on_peak_time_contour(...
    lookup_table_name, alpha_range, beta_range, alpha_val, peak_time)
%INFER_BETA_ON_PEAK_TIME_CONTOUR
%loads up a lookup table of growth rates for alpha = alpha_range and beta =
%beta_range and infers the value of alpha which preserves the given peak
%time

handle = load(lookup_table_name);
lookup_table = handle.peak_times_mean;  


%infer what the column of growth_rates vs beta would be for this value of
%alpha
constructed_beta_vals = zeros(size(lookup_table,2),1);
for beta_ind = 1:size(lookup_table, 2)
    constructed_beta_vals(beta_ind) = spline(alpha_range, lookup_table(:,beta_ind), alpha_val);
end

%interpolate to piecewise polynomial form of the above array
beta_poly = spline(beta_range, constructed_beta_vals);

%the inferred alpha value solves beta_poly - growth_rate = 0
[~,beta_0_ind] = min(abs(constructed_beta_vals - peak_time));
beta_0 = beta_range(beta_0_ind);
solve_opts = optimset('display', 'off');
beta_val = fsolve(@(x) ppval(beta_poly, x) - peak_time, beta_0, solve_opts);


