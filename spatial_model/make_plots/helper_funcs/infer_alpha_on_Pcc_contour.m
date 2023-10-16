function alpha_val = infer_alpha_on_Pcc_contour(...
    lookup_table_name, alpha_range, beta_range, beta_val, true_Pcc)
%INFER_ALPHA_ON_PEAK_TIME_CONTOUR
%loads up a lookup table of growth rates for alpha = alpha_range and beta =
%beta_range and infers the value of alpha which preserves the given peak
%time

handle = load(lookup_table_name);
lookup_table = handle.prop_CC_data;


%find the NaNs first
lookup_nans = isnan(lookup_table);
min_non_nan_alpha = 1;
while sum(sum(lookup_nans(min_non_nan_alpha:end,:)))
    min_non_nan_alpha = min_non_nan_alpha + 1;
end


%infer what the column of growth_rates vs alpha would be for this value of
%beta
constructed_alpha_vals = zeros(size(lookup_table,1),1);
for alpha_ind = min_non_nan_alpha:size(lookup_table, 1)
    constructed_alpha_vals(alpha_ind) = spline(beta_range, lookup_table(alpha_ind,:), beta_val);
end

%interpolate to piecewise polynomial form of the above array
alpha_poly = spline(alpha_range, constructed_alpha_vals);

%the inferred alpha value solves alpha_poly - growth_rate = 0
[~,alpha_0_ind] = min(abs(constructed_alpha_vals - true_Pcc));
alpha_0 = alpha_range(alpha_0_ind);
solve_opts = optimset('display', 'off');
alpha_val = fsolve(@(x) ppval(alpha_poly, x) - true_Pcc, alpha_0, solve_opts);


