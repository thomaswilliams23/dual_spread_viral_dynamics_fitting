function interp_val = query_peak_lookup_table(lookup_table_name, ...
    alpha_range, beta_range, alpha_query, beta_query)


%Uses interpolation to query a peak time lookup table at the query points.

handle = load(lookup_table_name);
lookup_table = handle.peak_times_mean;


%interpolate along the alpha query
constructed_vals = zeros(size(lookup_table,2),1);
for beta_ind = 1:size(lookup_table, 2)
    constructed_vals(beta_ind) = spline(alpha_range,lookup_table(:,beta_ind), alpha_query);
end

%build a spline out of interpolated values
beta_spline = spline(beta_range, constructed_vals);

%query spline at desired value
interp_val = ppval(beta_spline, beta_query);





