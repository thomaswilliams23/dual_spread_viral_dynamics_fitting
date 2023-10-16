function interp_val = query_Pcc_lookup_table(lookup_table_name, ...
    alpha_range, beta_range, alpha_query, beta_query)


%Uses interpolation to query a prop_CC lookup table at the query points.
%Allows for a singularity in the bottom left corner of the lookup table.

handle = load(lookup_table_name);
lookup_table = handle.prop_CC_data;



%detect NaNs
min_non_nan_alpha = 1;
min_non_nan_beta = 1;
while sum(sum(isnan(lookup_table(min_non_nan_alpha:end,:))))
    min_non_nan_alpha = min_non_nan_alpha + 1;
end
while sum(sum(isnan(lookup_table(:,min_non_nan_beta:end))))
    min_non_nan_beta = min_non_nan_beta + 1;
end


%bad query
if alpha_query<alpha_range(min_non_nan_alpha) && beta_query<beta_range(min_non_nan_beta)
    disp('ERROR: invalid query points for query_prop_CC_lookup_table');
    interp_val = NaN;
    return
end


%interpolate along the alpha query (follows structure of
%query_peak_lookup_table)
if alpha_query>=alpha_range(min_non_nan_alpha)
    constructed_vals = zeros(size(lookup_table,2),1);
    for beta_ind = 1:size(lookup_table, 2)
        constructed_vals(beta_ind) = spline(alpha_range(min_non_nan_alpha:end),...
            lookup_table(min_non_nan_alpha:end,beta_ind), alpha_query);
    end

    beta_spline = spline(beta_range, constructed_vals);

    interp_val = ppval(beta_spline, beta_query);

else %%%% can't interpolate along alpha query => can interpolate along beta query

    constructed_vals = zeros(size(lookup_table,1),1);
    for alpha_ind = 1:size(lookup_table, 1)
        constructed_vals(alpha_ind) = spline(beta_range(min_non_nan_beta:end),...
            lookup_table(alpha_ind,min_non_nan_beta:end), beta_query);
    end

    alpha_spline = spline(alpha_range, constructed_vals);

    interp_val = ppval(alpha_spline, alpha_query);

end







