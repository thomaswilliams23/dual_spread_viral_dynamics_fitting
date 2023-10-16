function rsse = data_to_sim_dist(data, model)


%computes mean root sum squared errors between observed and model data

rsse = norm(data-model);
