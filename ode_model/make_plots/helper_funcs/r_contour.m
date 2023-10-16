function beta_val = r_contour(alpha_val, params, stages)
%retrieves the corresponding beta value for a specified alpha value on the
%r contour, given values of the other model parameters
    
    %plug these parameters into the characteristic equation
    char_eq_these_params = @(x) char_eq_staged(x, params, stages);

    %the true r value is the largest solution of the characteristic
    %equation
    r_ref = max(fsolve(char_eq_these_params, 0.4));

    %unpack parameters
    gamma_ref = params(3);
    delta_ref = params(4);
    p_ref = params(5);
    c_ref = params(6);

    %and also number of latent stages
    latent_stages = stages(1);

    %formula for beta given other parameters
    beta_val = (1/p_ref) * (((1 + r_ref/(gamma_ref * latent_stages))^latent_stages)*...
        (c_ref + r_ref)*(delta_ref + r_ref) - alpha_val*(c_ref + r_ref));
end




function y = char_eq_staged(x, params, stages)
%characteristic equation defining r

    alpha_ref = params(1);
    beta_ref = params(2);
    gamma_ref = params(3);
    delta_ref = params(4);
    p_ref = params(5);
    c_ref = params(6);

    latent_stages = stages(1);

    y = ((gamma_ref * latent_stages + x)^latent_stages)*(c_ref + x)*(delta_ref + x) - ...
        (alpha_ref*(c_ref + x) + beta_ref*p_ref)*(latent_stages*gamma_ref)^latent_stages;

end