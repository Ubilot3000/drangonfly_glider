function [Cm_alpha] = stabilize_computeCma(x_cg, c, v, alpha)
    % Establishing glider model parameters
    p = getParameterVector(c);

    % Because only moment analysis is being carried out, we can just set
    % c.g. where we want it and ignore what the physical requirement for
    % its placement would be
    p.x_cg = x_cg;

    % --- Initial Conditions ---
    % Set flight
    alpha_delta = deg2rad(0.1);

    % Calculate Cm_alpha
    M_1 = calculateMoment(p, v, alpha);
    M_2 = calculateMoment(p, v, alpha + alpha_delta);

    % Moment coefficients
    Q = 0.5 * p.rho * v^2; % Dynamic pressure
    Cm_1 = M_1 / (Q * p.S_tot * p.c_bar);
    Cm_2 = M_2 / (Q * p.S_tot * p.c_bar);

    % Moment alpha coefficient
    Cm_alpha = (Cm_2 - Cm_1) / alpha_delta;
end