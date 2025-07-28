function [M_total] = calculateMoment(p, V, alpha)
    % Calculates the total pitching moment for a given velocity and alpha
    % --- State Variables ---
    % Front Wing
    alpha_f = alpha + p.i_f;

    Cl_f = p.Cl_a * alpha_f;
    Cd_f = p.Cd_0 + p.k_f * Cl_f^2;
    L_f = 0.5 * p.rho * V^2 * p.S_f * Cl_f;
    D_f = 0.5 * p.rho * V^2 * p.S_f * Cd_f;
    
    % Downwash
    epsilon = (2 * Cl_f) / (pi * p.AR_f) * (1 / (1 + p.k_wings * (p.delta_wings / p.span_f)^2));
    
    % Rear Wing
    alpha_r_eff = alpha + p.i_r - epsilon;
    Cl_r = p.Cl_a * alpha_r_eff;

    Cd_r = p.Cd_0 + p.k_r * Cl_r^2;
    L_r = 0.5 * p.rho * V^2 * p.S_r * Cl_r;
    D_r = 0.5 * p.rho * V^2 * p.S_r * Cd_r;
    
    % Summation
    N_f = L_f * cos(alpha) + D_f * sin(alpha);
    N_r = L_r * cos(alpha) + D_r * sin(alpha);
    M_total = (N_f * (p.x_cg - p.x_f_wing)) + (N_r * (p.x_cg - p.x_r_wing));
end
