function [L_total, D_total] = calculateAeroForces(alpha, v, p)
    % This local function calculates the Lift and Drag
    % Output
    %   L_total: lift
    %   D_total: drag


    % Front Wing
    alpha_f = alpha + p.i_f;

    Cl_f = p.Cl_a * alpha_f;
    Cd_f = p.Cd_0 + p.k_f * Cl_f^2;
    L_f = 0.5 * p.rho * v^2 * p.S_f * Cl_f;
    D_f = 0.5 * p.rho * v^2 * p.S_f * Cd_f;
    
    % Downwash
    epsilon = (2 * Cl_f) / (pi * p.AR_f) * (1 / (1 + p.k_wings * (p.delta_wings / p.span_f)^2));
    
    % Rear Wing
    alpha_r_eff = alpha + p.i_r - epsilon;
    Cl_r = p.Cl_a * alpha_r_eff;

    Cd_r = p.Cd_0 + p.k_r * Cl_r^2;
    L_r = 0.5 * p.rho * v^2 * p.S_r * Cl_r;
    D_r = 0.5 * p.rho * v^2 * p.S_r * Cd_r;
    
    L_total = L_f + L_r;
    D_total = D_f + D_r;
end