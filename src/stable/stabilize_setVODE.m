function dz = linearize_setVODE(~, z, p, v0)
    % z_guess = [v0; gamma0; theta0; q0];
    % --- Parameters ---
    rho = p.rho;
    g = p.g;
    m = p.m_total;
    I = p.I_total;
    x_cg = p.x_cg;

    % --- State Variables ---
    gamma = z(1); % velocity angle from horizontal
    theta = z(2); % heading angle from horinzontal    
    alpha = theta - gamma; % angle of attack
    
    % Front Wing
    alpha_f = alpha + p.i_f;

    Cl_f = p.Cl_a * alpha_f;
    Cd_f = p.Cd_0 + p.k_f * Cl_f^2;
    L_f = 0.5 * rho * v0^2 * p.S_f * Cl_f;
    D_f = 0.5 * rho * v0^2 * p.S_f * Cd_f;
    
    % Downwash
    epsilon = (2 * Cl_f) / (pi * p.AR_f);
    
    % Rear Wing
    alpha_r_eff = alpha + p.i_r - epsilon;
    Cl_r = p.Cl_a * alpha_r_eff;

    Cd_r = p.Cd_0 + p.k_r * Cl_r^2;
    L_r = 0.5 * rho * v0^2 * p.S_r * Cl_r;
    D_r = 0.5 * rho * v0^2 * p.S_r * Cd_r;
    
    % Summation
    N_f = L_f * cos(alpha) + D_f * sin(alpha);
    N_r = L_r * cos(alpha) + D_r * sin(alpha);
    M_total = (N_f * (x_cg - p.x_f_wing)) + (N_r * (x_cg - p.x_r_wing));
    
    L_total = L_f + L_r;
    D_total = D_f + D_r;

    % E.o.M.
    dgamma = (L_total / (m * v0)) - (g * cos(gamma) / v0);
    dtheta = z(3); % q
    dq = M_total / I;
    
    dz = [dgamma; dtheta; dq];
end