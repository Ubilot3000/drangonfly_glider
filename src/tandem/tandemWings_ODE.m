function dz = tandemWings_ODE(~, z, p)
    % --- Parameters ---
    rho = p.rho;
    g = p.g;
    m = p.m_total;
    I = p.I_total;
    x_cg = p.x_cg;

    % --- State Variables ---
    v = z(1);
    gamma = z(2); % velocity angle from horizontal
    theta = z(5); % heading angle from horinzontal    
    alpha = theta - gamma; % angle of attack
    
    % Front Wing
    alpha_f = alpha + p.i_f;

    Cl_f = p.Cl_a * alpha_f;
    Cd_f = p.Cd_0 + p.k_f * Cl_f^2;
    L_f = 0.5 * rho * v^2 * p.S_f * Cl_f;
    D_f = 0.5 * rho * v^2 * p.S_f * Cd_f;
    
    % Downwash
    epsilon = (2 * Cl_f) / (pi * p.AR_f);
    
    % Rear Wing
    alpha_r_eff = alpha + p.i_r - epsilon;
    Cl_r = p.Cl_a * alpha_r_eff;

    Cd_r = p.Cd_0 + p.k_r * Cl_r^2;
    L_r = 0.5 * rho * v^2 * p.S_r * Cl_r;
    D_r = 0.5 * rho * v^2 * p.S_r * Cd_r;
    
    % Summation
    N_f = L_f * cos(alpha) + D_f * sin(alpha);
    N_r = L_r * cos(alpha) + D_r * sin(alpha);
    M_total = (N_f * (x_cg - p.x_f_wing)) + (N_r * (x_cg - p.x_r_wing));
    
    L_total = L_f + L_r;
    D_total = D_f + D_r;

    % E.o.M.
    dv = -g * sin(gamma) - D_total / m;
    dgamma = (L_total / (m * v)) - (g * cos(gamma) / v);
    dy = v * sin(gamma);
    dx = v * cos(gamma);
    dtheta = z(6); % q
    dq = M_total / I;
    
    dz = [dv; dgamma; dx; dy; dtheta; dq];
end