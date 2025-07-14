function dz = pidControl_ODE(~, z, p, c, ctrl)
    % --- Parameters ---
    rho = p.rho;
    g = p.g;

    % --- State Variables ---
    v = z(1);
    gamma = z(2); % velocity angle from horizontal
    theta = z(5); % heading angle from horinzontal    
    alpha = theta - gamma; % angle of attack
    q = z(6);
    error_int = z(7);

    % --- PID Controller Logic ---
    % Calculate the error terms
    error_prp = ctrl.theta_goal - theta;
    error_dot = 0 - q;

    % Calculate the PID command for the mass position
    x_mass_cmd = ctrl.Kp * error_prp + ctrl.Ki * error_int + ctrl.Kd * error_dot;
   
    % The mass cannot move beyond the ends of its track.
    x_mass_cmd = max(ctrl.x_mass_min, min(ctrl.x_mass_max, x_mass_cmd));

    c.x_mass = x_mass_cmd;

    % Recalculating plane dynamics
    [m_total, x_cg, I] = tandemWings_Mech(c);
    
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
    dv = -g * sin(gamma) - D_total / m_total;
    dgamma = (L_total / (m_total * v)) - (g * cos(gamma) / v);
    dy = v * sin(gamma);
    dx = v * cos(gamma);
    dtheta = z(6); % q
    dq = M_total / I;
    derror_int = error_prp;

    dz = [dv; dgamma; dy; dx; dtheta; dq; derror_int];
end