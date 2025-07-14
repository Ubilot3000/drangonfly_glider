function dz = rampingMass_ODE(t, z, p, c, start_time, end_time, x_mass_start, x_mass_end, v_mass)
    % This is the core 6-state ODE function for a ramp mass input.

    % --- Step 1: Determine current mass position based on time and velocity ---
    
    if t < start_time
        % Phase 1: Before movement
        x_mass = x_mass_start;
    elseif t < end_time
        % Phase 2: During movement (ramping)
        x_mass = x_mass_start + v_mass * (t - start_time);
    else
        % Phase 3: After movement
        x_mass = x_mass_end;
    end
    % --- Parameters ---
    rho = p.rho;
    g = p.g;
    c.x_mass = x_mass;
    [m, x_cg, I] = tandemWings_Mech(c);

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