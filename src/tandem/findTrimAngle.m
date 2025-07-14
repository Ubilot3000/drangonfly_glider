function [alpha_trim] = findTrimAngle(V, gamma)
        % --- Environmental Parameters ---
    p.g = 9.81;
    p.rho = 1.225;

    % --- Component Parameters ---
    % Nose
    c.m_nose = 0.004;
    c.x_nose = 0;

    % Mass
    c.m_mass = 0.002;
    c.x_mass = 0;

    % Front wing
    c.m_f_wing = 0.004;
    c.x_f_wing = 0.085;

    % Back wing
    c.m_r_wing = 0.006;            
    c.x_r_wing = 0.140;

    % Rudder
    c.m_rudder = 0.002;            
    c.x_rudder = 0.230;
    
    % CF rod
    c.m_rod = 0.005;
    c.L_rod = 0.230;
    c.x_rod = c.L_rod / 2;

    % --- Mass ---
    [m_total, x_cg, I_total] = tandemWings_Mech(c);
    p.m_total = m_total;
    p.x_cg = x_cg;
    p.I_total = I_total;


    % --- Aerodynamic Parameters ---
    % Aerodynamic Coefficients
    p.Cl_a = 2 * pi;     % Lift curve slope (per radian)
    p.Cd_0 = 0.02;       % Parasitic drag
    oswald_e = 0.9;     % Oswald efficiency factor 

    % Front Wing
    p.S_f = 0.007;        % Wing area (m^2)
    p.x_f_wing = c.x_f_wing;
    p.i_f = deg2rad(6);  % Incidence angle (rad)
    p.AR_f = 11.2;
    p.k_f = 1 / (pi * oswald_e * p.AR_f);

    % Rear Wing
    p.S_r = 0.0084;        % Rear wing is larger
    p.x_r_wing = c.x_r_wing;
    p.i_r = deg2rad(3);
    p.AR_r = 9.33;
    p.k_r = 1 / (pi * oswald_e * p.AR_r);


    % Convert flight path angle to radians for calculations
    gamma_rad = deg2rad(gamma);

    % Define the moment function to be passed to fzero
    moment_function = @(alpha) calculateMoment(p, V, alpha);
    
    % Find the alpha (in radians) that makes the moment function zero
    initial_guess_alpha = 0;
    alpha_trim_rad = fzero(moment_function, initial_guess_alpha);

    alpha_trim = rad2deg(alpha_trim_rad);
end