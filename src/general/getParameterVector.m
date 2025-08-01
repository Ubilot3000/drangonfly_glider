function p = getParameterVector(c)
    % Establishes all parameters to due with glider flight dyanmics
    % --- Environmental Parameters ---
    p.g = 9.81;
    p.rho = 1.225;
    
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
    p.i_f = c.i_f;  % Incidence angle (rad)
    p.AR_f = 11.2;
    p.k_f = 1 / (pi * oswald_e * p.AR_f);
    p.span_f = sqrt(p.AR_f * p.S_f);

    % Rear Wing
    p.S_r = 0.0084;        % Rear wing is larger
    p.x_r_wing = c.x_r_wing;
    p.i_r = c.i_r;
    p.AR_r = 9.33;
    p.k_r = 1 / (pi * oswald_e * p.AR_r);

    % Glider properties
    p.S_tot = p.S_f + p.S_r;
    p.c_bar = (sqrt(p.S_f / p.AR_f) * p.S_f + sqrt(p.S_r / p.AR_r) * p.S_r) / p.S_tot;
    
    % Wing interactions
    p.delta_wings = c.x_r_wing - c.x_f_wing; 
    p.k_wings = 8;
end