function [m_total, x_cg, I_total] = tandemWings_Mech(c) 
    % --- Component Parameters ---
    % Nose
    m_nose = c.m_nose;
    x_nose = c.x_nose;
    
    % Mass
    m_mass = c.m_mass;
    x_mass = c.x_mass;

    % Front wing
    m_f_wing = c.m_f_wing;
    x_f_wing = c.x_f_wing;

    % Back wing
    m_r_wing = c.m_r_wing;            
    x_r_wing = c.x_r_wing;

    % Rudder
    m_rudder = c.m_rudder;            
    x_rudder = c.x_rudder;
    
    % CF rod
    m_rod = c.m_rod;
    L_rod = c.L_rod;
    x_rod = L_rod / 2;
    
    
    % --- Total mass & cg ---
    cg_moments = m_nose * x_nose + m_mass * x_mass + m_f_wing * x_f_wing + m_r_wing * x_r_wing + m_rudder * x_rudder + m_rod * x_rod;
    m_total = m_nose + m_mass + m_f_wing + m_r_wing + m_rudder + m_rod;

    x_cg = cg_moments / m_total;


    % --- Moments of Inertia ---
    % Distance from cg
    d_nose = x_nose   - x_cg;
    d_mass = x_mass - x_cg;
    d_f_wing = x_f_wing - x_cg;
    d_r_wing = x_r_wing - x_cg;
    d_rudder = x_rudder - x_cg;
    d_rod = x_rod - x_cg;
    
    % I_total = sum(I_local + m*d^2) for each part
    I_nose = m_nose   * d_nose^2;
    I_mass = m_mass * d_mass^2;
    I_f_wing = m_f_wing * d_f_wing^2;
    I_r_wing = m_r_wing * d_r_wing^2;
    I_rudder = m_rudder * d_rudder^2;
    
    % Rod I
    I_local_rod = (1/12) * m_rod * L_rod^2;
    I_rod = I_local_rod + m_rod * d_rod^2;
    
    % Sum the contributions from all parts
    I_total = I_nose + I_f_wing + I_r_wing + I_rudder + I_rod + I_mass;
end