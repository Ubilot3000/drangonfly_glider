function c = getConstructionVector()
    % Estabilishes the vector detailing the construction of the glider
    % --- Component Parameters ---
    % Nose
    c.m_nose = 0.004;
    c.x_nose = 0;

    % Mass
    c.m_mass = 0.004;
    c.x_mass = 0.01848;

    % Front wing
    c.m_f_wing = 0.004;
    c.x_f_wing = 0.085;
    c.i_f = deg2rad(6);

    % Back wing
    c.m_r_wing = 0.006;            
    c.x_r_wing = 0.140;
    c.i_r = deg2rad(3);

    % Rudder
    c.m_rudder = 0.002;            
    c.x_rudder = 0.230;
    
    % CF rod
    c.m_rod = 0.005;
    c.L_rod = 0.230;
    c.x_rod = c.L_rod / 2;
end