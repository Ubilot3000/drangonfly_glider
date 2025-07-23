function linearize_wingPosSweeper()
    function [eigens, omega_nats, damping_ratios, x_cg] = findEigenvalues(in_x_mass, delta_wing)
        % --- Environmental Parameters ---
        p.g = 9.81;
        p.rho = 1.225;
    
        % --- Component Parameters ---
        % Nose
        c.m_nose = 0.004;
        c.x_nose = 0;
    
        % Mass
        c.m_mass = 0.009;
        c.x_mass = in_x_mass;
    
        % Front wing
        c.m_f_wing = 0.004;
        c.x_f_wing = 0.085;% in_x_f_wing; % 0.085;
    
        % Back wing
        c.m_r_wing = 0.006;            
        c.x_r_wing = c.x_f_wing + delta_wing; %in_x_f_wing + 0.055; % 0.140
    
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
    
        % --- Flight Parameters ---
        % Guess equilibrium conditions
    
        v0 = 6;
        gamma0 = deg2rad(-5);
        theta0 = deg2rad(-6);
        q0 = 0;
        

    end
    
    % --- Generate surface plot ---
    [X, Y] = meshgrid(all_x_cg, all_delta_wing);  % X = mass, Y = delta_wing
    
    figure('Position', [100, 100, 1000, 500]);
    hold on;
    surf(X, Y, Z);                 % 3D surface
    % Z_zero = zeros(size(Z));
    % surf(X, Y, Z_zero,'FaceColor', [1, 0.2, 0.2], 'EdgeColor', 'k', 'FaceAlpha', 0.5, 'DisplayName', 'z = 0 plane');
    hold off;
    shading flat;                % smooth shading
    colorbar;                      % adds color scale bar
    clim([min(Z(~isnan(Z))), max(Z(~isnan(Z)))]);                 % optional: fix color scale
    
    xlabel('x_{mass} (m)');
    ylabel('\delta_{wing} (m)');
    zlabel('\zeta_{phugoid}');
    title('Phugoid Damping vs Mass and Wing Separation');
    view(135, 30);                  % adjust view angle
    grid on;
end