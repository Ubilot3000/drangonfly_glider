function [x, y] = tandemWings_Simulate(IN_x_mass)
    % --- Environmental Parameters ---
    p.g = 9.81;
    p.rho = 1.225;

    % --- Component Parameters ---
    % Nose
    c.m_nose = 0.004;
    c.x_nose = 0;

    % Mass
    c.m_mass = IN_x_mass;
    c.x_mass = 0.01848;

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

    % --- Flight Parameters ---
    % Initial conditions
    v0 = 9;
    gamma0 = deg2rad(2);
    x0 = 0;
    y0 = 100;
    theta0 = deg2rad(-2);
    q0 = 0;
    
    z0 = [v0; gamma0; x0; y0; theta0; q0];
    tspan = [0, 60];

    [t, z] = ode45(@(t, z) tandemWings_ODE(t, z, p), tspan, z0);

    % disp([x_cg, IN_x_mass]);

    % --- Plotting ---
    % Plotting setup
    figure("Position", [100, 100, 1000, 800]);
    sgtitle('Tandem-Wing Glider Flight Analysis', 'FontSize', 16);
    v=z(:,1); gamma=z(:,2); x=z(:,3); y=z(:,4); theta=z(:,5);
    alpha = theta - gamma;
    alpha_f = alpha + p.i_f;
    Cl_f = p.Cl_a * alpha_f;
    epsilon = (2 * Cl_f) / (pi * p.AR_f);
    alpha_r_eff = alpha + p.i_r - epsilon;

    % Extra plotting calculations
    p_fit = polyfit(t, y, 2);
    func_y_vs_t = @ (t) p_fit(1) .* t.^2 + p_fit(2) .* t + p_fit(3);
    t_y_0 = fzero(func_y_vs_t, 0.5);
    fprintf("Intercepts floor at t = %.3f\n", t_y_0);

    subplot(3,2,1); 
    plot(t, v, 'b', 'LineWidth', 1.5); 
    title('Speed'); 
    ylabel('V (m/s)'); 
    xlabel("Time (s)")
    grid on;

    subplot(3,2,2); 
    plot(t, rad2deg(theta), 'r', 'LineWidth', 1.5); 
    title('Pitch Angle'); 
    ylabel('θ (deg)'); 
    xlabel("Time (s)")
    grid on;

    subplot(3,2,3);
    % hold on;
    plot(t, y, 'k', 'LineWidth', 1.5); 
    % plot(t, func_y_vs_t(t), "LineWidth", 1.5);
    % hold off;
    title('Altitude'); 
    ylabel('y (m)'); 
    xlabel('Time (s)'); 
    grid on;

    subplot(3,2,4); 
    plot(x, y, 'g', 'LineWidth', 1.5); 
    title('Trajectory'); 
    ylabel("Altitude (m)"); 
    xlabel("Position (m)"); 
    axis equal; 
    grid on;

    subplot(3,2,5); 
    plot(t, rad2deg(alpha_f), 'LineWidth', 1.5); 
    hold on; 
    plot(t, rad2deg(alpha_r_eff), 'LineWidth', 1.5); 
    title('Effective Wing AoA'); 
    ylabel('α (deg)'); 
    xlabel('Time (s)'); 
    legend('Front Wing', 'Rear Wing', 'Location', 'best'); 
    grid on;

    subplot(3,2,6); 
    plot(t, rad2deg(epsilon), 'LineWidth', 1.5); 
    title('Downwash Angle'); 
    ylabel('\epsilon (deg)'); 
    xlabel('Time (s)'); 
    grid on;

    x = 0; y = 0;


end