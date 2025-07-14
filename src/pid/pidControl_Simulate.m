function pidControl_Simulate()
    % --- Environmental Parameters ---
    p.g = 9.81;
    p.rho = 1.225;

    % --- Component Parameters ---
    % Nose
    c.m_nose = 0.004;
    c.x_nose = 0;

    % Mass
    c.m_mass = 0.004;
    % c.x_mass = 0.01848;

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

    % --- Controller Parameters ---
    % Gains and command
    ctrl.theta_goal = deg2rad(-30); % goal theta
    ctrl.Kp = 0.7;
    ctrl.Ki = 0.45;
    ctrl.Kd = 1.4;

    % Control limits
    ctrl.x_mass_max = 0.100;
    ctrl.x_mass_min = 0.000;

    % --- Flight Parameters ---
    % Initial conditions

    v0 = 8;
    gamma0 = deg2rad(0);
    x0 = 0;
    y0 = 16.5;
    theta0 = deg2rad(0);
    q0 = 0;
    int_error0 = 0;

    % Simluation inputs
    z0 = [v0; gamma0; x0; y0; theta0; q0; int_error0];
    tspan = [0, 30];

    % --- ODE Solver ---
    ode_func = @(t, z) pidControl_ODE(t, z, p, c, ctrl);
    [t, z] = ode45(ode_func, tspan, z0);

        % --- Unpack and Plot Results ---
    theta=z(:,5);
    
    % Recreate the mass position signal by running the controller logic again
    x_mass_signal = zeros(size(t));
    for i = 1:length(t)
        error_prp = ctrl.theta_goal - z(i,5); % theta
        error_dot = 0 - z(i,6); % q is the rate of change of theta
        error_int = z(i,7);
        x_mass_cmd = ctrl.Kp * error_prp + ctrl.Ki * error_int + ctrl.Kd * error_dot;
        x_mass_signal(i) = max(ctrl.x_mass_min, min(ctrl.x_mass_max, x_mass_cmd)); % Saturate
    end

    figure('Position', [100, 100, 800, 600]);
    
    subplot(2,1,1);
    plot(t, x_mass_signal * 1000, 'LineWidth', 2);
    title('PID Controller Output');
    xlabel('Time (s)'); ylabel('Commanded Mass Position (mm)'); grid on;

    subplot(2,1,2);
    plot(t, rad2deg(theta), 'LineWidth', 2);
    hold on;
    plot(t, ones(size(t))*rad2deg(ctrl.theta_goal), '--', 'LineWidth', 1.5);
    title('Glider Pitch Response with PID Control');
    xlabel('Time (s)'); ylabel('Pitch Angle \theta (deg)');
    legend('Actual Pitch', 'Desired Pitch'); grid on;
    
    sgtitle('Closed-Loop PID Control', 'FontSize', 16);
end