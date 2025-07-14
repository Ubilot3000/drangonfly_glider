function rampingMass_Simulate()

    % --- Ramping Parameters ---
    start_time = 5.0;
    start_time_2 = 100;
    x_mass_start = 0.00;
    x_mass_end = 0.050;
    v_mass = 0.02;
    x_mass_end_neg = -x_mass_end;
    v_mass_neg = -v_mass;

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

    % --- Flight Parameters ---
    % Initial conditions

    v0 = 8;
    gamma0 = deg2rad(-5);
    x0 = 0;
    y0 = 16.5;
    theta0 = deg2rad(-1.35 + gamma0);
    q0 = 0;
    
    z0 = [v0; gamma0; x0; y0; theta0; q0];
    tspan = [0, 10];

    distance = abs(x_mass_end - x_mass_start);
    duration = abs(distance / v_mass);
    end_time = start_time + duration;

    % --- Run the ODE Solver ---
    ode_func = @(t, z) rampingMass_ODE(t, z, p, c, start_time, end_time, x_mass_start, x_mass_end, v_mass);
    ode_func_2 = @(t, z) rampingMass_ODE(t, z, p, c, start_time_2, end_time, x_mass_start, x_mass_end, v_mass);
    ode_func_neg = @(t, z) rampingMass_ODE(t, z, p, c, start_time, end_time, x_mass_start, x_mass_end_neg, v_mass_neg);
    [t, z] = ode45(ode_func, tspan, z0);
    [t_2, z_2] = ode45(ode_func_2, tspan, z0);
    [t_neg, z_neg] = ode45(ode_func_neg, tspan, z0);

    % --- Unpack and Plot Results ---
    % Control pitch, mass backwards
    theta=z(:,5);
    x = z (:, 3);
    y = z(:, 4);

    % No pitch control
    theta_2 = z_2(:, 5);
    x_2 = z_2(:, 3);
    y_2 = z_2(:, 4);

    % Control pitch, mass fowards
    theta_neg = z_neg(:, 5);
    x_neg = z_neg(:, 3);
    y_neg = z_neg(:, 4);
    
    
    % Recreate the commanded mass position signal for plotting
    distance = abs(x_mass_end - x_mass_start);
    duration = abs(distance / v_mass);
    end_time = start_time + duration;
    
    x_mass_signal = zeros(size(t));
    for i = 1:length(t)
        if t(i) < start_time
            x_mass_signal(i) = x_mass_start;
        elseif t(i) < end_time
            x_mass_signal(i) = x_mass_start + v_mass * (t(i) - start_time);
        else
            x_mass_signal(i) = x_mass_end;
        end
    end

    figure('Position', [100, 100, 800, 900]);
    
    subplot(3,1,1);
    hold on;
    plot(t, x_mass_signal * 1000, 'LineWidth', 2);
    plot(t, x_mass_signal * -1000, "LineWidth", 2);
    yline(0,"o-", 'LineWidth', 2, "Color", "#EEB220");
    hold off;
    title('Commanded Mass Position (Ramp Input)');
    xlabel('Time (s)');
    ylabel('Mass Position (mm)');
    legend(["Bwd. Pitch Control", "Fwd. Pitch Control", "No Pitch Control"], Orientation="horizontal", Location="north");
    grid on;

    subplot(3,1,2);
    hold on;
    plot(t, rad2deg(theta), 'LineWidth', 2);
    plot(t_neg, rad2deg(theta_neg), "LineWidth", 2);
    plot(t_2, rad2deg(theta_2), 'LineWidth', 2);
    hold off;
    title('Glider Pitch Response to Ramp Input');
    xlabel('Time (s)');
    ylabel('Pitch Angle \theta (deg)');
    grid on;
    legend(["Bwd. Pitch Control", "Fwd. Pitch Control", "No Pitch Control"], Orientation="horizontal", Location="north");

    subplot(3, 1, 3);
    hold on;
    plot(x, y, "LineWidth", 2);
    plot(x_neg, y_neg, "LineWidth", 2);
    plot(x_2, y_2, "LineWidth", 2);
    hold off;
    title("Glider Flight");
    xlabel("X-Position (m)");
    ylabel("Y-Position (m)");
    grid on;
    legend(["Bwd. Pitch Control", "Fwd. Pitch Control", "No Pitch Control"], Orientation="horizontal", Location="north");
    
    sgtitle('Open-Loop Response (Ramp Actuator)', 'FontSize', 16);
end