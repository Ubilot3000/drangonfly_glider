function pointMass_Simulate()
    % --- Parameters ---
    % Using a struct for parameters for better readability and maintenance.
    p.rho = 1.225;      % air density (kg/m^3)
    p.S = 0.5;          % wing area (m^2)
    p.m = 0.3;          % mass (kg) - Reduced for a "small, light" glider
    p.g = 9.81;         % gravity (m/s^2)
    p.Cl = 0.8;         % constant lift coefficient
    p.Cd = 0.08;        % constant drag coefficient - Increased for a more realistic small model

    % --- Initial Conditions ---
    v0 = 0.001;               % initial speed (m/s)
    theta0 = deg2rad(5);   % small positive initial angle
    y0 = 100;              % initial height (m)
    x0 = 0;                % initial x-position (m)
    z0 = [v0; theta0; y0; x0];
    tspan = [0, 15];      

    % --- Solve ODE ---
    % The anonymous function now passes the struct 'p'
    [t, z] = ode45(@(t, z) pointMass_ODE(t, z, p), tspan, z0);

    % --- Plotting ---
    figure("Position", [100, 100, 1400, 400]);
    sgtitle("Point-Mass Approximation of Glider Flight")

    subplot(1,4,1);
    plot(t, z(:,1), 'b', 'LineWidth', 1.5, "Color", "#0072BE"); 
    title('Speed vs. Time');
    ylabel('Speed V (m/s)'); 
    xlabel('Time (s)');
    grid on;

    subplot(1,4,2);
    plot(t, rad2deg(z(:,2)), 'r', 'LineWidth', 1.5, "Color", "#DA5319");
    title('Flight Path Angle vs. Time');
    ylabel('Flight path angle γ (deg)'); 
    xlabel('Time (s)');
    grid on;

    subplot(1,4,3);
    plot(t, z(:,3), 'k', 'LineWidth', 1.5, "Color", "#EEB220"); 
    title('Altitude vs. Time');
    ylabel('Altitude y (m)'); 
    xlabel('Time (s)'); 
    grid on;
    
    subplot(1, 4, 4);
    plot(z(:, 4), z(:, 3), "g", 'LineWidth', 1.5, "Color", "#7E2F8E");
    title('Glider Trajectory');
    ylabel("Altitude y (m)");
    xlabel("Position x (m)");
    grid on;
end