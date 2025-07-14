function rigidBody_Simulate
    % --- Parameters ---
    p.rho = 1.225;      % air density (kg/m^3)
    p.S = 0.5;          % wing area (m^2)
    p.m = 0.5;          % mass (kg)
    p.g = 9.81;         % gravity (m/s^2)
    
    % New Rigid-Body and Aerodynamic Parameters
    p.I_y = 0.05;       % Pitch moment of inertia (kg*m^2), estimate for a small glider
    p.c = 0.25;         % Mean aerodynamic chord (m)
    
    % Aerodynamic Coefficients (example for a stable conventional glider)
    p.Cl0 = 0.2;        % Lift coefficient at alpha=0
    p.Cla = 4.5;        % Lift curve slope (per radian)
    p.Cd0 = 0.04;       % Parasitic drag coefficient
    p.k = 0.05;         % Induced drag factor
    p.Cm0 = 0.05;       % Pitching moment at alpha=0 (typically slightly positive)
    p.Cma = -0.8;       % Pitch stability derivative (must be negative for stability!)
    p.Cmq = -11.0;      % Pitch damping derivative (negative for damping)

    % --- Initial Conditions ---
    % Find a trim state (optional but good practice)
    % For simplicity, we'll start with a small disturbance from level flight.
    v0 = 20;               % initial speed (m/s)
    gamma0 = deg2rad(-2);  % Start in a shallow glide
    y0 = 100;              % initial height (m)
    x0 = 0;                % initial x-position (m)
    theta0 = deg2rad(2);   % Initial pitch angle. alpha0 = theta0 - gamma0 = 4 deg.
    q0 = 0;                % Initial pitch rate (rad/s)
    
    z0 = [v0; gamma0; y0; x0; theta0; q0];
    tspan = [0, 60];

    % --- Solve ODE ---
    % Call the new rigidBodyODE function
    options = odeset('RelTol',1e-6); % Use tighter tolerance for better accuracy
    [t, z] = ode45(@(t, z) rigidBodyODE(t, z, p), tspan, z0, options);

    % --- Post-processing: Calculate AoA ---
    gamma = z(:,2); 
    theta = z(:,5);
    alpha_deg = rad2deg(theta - gamma);

    % --- Plotting ---
    figure("Position", [100, 100, 1200, 800]);

    subplot(3,2,1);
    plot(t, z(:,1), 'b', 'LineWidth', 1.5); 
    title('Speed vs. Time'); ylabel('Speed V (m/s)'); grid on;

    subplot(3,2,2);
    plot(t, rad2deg(z(:,2)), 'r', 'LineWidth', 1.5);
    title('Flight Path Angle vs. Time'); ylabel('γ (deg)'); grid on;

    subplot(3,2,3);
    plot(t, z(:,3), 'k', 'LineWidth', 1.5); 
    title('Altitude vs. Time'); ylabel('Altitude y (m)'); xlabel('Time (s)'); grid on;
    
    subplot(3,2,4);
    plot(z(:,4), z(:,3), 'g', 'LineWidth', 1.5);
    title('Glider Trajectory'); ylabel("Altitude y (m)"); xlabel("Position x (m)");
    axis equal; grid on;

    subplot(3,2,5);
    plot(t, alpha_deg, 'm', 'LineWidth', 1.5);
    title('Angle of Attack vs. Time'); ylabel('α (deg)'); xlabel('Time (s)'); grid on;

    subplot(3,2,6);
    plot(t, rad2deg(z(:,5)), 'c', 'LineWidth', 1.5);
    title('Pitch Angle vs. Time'); ylabel('θ (deg)'); xlabel('Time (s)'); grid on;
end