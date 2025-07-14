function dz = pointMass_ODE(~, z, p)
    % Unpack state variables for clarity
    v = z(1);       % speed (m/s)
    gamma = z(2);   % flight path angle (rad)
    % z(3) is y (altitude)
    % z(4) is x (position)

    % Unpack parameters from the struct for clarity and robustness
    rho = p.rho; 
    S = p.S; 
    m = p.m; 
    g = p.g;
    Cl = p.Cl; 
    Cd = p.Cd;

    % Aerodynamic forces
    L = 0.5 * rho * v^2 * S * Cl;
    D = 0.5 * rho * v^2 * S * Cd;

    % Equations of motion
    dv = -g * sin(gamma) - D / m;
    dgamma = (L / (m * v)) - (g * cos(gamma) / v);
    dy = v * sin(gamma);
    dx = v * cos(gamma);

    dz = [dv; dgamma; dy; dx];
end