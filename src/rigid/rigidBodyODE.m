% We rename the file to reflect the new model
function dz = rigidBodyODE(~, z, p)
    % Unpack state variables for clarity
    v = z(1);       % speed (m/s)
    gamma = z(2);   % flight path angle (rad)
    % y = z(3)      % altitude (m)
    % x = z(4)      % position (m)
    theta = z(5);   % pitch angle (rad)
    q = z(6);       % pitch rate (rad/s)

    % Unpack parameters from the struct
    rho = p.rho; 
    S = p.S; 
    m = p.m; 
    g = p.g;
    I_y = p.I_y; 
    c = p.c;
    Cl0 = p.Cl0; 
    Cla = p.Cla;
    Cd0 = p.Cd0; 
    k = p.k;
    Cm0 = p.Cm0; 
    Cma = p.Cma; 
    Cmq = p.Cmq;

    % --- Core Calculation of the Rigid-Body Model ---
    % Calculate Angle of Attack (the crucial link)
    alpha = theta - gamma;

    % Calculate aerodynamic coefficients based on the current state
    Cl = Cl0 + Cla * alpha;
    Cd = Cd0 + k * Cl^2;
    
    % The q_dimless term is a non-dimensional pitch rate.
    q_dimless = q * c / (2 * v);
    Cm = Cm0 + Cma * alpha + Cmq * q_dimless;

    % Calculate Forces and Moments
    L = 0.5 * rho * v^2 * S * Cl;
    D = 0.5 * rho * v^2 * S * Cd;
    M = 0.5 * rho * v^2 * S * c * Cm;

    % --- Equations of Motion ---
    % Original force equations (now with non-constant L and D)
    dv = -g * sin(gamma) - D / m;
    dgamma = (L / (m * v)) - (g * cos(gamma) / v);
    dy = v * sin(gamma);
    dx = v * cos(gamma);

    % New rotational equations
    dtheta = q;
    dq = M / I_y;

    dz = [dv; dgamma; dy; dx; dtheta; dq];
end