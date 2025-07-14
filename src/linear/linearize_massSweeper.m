function linearize_massSweeper()
    function [eigens, omega_nats, damping_ratios] = findEigenvalues(in_x_mass)
        % --- Environmental Parameters ---
        p.g = 9.81;
        p.rho = 1.225;
    
        % --- Component Parameters ---
        % Nose
        c.m_nose = 0.004;
        c.x_nose = 0;
    
        % Mass
        c.m_mass = 0.004;
        c.x_mass = in_x_mass;
    
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
        % Guess equilibrium conditions
    
        v0 = 6;
        gamma0 = deg2rad(-5);
        theta0 = deg2rad(-6);
        q0 = 0;
        
        % --- Equilibrium ---
        z_guess = [v0; gamma0; theta0; q0];
    
        ode_func = @(z) linearize_ODE(0, z, p);
        options = optimoptions('fsolve','FunctionTolerance',1e-8, 'Display',"none"); % 
        z_eq = fsolve(ode_func, z_guess, options);
    
        % --- Computing Jacobian ---
        % Isolating important states
        % Setup
        n = length(z_eq);
        A = zeros(n);
        dz = 1e-6;
    
        for i = 1:n
            zd1 = z_eq;
            zd2 = z_eq;
            zd1(i) = zd1(i) + dz;
            zd2(i) = zd2(i) - dz;
    
            z1 = zd1;
            z2 = zd2;
    
            f1 = ode_func(z1);
            f2 = ode_func(z2);
    
            A(:, i) = (f1 - f2) / (2 * dz);
        end
    
        % --- Eigenvalue Analysis ---
        eigs_A = eig(A);
        % disp('Eigenvalues:');
        % disp(eigs_A);
    
        select_idxs = find(imag(eigs_A) > 0);
        omega_nats = abs(eigs_A);
        damping_ratios = -real(eigs_A) ./ omega_nats;  % Element-wise division
        
        % Collecting return values
        eigens = eigs_A(select_idxs);
        omega_nats = omega_nats(select_idxs);
        damping_ratios = damping_ratios(select_idxs);
    end

    % Setup iterations
    x_max = 0.200;
    x_min = 0.000;
    num_stations = 45;
    all_x_mass = linspace(x_min, x_max, num_stations);

    % Initialize plot
    figure('Position', [100, 100, 1000, 500]);
    subplot(1, 2, 1); hold on;
    title('Natural Frequency vs Mass Position');
    ylabel('\omega_n (rad/s)'); xlabel('x_{mass} (mm)');
    grid on;

    subplot(1, 2, 2); hold on;
    title('Damping Ratio vs Mass Position');
    ylabel('\zeta'); xlabel('x_{mass} (mm)');
    grid on;

    for j = 1:num_stations
        in_x_mass = all_x_mass(j);
            [eigens, omega_nats, damping_ratios] = findEigenvalues(in_x_mass);
            
            % Plot each eigenvalue mode for this x_mass
            subplot(1, 2, 1);
            scatter(in_x_mass * ones(size(omega_nats)) * 1000, omega_nats, 'filled');

            subplot(1, 2, 2);
            scatter(in_x_mass * ones(size(damping_ratios)) * 1000, damping_ratios, 'filled');
    end
end