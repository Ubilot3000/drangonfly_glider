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
        
        % Sorting values for consistency
        [omega_nats, sort_idx] = sort(omega_nats);
        eigs_A = eigs_A(sort_idx);
        damping_ratios = damping_ratios(sort_idx);

        % Collecting return values
        eigens = eigs_A(select_idxs);
        omega_nats = omega_nats(select_idxs);
        damping_ratios = damping_ratios(select_idxs);
    end

    % Setup mass iterations
    x_mass_max = 0.100;
    x_mass_min = 0.000;
    num_stations_mass = 30;
    all_x_mass = linspace(x_mass_min, x_mass_max, num_stations_mass);
    all_x_cg = zeros(1, num_stations_mass);

    % Setup wing iterations
    delta_wing_max = 0.010;
    delta_wing_min = 0.150;
    num_stations_wing = num_stations_mass;
    all_delta_wing = linspace(delta_wing_min, delta_wing_max, num_stations_wing);

    % --- Preallocate matrix for phugoid damping ratios ---
    Z = NaN(num_stations_wing, num_stations_mass);  % Z(row: delta_wing, col: x_mass)
    
    % --- Fill Z matrix with phugoid damping values ---
    for r = 1:num_stations_wing
        in_delta_wing = all_delta_wing(r);
        for j = 1:num_stations_mass
            in_x_mass = all_x_mass(j);
            [~, omega_nats, damping_ratios, x_cg] = findEigenvalues(in_x_mass, in_delta_wing);
            all_x_cg(j) = x_cg;
            
            % Check if valid result
            if length(omega_nats) ~= 2
                continue;
            end
    
            % Assign mode indices based on frequency
            if omega_nats(1) > omega_nats(2)
                short_idx = 1;
                phugoid_idx = 2;
            else
                short_idx = 2;
                phugoid_idx = 1;
            end
    
            % Store damping ratio into matrix
            Z(r, j) = damping_ratios(phugoid_idx);
        end
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