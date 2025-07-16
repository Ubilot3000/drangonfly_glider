function linearize_massPosSweeper()
    function [eigens, omega_nats, damping_ratios, x_cg] = findEigenvalues(in_x_mass)
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
        
        % Sorting values for consistency
        [omega_nats, sort_idx] = sort(omega_nats);
        eigs_A = eigs_A(sort_idx);
        damping_ratios = damping_ratios(sort_idx);

        % Collecting return values
        eigens = eigs_A(select_idxs);
        omega_nats = omega_nats(select_idxs);
        damping_ratios = damping_ratios(select_idxs);
    end


    % ===========================
    % --- MAIN FUNCTION START ---
    % ===========================
    graph1 = false;
    graph2 = true;

    % Setup iterations
    x_max = 0.100;
    x_min = 0.000;
    num_stations = 10;
    all_x_mass = linspace(x_min, x_max, num_stations);
    all_x_cg = zeros(1, num_stations);

    % Data collection arrays
    omega_short = zeros(1, num_stations);
    omega_phugoid = zeros(1, num_stations);

    damp_short = zeros(1, num_stations);
    damp_phugoid = zeros(1, num_stations);

    eigens_short = zeros(1, num_stations);
    eigens_phugoid = zeros(1, num_stations);
    % Sweep across mass positions
    for j = 1:num_stations
        in_x_mass = all_x_mass(j);
        [eigens, omega_nats, damping_ratios, x_cg] = findEigenvalues(in_x_mass);

        if length(omega_nats) ~= 2
            fprintf("Skipping index %d, invalid eigenvalue count.\n", j);
            continue;
        end

        % Assign based on frequency
        if omega_nats(1) > omega_nats(2)
            short_idx = 1;
            phugoid_idx = 2;
        else
            short_idx = 2;
            phugoid_idx = 1;
        end

        % Data saving
        all_x_cg(j) = x_cg;
        omega_short(j) = omega_nats(short_idx);
        omega_phugoid(j) = omega_nats(phugoid_idx);
        damp_short(j) = damping_ratios(short_idx);
        damp_phugoid(j) = damping_ratios(phugoid_idx);
        eigens_short(j) = eigens(short_idx);
        eigens_phugoid(j) = eigens(phugoid_idx);
    end


    if graph1
        % Find instability boundary (where phugoid goes unstable)
        instability_idx = find(damp_phugoid < 0, 1, 'first');
        x_unstable_start = NaN;
        if ~isempty(instability_idx)
            x_unstable_start = all_x_cg(instability_idx);
        end

        % --- Plotting ---
        figure('Position', [100, 100, 1000, 400]);
        sgtitle("Flight Stability Analyses");
    
        % Frequency Plot
        subplot(1, 2, 1); 
        hold on;
        plot(all_x_cg, omega_short, 'LineWidth', 1.5);
        plot(all_x_cg, omega_phugoid, 'LineWidth', 1.5);
        hold off;
        title('Natural Frequency vs Mass Position');
        ylabel('\omega_n (rad/s)'); xlabel('x_{mass} (m)');
        legend("Short Mode", "Phugoid Mode", 'Location', 'northeast');
        grid on;
    
        % Damping Ratio Plot
        subplot(1, 2, 2); 
        hold on;
    
        % Shade unstable phugoid region
        if ~isnan(x_unstable_start)
            xfill = [x_unstable_start, max(all_x_cg)];
            fill([xfill(1), xfill(2), xfill(2), xfill(1)], [-1, -1, 1, 1], [1, 0.8, 0.8], 'EdgeColor', 'none', "FaceAlpha",0.5);
        end
    
        % Plot damping curves
        plot(all_x_cg, damp_short, 'LineWidth', 1.5);
        plot(all_x_cg, damp_phugoid, 'LineWidth', 1.5);
    
        title('Damping Ratio vs Mass Position');
        ylabel('\zeta'); xlabel('x_{mass} (m)');
        legend("Instability Region", "Short Mode", "Phugoid Mode", 'Location', 'southeast');
        grid on;
        hold off;
    end

    if graph2
        figure('Position', [100, 100, 1000, 400]);
        subplot(1, 2, 1);
        hold on;
        plot(real(eigens_short), imag(eigens_short), 'o-', 'LineWidth', 1.5);
        xline(0, '--k');
        yline(0, '--k');
        hold off;
        axis equal;
        xlabel('Real Part'); 
        ylabel('Imaginary Part');
        title('Short Mode Trajectories');
        sgrid; 

        subplot(1, 2, 2);
        hold on;
        plot(real(eigens_phugoid), imag(eigens_phugoid), 'x-', 'LineWidth', 1.5);
        xline(0, '--k');
        yline(0, '--k');
        hold off;
        axis equal;
        xlabel('Real Part'); 
        ylabel('Imaginary Part');
        title('Phugoid Mode Trajectories');
        sgrid; 
    end
end
