function linearize_massPosSweeper()
    % ===========================
    % --- MAIN FUNCTION START ---
    % ===========================
    graph1 = false;
    graph2 = true;

    x_min = 0; 
    x_max = 0.110; 
    num_stations = 50;
    all_x_cg = linspace(x_min, x_max, num_stations);
    
    
    
    % Preallocating arrays
    omega_short    = NaN(1,num_stations);
    omega_phugoid  = NaN(1,num_stations);
    damp_short     = NaN(1,num_stations);
    damp_phugoid   = NaN(1,num_stations);
    eigens_short   = NaN(1,num_stations);
    eigens_phugoid = NaN(1,num_stations);

    % Model construction
    c = getConstructionVector();
    c.x_r_wing = 0.200;
    v0 = 4;    % a realistic trim speed


    for j = 1:num_stations
        target_x_cg = all_x_cg(j);
        c.x_mass = xMassForTargetCG(target_x_cg, c);
        p = getParameterVector(c);

        %--- Call the trimmed‐based eigen solver --- 
        [eigens, omega_nats, damping_ratios] = findEigenvalues(p, v0);

        % If it failed to return exactly two modes, skip
        if numel(omega_nats) ~= 2
            fprintf("Skipping idx %d: got %d modes\n", j, numel(omega_nats));
            continue;
        end

        [~, short_idx]   = max(damping_ratios);
        phugoid_idx     = 3 - short_idx;  % the other one
    
        % store in your arrays
        omega_short(j)    = omega_nats(short_idx);
        damp_short(j)     = damping_ratios(short_idx);
        eigens_short(j)   = eigens(short_idx);
    
        omega_phugoid(j)  = omega_nats(phugoid_idx);
        damp_phugoid(j)   = damping_ratios(phugoid_idx);
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
        cmap = jet(length(eigens_short));
        for k = 1:length(eigens_short)-1
          plot(real(eigens_short(k:k+1)), ...
               imag(eigens_short(k:k+1)), ...
               'LineWidth',2, 'Color',cmap(k,:));
        end
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
        cmap = jet(length(eigens_phugoid));
        hold on
        for k = 1:length(eigens_phugoid)-1
            plot(real(eigens_phugoid(k:k+1)), ...
                imag(eigens_phugoid(k:k+1)), ...
                'LineWidth',2, 'Color',cmap(k,:));
        end
        colormap(jet);       % set the colormap
        c = colorbar;        % draw the colorbar
        c.Label.String = 'Sweep Step';  % optional label
        clim([1, numel(eigens_short)]);
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

 % function [eigens, omega_nats, damping_ratios] = findEigenvalues_2(p)
    %     % --- Flight Parameters ---
    %     % Guess equilibrium conditions
    % 
    %     v0 = 6;
    %     gamma0 = deg2rad(-5);
    %     theta0 = deg2rad(-6);
    %     q0 = 0;
    % 
    %     % --- Equilibrium ---
    %     z_guess = [v0; gamma0; theta0; q0];
    % 
    %     ode_func = @(z) linearize_ODE(0, z, p);
    %     options = optimoptions('fsolve','FunctionTolerance',1e-8, 'Display',"none"); % 
    %     z_eq = fsolve(ode_func, z_guess, options);
    % 
    %     % --- Computing Jacobian ---
    %     % Isolating important states
    %     % Setup
    %     n = length(z_eq);
    %     A = zeros(n);
    %     dz = 1e-6;
    % 
    %     for i = 1:n
    %         zd1 = z_eq;
    %         zd2 = z_eq;
    %         zd1(i) = zd1(i) + dz;
    %         zd2(i) = zd2(i) - dz;
    % 
    %         z1 = zd1;
    %         z2 = zd2;
    % 
    %         f1 = ode_func(z1);
    %         f2 = ode_func(z2);
    % 
    %         A(:, i) = (f1 - f2) / (2 * dz);
    %     end
    % 
    %     % --- Eigenvalue Analysis ---
    %     eigs_A = eig(A);
    %     % disp('Eigenvalues:');
    %     % disp(eigs_A);
    % 
    %     select_idxs = find(imag(eigs_A) > 0);
    %     omega_nats = abs(eigs_A);
    %     damping_ratios = -real(eigs_A) ./ omega_nats;  % Element-wise division
    % 
    %     % Sorting values for consistency
    %     [omega_nats, sort_idx] = sort(omega_nats);
    %     eigs_A = eigs_A(sort_idx);
    %     damping_ratios = damping_ratios(sort_idx);
    % 
    %     % Collecting return values
    %     eigens = eigs_A(select_idxs);
    %     omega_nats = omega_nats(select_idxs);
    %     damping_ratios = damping_ratios(select_idxs);
    % end


