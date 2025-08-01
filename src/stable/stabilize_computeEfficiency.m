function [L_D_best, V_best] = stabilize_computeEfficiency(p)

    % Variable setup
    num_stations = 20; % A reasonable number for speed
    v_range = linspace(0.01, 10, num_stations); % A more realistic speed range

    L_D_array = zeros(1, num_stations);

    for i = 1:num_stations
        v = v_range(i);

        % --- THIS IS THE CORRECT METHOD ---
        % Find the true trim state (alpha, gamma) for this velocity 'v'
        alpha_stable = findTrimCondition(v, p);

        % If no trim state exists, this design is not flyable at this speed
        if isnan(alpha_stable)
            L_D_array(i) = 0; % Assign a poor value
            continue;
        end
        
        % Calculate total Lift and Drag at the trim condition
        alpha_f = alpha_stable + p.i_f;

        Cl_f = p.Cl_a * alpha_f;
        Cd_f = p.Cd_0 + p.k_f * Cl_f^2;
        L_f = 0.5 * p.rho * v^2 * p.S_f * Cl_f;
        D_f = 0.5 * p.rho * v^2 * p.S_f * Cd_f;
        
        % Downwash
        epsilon = (2 * Cl_f) / (pi * p.AR_f);
        
        % Rear Wing
        alpha_r_eff = alpha_stable + p.i_r - epsilon;
        Cl_r = p.Cl_a * alpha_r_eff;
    
        Cd_r = p.Cd_0 + p.k_r * Cl_r^2;
        L_r = 0.5 * p.rho * v^2 * p.S_r * Cl_r;
        D_r = 0.5 * p.rho * v^2 * p.S_r * Cd_r;

        % Total forces
        L = L_f + L_r;
        D = D_f + D_r;

        % Saving values
        L_D_array(i) = L / D;

        % if L < 0
        %     fprintf("Lift is negative: %.3f\n", L);
        %     fprintf("Crucial Values --> L_f: %.3f, L_r: %.3f, Cl_f: %.3f, Cl_r: %.3f, Alpha_f: %.3f, Alpha_r: %.3f\n", L_f, L_r, Cl_f, Cl_r, alpha_f, alpha_r_eff);
        % end
    end
    
    % Find the best L/D and the velocity at which it occurred
    [L_D_best, best_idx] = max(L_D_array);
    V_best = v_range(best_idx);
end