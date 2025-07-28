function designSpace_runSweep()
    % Load base construction & parameters
    c = getConstructionVector();

    % Design sweep ranges
    i_f_vals   = deg2rad(linspace(0, 10, 4));     % Front wing AoA
    i_r_vals   = deg2rad(linspace(0, 10, 4));         % Rear wing AoA
    x_f_vals   = linspace(0.063, 0.072, 12);    % 15      % Front wing position
    x_r_vals   = linspace(0.100, 0.130, 12);    % 15      % Rear wing position
    x_cg_vals  = linspace(0.00, 0.140, 5);          % CG position
    V_vals     = linspace(0.5, 2.0, 8);      % Test velocities (m/s)

    % Output storage
    results = [];

    for i = 1:length(i_f_vals)
        fprintf("Percent Complete: %.1f\n", i / length(i_f_vals) * 100);
        for j = 1:length(i_r_vals)
            for k = 1:length(x_f_vals)
                for m = 1:length(x_r_vals)
                    for n = 1:length(x_cg_vals)

                        % Build config
                        c_mod = c;
                        c_mod.i_f = i_f_vals(i);
                        c_mod.i_r = i_r_vals(j);
                        c_mod.x_f_wing = x_f_vals(k);
                        c_mod.x_r_wing = x_r_vals(m);
                        c_mod.x_cg = x_cg_vals(n);
                        p = getParameterVector(c_mod);

                        best_L_D = -inf;
                        best_result = [];

                        % ✅ NOW start velocity sweep
                        for v = 1:length(V_vals)
                            V_test = V_vals(v);

                            % Just for debug:
                            % fprintf('Trying config: i_f = %.1f°, i_r = %.1f°, x_f = %.2f, x_r = %.2f, x_cg = %.2f, V = %.2f\n', ...
                                % rad2deg(c_mod.i_f), rad2deg(c_mod.i_r), c_mod.x_f_wing, c_mod.x_r_wing, c_mod.x_cg, V_test);

                            try
                                alpha_trim = findTrimCondition(V_test, p);
                                success = true;
                                % alpha_trim = deg2rad(5);
                                if success
                                    % fprintf('Trim success at alpha_trim = %.2f°\n', rad2deg(alpha_trim));

                                    % Front Wing
                                    alpha_f = alpha_trim + p.i_f;
                                    Cl_f = p.Cl_a * alpha_f;
                                    Cd_f = p.Cd_0 + p.k_f * Cl_f^2;
                                    L_f = 0.5 * p.rho * V_test^2 * p.S_f * Cl_f;
                                    D_f = 0.5 * p.rho * V_test^2 * p.S_f * Cd_f;

                                    % Downwash with spacing dependency
                                    epsilon = (2 * Cl_f) / (pi * p.AR_f) * ...
                                        (1 / (1 + p.k_wings * (p.delta_wings / p.span_f)^2));

                                    % Rear Wing
                                    alpha_r_eff = alpha_trim + p.i_r - epsilon;
                                    Cl_r = p.Cl_a * alpha_r_eff;
                                    Cd_r = p.Cd_0 + p.k_r * Cl_r^2;
                                    L_r = 0.5 * p.rho * V_test^2 * p.S_r * Cl_r;
                                    D_r = 0.5 * p.rho * V_test^2 * p.S_r * Cd_r;

                                    % Total Lift and Drag
                                    L = L_f + L_r;
                                    D = D_f + D_r;

                                    if D <= 0
                                        fprintf("Triggedred skip condition")
                                        continue
                                    end

                                    % L/D Ratio
                                    L_D = L / D;

                                    if L_D > best_L_D
                                        best_L_D = L_D;
                                        best_result = [ ...
                                            rad2deg(c_mod.i_f), ...
                                            rad2deg(c_mod.i_r), ...
                                            c_mod.x_f_wing, ...
                                            c_mod.x_r_wing, ...
                                            c_mod.x_cg, ...
                                            V_test, ...
                                            rad2deg(alpha_trim), ...
                                            L_D];
                                    end
                                end
                            catch ME
                                fprintf('\n⚠️ Error in findTrimConditions at V = %.2f:\n%s\n', V_test, ME.message);
                                continue
                            end
                        end

                        % Save best if found
                        if ~isempty(best_result)
                            results(end+1,:) = best_result;
                        end
                    end
                end
            end
        end
    end

    % Convert to table
    results_table = array2table(results, ...
        'VariableNames', {'i_f','i_r','x_f','x_r','x_cg','V','alpha','L_D'});
    % Save to files
    % writetable(results_table, 'designSweep_results.csv');
    % save('designSweep_results.mat', 'results_table');

    % Cutoff
    cutoff = prctile(results_table.L_D, 90);

    % Best L_D
    fprintf("Best Efficiency: %.3f\n", max(results_table.L_D));

    % Filter to only those rows whose L_D is ≥ the cutoff
    top10_table = results_table(results_table.L_D >= cutoff, :);

    % Plot Parallel Coordinates
    figure;
    parallelplot(top10_table, 'GroupVariable', 'L_D');
    title('Design Sweep: Parallel Coordinates (Colored by L/D)');
end
