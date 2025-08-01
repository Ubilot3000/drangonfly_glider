function designSpace_runSweep()
    % Load base construction & parameters.
    c = getConstructionVector();
    c.m_mass = 0.014;

    % --- Corrected Design Sweep Ranges ---
    % We iterate on PHYSICAL parameters we can actually change.
    i_f_vals     = deg2rad(linspace(0, 6, 4)); % Front wing incidence
    delta_i_vals = deg2rad(round(linspace(-5, -0.5, 20), 1)); % DECALAGE: i_r relative to i_f
    x_f_vals     = round(linspace(0.020, 0.080, 20), 3);
    x_r_vals     = round(linspace(0.200, 0.225, 15), 3);
    x_mass_vals  = round(linspace(0.000, 0.100, 2), 3);  % Position of movable mass
    V_vals       = round(linspace(0.1, 1.2, 8), 3);

    %     % Design sweep ranges
%     i_f_vals   = deg2rad(linspace(0, 10, 5));     % Front wing AoA
%     i_r_vals   = deg2rad(linspace(0, 10, 5));         % Rear wing AoA
%     x_f_vals   = linspace(0.062, 0.067, 15);    % 15      % Front wing position
%     x_r_vals   = linspace(0.210, 0.225, 15);    % 15      % Rear wing position
%     x_cg_vals  = linspace(0.00, 0.140, 5);          % CG position
%     V_vals     = linspace(0.01, 1.2, 5);      % Test velocities (m/s)
% 
%     % Broad test
%     % i_f_vals   = deg2rad(linspace(0, 10, 5));     % Front wing AoA
%     % i_r_vals   = deg2rad(linspace(0, 10, 5));         % Rear wing AoA
%     % x_f_vals   = linspace(0.050, 0.090, 4);    % 15      % Front wing position
%     % x_r_vals   = linspace(0.100, 0.220, 4);    % 15      % Rear wing position
%     % x_cg_vals  = linspace(0.00, 0.140, 5);          % CG position
%     % V_vals     = linspace(0.01, 2, 5);      % Test velocities (m/s)

    % Output storage
    results = [];
    iter_tot = length(i_f_vals) * length(delta_i_vals) * length(x_f_vals) * length(x_r_vals) * length(x_mass_vals);
    counter = 0;
    fprintf("Total Configurations tested: %.0f\n", iter_tot * length(V_vals));

    for i_f = i_f_vals
        for delta_i = delta_i_vals
            fprintf("Percent Complete: %.1f%%\n", counter / iter_tot * 100);
            for x_f = x_f_vals
                for x_r = x_r_vals
                    if x_r <= x_f + 0.05, continue; end % Ensure meaningful separation
                    
                    for x_m = x_mass_vals
                        counter = counter + 1;
                        
                        % --- Build the physical configuration ---
                        c_mod = c;
                        c_mod.i_f = i_f;
                        c_mod.i_r = i_f + delta_i; % i_r is now DEPENDENT on i_f
                        c_mod.x_f_wing = x_f;
                        c_mod.x_r_wing = x_r;
                        c_mod.x_mass = x_m; % Set the physical mass position

                        % Create the parameter vector FROM the physical layout.
                        % This correctly calculates the true x_cg.
                        p = getParameterVector(c_mod);
                        p.i_f = c_mod.i_f;
                        p.i_r = c_mod.i_r;

                        best_L_D = -inf;
                        best_result = [];

                        % --- Velocity Sweep ---
                        for V_test = V_vals
                            try
                                alpha_trim = findTrimCondition(V_test, p);
                                
                                if ~isnan(alpha_trim) % Check for valid trim
                                    % --- Aero Calculation Block ---
                                    alpha_f = alpha_trim + p.i_f;
                                    Cl_f = p.Cl_a * alpha_f;
                                    L_f = 0.5 * p.rho * V_test^2 * p.S_f * Cl_f;
                                    D_f = 0.5 * p.rho * V_test^2 * p.S_f * (p.Cd_0 + p.k_f * Cl_f^2);
                                    
                                    epsilon = (2 * Cl_f) / (pi * p.AR_f) * (1 / (1 + p.k_wings * (p.delta_wings / p.span_f)^2));
                                    
                                    alpha_r_eff = alpha_trim + p.i_r - epsilon;
                                    Cl_r = p.Cl_a * alpha_r_eff;
                                    L_r = 0.5 * p.rho * V_test^2 * p.S_r * Cl_r;
                                    D_r = 0.5 * p.rho * V_test^2 * p.S_r * (p.Cd_0 + p.k_r * Cl_r^2);
                                    
                                    L = L_f + L_r;
                                    D = D_f + D_r;

                                    % --- End Aero Block ---

                                    W = p.m_total * p.g;
                                    if D <= 0 || D >= W
                                        fprintf(...
                                            '  Skipping: V=%.2f  α=%.2f°  L=%.3f N  D=%.3f N  W=%.3f N\n',...
                                            V_test, rad2deg(alpha_trim), L, D, W);
                                        continue;
                                    end
                                    L_D = L / D;

                                    if L_D > best_L_D
                                        best_L_D = L_D;
                                        best_result = [ ...
                                            rad2deg(p.i_f), ...
                                            rad2deg(p.i_r - p.i_f), ...
                                            p.x_f_wing, ...
                                            p.x_r_wing, ...
                                            p.x_cg, ...
                                            V_test, ...
                                            rad2deg(alpha_trim), ...
                                            L_D];
                                    end
                                end
                            catch ME
                                % This catch block is still good practice
                                continue;
                            end
                        end

                        if ~isempty(best_result)
                            results(end+1,:) = best_result;
                        end
                    end
                end
            end
        end
    end

    % --- Analysis and Plotting (similar to your original code) ---
    results_table = array2table(results, ...
        'VariableNames', {'i_f','decalage','x_f','x_r','x_cg','V','alpha','L_D'});

    % =======================
    % === TOP "K" RESULTS ===
    % =======================
    % K = 7;
    % T = sortrows(results_table, 'L_D', 'descend');
    % uLD = unique(T.L_D, 'stable');
    % cut = uLD(min(K, numel(uLD)));
    % topT = T(T.L_D >= cut, :);
    % topT = sortrows(topT, 'L_D', 'ascend');
    % topT.LD_cat = categorical(string(round(topT.L_D, 4)));
    % 
    % figure;
    % parallelplot(topT, ...
    %     'GroupVariable', 'LD_cat', ...
    %     'CoordinateVariables', {'i_f','decalage','x_cg','x_f','x_r','V', 'alpha', 'L_D'});
    % title(sprintf('Top %d Unique L/D Configs', K));


    % ===============================
    % === GROUPING TOP PERCENTILE ===
    % ===============================
    % K = 7;
    % T = sortrows(results_table,'L_D','descend');
    % 
    % % 1) pick your “top 10%” (or “top bin of the full range”) however you like:
    % cutoff = prctile(T.L_D,98);
    % topAll = T(T.L_D>=cutoff,:);
    % 
    % % 2) now create 7 *new* bins just on topAll.L_D
    % edgesTop = linspace(min(topAll.L_D),max(topAll.L_D),K+1);
    % % edgesTop = quantile(topAll.L_D, linspace(0,1,K+1));
    % topAll.LD_bin_sub = discretize(topAll.L_D,edgesTop);
    % topAll = sortrows(topAll, 'L_D', 'ascend');
    % 
    % % 3) plot, grouping by the *sub*‐bin:
    % figure
    % parallelplot(topAll, ...
    %     'GroupVariable','LD_bin_sub', ...
    %     'CoordinateVariables',{'i_f','decalage','x_f','x_r','x_cg','L_D'},...
    %     'Color', "parula"); % removed: 'V','alpha'
    % title('Top‐bin broken into 7 colours')

   
    % =============================================
    % === COMPARING DIFFERENT PERCENTILES PLOTS ===
    % =============================================
    % --- 1) Compute cutoffs for L/D percentiles ---
    p_top    = 99;    % top 1%
    p_mid_lo = 45;    % bottom edge of middle 50%
    p_mid_hi = 55;    % top edge of middle 50%
    p_bot    = 1;     % bottom 1%
    
    L = results_table.L_D;
    cut_top   = prctile(L, p_top);
    cut_midLo = prctile(L, p_mid_lo);
    cut_midHi = prctile(L, p_mid_hi);
    cut_bot   = prctile(L, p_bot);
    
    % --- 2) Grab each subset ---
    T = results_table;
    idx_top  = T.L_D >= cut_top;
    idx_mid  = T.L_D >= cut_midLo  &  T.L_D <= cut_midHi;
    idx_bot  = T.L_D <= cut_bot;
    
    topTbl = T(idx_top, :);
    midTbl = T(idx_mid, :);
    botTbl = T(idx_bot, :);
    
    maxShow = 75;

    % down-sample the middle 50% if it’s too big
    if height(topTbl) > maxShow
        idx = randperm(height(topTbl), maxShow);
        topTbl = topTbl(idx, :);
    end
    
    % down-sample the middle 50% if it’s too big
    if height(midTbl) > maxShow
        idx = randperm(height(midTbl), maxShow);
        midTbl = midTbl(idx, :);
    end
    
    % down-sample the bottom 1% likewise
    if height(botTbl) > maxShow
        idx = randperm(height(botTbl), maxShow);
        botTbl = botTbl(idx, :);
    end
    fprintf("Sizes: %d, %d, %d\n", height(topTbl), height(midTbl), height(botTbl));
    
    % stitch back together and add grouping
    S = [ midTbl; botTbl; topTbl ];
    G = [ repmat("Bot 1%",     height(botTbl),1);
          repmat("Middle 50%", height(midTbl),1);
          repmat("Top 1%",    height(topTbl),1)];
    S.Group = categorical(G);
    
    % final plot
    figure("Position",[100, 100, 1400, 800]);
    parallelplot(S, ...
        'GroupVariable','Group', ...
        'CoordinateVariables',{'i_f','x_r','V','decalage','x_f','x_cg','alpha','L_D'});
    title('Configuration Percentile Efficiency Comparisons');
end