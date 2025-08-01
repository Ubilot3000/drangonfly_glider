function designSpace_wingDecalageSweep()
    % --- Top-Level Sweep Setup ---
    % Define the array of target C.G. values you want to test
    x_cg_target_array = [0.080, 0.095, 0.110, 0.130, 0.150, 0.170];
    
    % Create the figure window once, before all loops begin
    figure('Position',[200 200 1400 500]);
    
    % Calculate the subplot layout
    num_cgs = length(x_cg_target_array);
    num_cols = min(num_cgs, 3); % At most 3 columns
    num_rows = ceil(num_cgs / num_cols); % Calculate rows needed

    % --- Main Loop to Iterate Through Each Target C.G. ---
    for k = 1:length(x_cg_target_array)
        
        x_cg_target = x_cg_target_array(k); % Select the C.G. for this iteration
        fprintf("Analyzing C.G. %dmm\n", x_cg_target * 1000);
        % --- Per-Analysis Setup ---
        c = getConstructionVector();
        c.x_r_wing = 0.210; % Fix the rear wing position for this study
        
        num_stations = 25;
        x_f_vals = linspace(0.030, 0.120, num_stations);
        decalage_vals = deg2rad(linspace(-5, 2, num_stations));

        % Preallocating arrays for the current C.G. analysis
        L_Ds = zeros(num_stations, num_stations);

        % --- Inner Loops for 2D Sweep (x_f vs. decalage) ---
        for i = 1:num_stations
            x_f = x_f_vals(i);
            for j = 1:num_stations
                decalage = decalage_vals(j);

                % Configuration edits
                c.x_f_wing = x_f;
                c.x_mass = xMassForTargetCG(x_cg_target, c);
                
                p = getParameterVector(c);
                p.i_r = p.i_f + decalage;

                % Finding efficiency
                [L_D_best, ~] = stabilize_computeEfficiency(p);

                % Saving results
                L_Ds(j, i) = L_D_best;
            end
        end

        % --- Plotting for the Current C.G. ---
        % Activate the correct subplot for this iteration
        subplot(num_rows, num_cols, k);
        
        contourf(x_f_vals * 1000, rad2deg(decalage_vals), L_Ds, 15);
        hold on;
        
        % Find and mark the location of the maximum L/D
        [max_LD, linear_idx] = max(L_Ds(:));
        if ~isempty(linear_idx) && isfinite(max_LD)
            [row_idx, col_idx] = ind2sub(size(L_Ds), linear_idx);
            best_x_f = x_f_vals(col_idx);
            best_decalage = rad2deg(decalage_vals(row_idx));
            plot(best_x_f * 1000, best_decalage, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
            text_label = sprintf('  Max L/D = %.2f', max_LD);
            text(best_x_f * 1000, best_decalage, text_label, 'Color', 'red', 'FontWeight', 'bold');
        end
        
        hold off;
        
        % Plot Formatting
        grid on;
        cb = colorbar;
        ylabel(cb, 'Max L/D Ratio');
        xlabel('Front Wing Position, x_f (mm)');
        ylabel('Decalage Angle (i_r - i_f) [degrees]');
        title(sprintf('L/D Performance Map for Fixed C.G. @ %.3f m', x_cg_target));
    end
end