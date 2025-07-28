function stabilize_efficiencyRange()
    num_stations = 25;

    x_mass_range = linspace(0.00, 0.190, num_stations);
    x_f_wing_range = linspace(0.030, 0.120, num_stations);
    x_r_wing_range = linspace(0.080, 0.170, num_stations);

    % Preallocate arrays.
    results_count = 0;
    total_iterations = num_stations^3;
    X_cg = zeros(total_iterations, 1);
    X_f = zeros(total_iterations, 1);
    X_r = zeros(total_iterations, 1);
    L_D = zeros(total_iterations, 1);
    
    fprintf('Starting analysis with %d total configurations...\n', total_iterations);
    
    for i = 1:num_stations
        x_m = x_mass_range(i);
        fprintf("%.1f Percent Complete\n", (i / num_stations) * 100);
        for j = 1:num_stations
            x_f = x_f_wing_range(j);
            for k = 1:num_stations
                x_r = x_r_wing_range(k);

                if x_r <= x_f + 0.015 % Ensure wings are separated
                    continue;
                end

                % Building specific plane
                c = getConstructionVector();
                c.x_mass = x_m;
                c.x_f_wing = x_f;
                c.x_r_wing = x_r;
                
                % THEN, create p from the fully updated c
                p = getParameterVector(c); 
                
                % Finding best efficiency
                [L_D_best, ~] = stabilize_computeEfficiency(c);
                
                % Store the results - p.x_cg is now correct for this iteration
                results_count = results_count + 1;
                X_cg(results_count) = p.x_cg;
                X_f(results_count) = x_f;
                X_r(results_count) = x_r;
                L_D(results_count) = L_D_best;
            end
        end
    end

    % Trim unused preallocated space
    X_cg = X_cg(1:results_count);
    X_f = X_f(1:results_count);
    X_r = X_r(1:results_count);
    L_D = L_D(1:results_count);
    fprintf("Ended with %d configurations tested.\n", length(L_D));
    save('myArrays.mat', 'X_cg', 'X_f', 'X_r', 'L_D');
    
    top_5_percent_threshold = prctile(L_D, 95); % Use 95 to get the top 5%

    % Logical mask for values in the top 5%
    mask = (L_D >= top_5_percent_threshold);
    % 
    % % Current values
    % c = getConstructionVector();
    % p = getParameterVector(c);
    
    % Plot the plane
    figure('Position', [100, 100, 1400, 600]); % Made figure wider for subplots
    
    % --- Subplot 1: All Data ---
    subplot(1, 2, 1);
    hold on;
    % Use a lighter alpha to see the distribution
    scatter3(X_cg*1000, X_f*1000, X_r*1000, 36, L_D, 'filled', 'MarkerFaceAlpha', 0.2); 
    hold off;
    xlabel('x_{cg} (mm)');
    ylabel('x_{front wing} (mm)');
    zlabel('x_{rear wing} (mm)');
    title('Efficiency Configurations Spectrum');
    colormap(gca, parula); % Apply colormap to current axes (gca)
    cb = colorbar;
    ylabel(cb, 'L/D Ratio');
    grid on;
    view(135, 30);
    axis tight;

    subplot(1, 2, 2);
    hold on;
    % Use scatter3 for a 3D plot
    scatter3(X_cg(mask)*1000, X_f(mask)*1000, X_r(mask)*1000, 36, L_D(mask), "filled");
    hold off;
    xlabel('x_{cg} (mm)');
    ylabel('x_{front wing} (mm)');
    zlabel('x_{rear wing} (mm)');
    title('Top 5% Efficiency Configurations');
    colormap(gca, hot); % Apply a different colormap to see the difference
    cb = colorbar;
    ylabel(cb, 'L/D Ratio');
    grid on;
    view(135, 30);
    axis tight;

end