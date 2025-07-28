function stabilize_staticMarginRange()
    num_stations = 35;

    x_mass_range = linspace(0.00, 0.190, num_stations);
    x_f_wing_range = linspace(0.030, 0.120, num_stations);
    x_r_wing_range = linspace(0.050, 0.170, num_stations);

    % Preallocate arrays.
    results_count = 0;
    total_iterations = num_stations^3;
    X_cg = zeros(total_iterations, 1);
    X_f = zeros(total_iterations, 1);
    X_r = zeros(total_iterations, 1);
    SM = zeros(total_iterations, 1);
    
    fprintf('Starting analysis with %d total configurations...\n', total_iterations);
    
    for i = 1:num_stations
        x_m = x_mass_range(i);
        fprintf("%.1f%% complete...\n", i / num_stations * 100);
        for j = 1:num_stations
            x_f = x_f_wing_range(j);
            for k = 1:num_stations
                x_r = x_r_wing_range(k);

                if x_r <= x_f + 0.015 % Ensure wings are separated
                    continue;
                end

                c = getConstructionVector();
                
                c.x_mass = x_m;
                c.x_f_wing = x_f;
                c.x_r_wing = x_r;
                
                [static_margin, ~, x_cg] = stabilize_staticStability(c);
                
                % Store the results
                results_count = results_count + 1;
                X_cg(results_count) = x_cg;
                X_f(results_count) = x_f;
                X_r(results_count) = x_r;
                SM(results_count) = static_margin;
            end
        end
    end

    % Trim unused preallocated space
    X_cg = X_cg(1:results_count);
    X_f = X_f(1:results_count);
    X_r = X_r(1:results_count);
    SM = SM(1:results_count);
    fprintf("Ended with %d configurations tested.\n", length(SM));
    
    % Define your target static margin and tolerance
    target_margin = 10;
    tolerance = 2;
    
    % Logical mask for values near the target
    S = SM(:);
    mask = abs(S - target_margin) < tolerance;

    % Current values
    c = getConstructionVector();
    p = getParameterVector(c);
    
    % Plot the plane
    figure('Position', [100, 100, 1000, 700]);
    hold on;
    scatter3(X_cg, X_f, X_r, 36, SM, 'filled', MarkerFaceAlpha=0.1);
    scatter3(X_cg(mask), X_f(mask), X_r(mask), 60, 'red', 'x', 'LineWidth', 1.5);
    scatter3(p.x_cg, c.x_f_wing, c.x_r_wing, 60, "black", "x", "LineWidth", 3);
    h1 = scatter3(NaN, NaN, NaN, 60, 'red', 'x', 'LineWidth', 1.5, 'DisplayName', 'Static Margin \approx 10%');
    h2 = scatter3(NaN, NaN, NaN, 60, 'black', 'x', 'LineWidth', 3, 'DisplayName', 'Current Static Margin');
    hold off;
    xlabel('x_{cg} (mm)');
    ylabel('x_{front wing} (mm)');
    zlabel('x_{rear wing} (mm)');
    title('Static Margin vs Component Positions');
    colormap parula;
    colorbar;
    legend([h1, h2]);
    grid on;
    view(135, 30);
end