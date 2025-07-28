function stabilize_epsilonDeltaWings()
    delta_range = linspace(0.01, 0.3, 200);  % Wing spacing in meters
    epsilon = zeros(size(delta_range));     % Preallocate
    epsilon_old = zeros(size(delta_range));
    epsilon_alt = zeros(size(delta_range));
    epsilon_alt_2 = zeros(size(delta_range));

    % Model parameters
    c = getConstructionVector();
    p = getParameterVector(c);

    alpha_deg = 4; 
    Cl_f = p.Cl_a * deg2rad(alpha_deg);
    k_wings_alt = 12;
    k_wings_alt_2 = 6;
    
    % Epsilon for each delta
    for i = 1:length(delta_range)
        delta = delta_range(i);
        epsilon(i) = (2 * Cl_f) / (pi * p.AR_f) * (1 / (1 + p.k_wings * (delta / p.span_f)^2));
        epsilon_old(i) = (2 * Cl_f) / (pi * p.AR_f);
        epsilon_alt(i) = (2 * Cl_f) / (pi * p.AR_f) * (1 / (1 + k_wings_alt * (delta / p.span_f)^2));
        epsilon_alt_2(i) = (2 * Cl_f) / (pi * p.AR_f) * (1 / (1 + k_wings_alt_2 * (delta / p.span_f)^2));
    end
    
    % --- Plotting ---
    figure;
    hold on;
    plot(delta_range, rad2deg(epsilon_old), "LineWidth", 2);
    plot(delta_range, rad2deg(epsilon_alt_2), "LineWidth", 2);
    plot(delta_range, rad2deg(epsilon), 'LineWidth', 2);
    plot(delta_range, rad2deg(epsilon_alt), "LineWidth", 2);
    hold off;
    xlabel('Wing Spacing \delta (m)');
    ylabel('Downwash Angle \epsilon (degrees)');
    title(sprintf('Downwash Angle vs Wing Spacing at \\alpha = %d^\\circ', alpha_deg));
    legend(["Old Model","K = 6", "K = 8", "K = 12"]);
    grid on;
end