function stabilize_staticMarginRange()
    function [static_margin, x_cg] = computeMargin(IN_x_mass, IN_x_f_wing, IN_x_r_wing)    
        % --- Component Parameters ---
        % Nose
        c.m_nose = 0.004;
        c.x_nose = 0;
    
        % Mass
        c.m_mass = 0.009;
        c.x_mass = IN_x_mass;
    
        % Front wing
        c.m_f_wing = 0.004;
        c.x_f_wing = IN_x_f_wing;%0.085;
    
        % Back wing
        c.m_r_wing = 0.006;            
        c.x_r_wing = IN_x_r_wing;%0.140;
    
        % Rudder
        c.m_rudder = 0.002;            
        c.x_rudder = 0.230;
        
        % CF rod
        c.m_rod = 0.005;
        c.L_rod = 0.230;
        c.x_rod = c.L_rod / 2;
    
        % --- Mass ---
        [~, x_cg, ~] = tandemWings_Mech(c);
    
    
        % --- Aerodynamic Parameters ---    
        % Front Wing
        S_f = 0.007;        % Wing area (m^2)
        x_f_wing = c.x_f_wing;
        AR_f = 11.2;
    
        % Rear Wing
        S_r = 0.0084;        % Rear wing is larger
        x_r_wing = c.x_r_wing;
        AR_r = 9.33;


        % --- Calculations ---
        % Center of pressure
        S_tot = S_f + S_r;
        x_cp = (S_f * x_f_wing + S_r * x_r_wing) / S_tot;
        
        % Mean aerodynamic chord
        c_bar = (sqrt(S_f / AR_f) * S_f + sqrt(S_r / AR_r) * S_r) / S_tot;

        % static margin
        static_margin = 100 * (x_cp - x_cg) / c_bar;
    end

    num_stations = 60;

    x_mass_range = linspace(0.00, 0.140, num_stations);
    x_f_wing_range = linspace(0.030, 0.120, num_stations);
    x_r_wing_range = linspace(0.050, 0.170, num_stations);

    % Preallocate arrays to hold valid points
    X_cg = [];
    X_f = [];
    X_r = [];
    SM_values = [];

    for i = 1:num_stations
        x_mass = x_mass_range(i);
        for j = 1:num_stations
            x_f = x_f_wing_range(j);
            for k = 1:num_stations
                x_r = x_r_wing_range(k);

                % Only continue if rear wing is behind front wing
                if x_r <= x_f
                    continue;
                end

                % Compute static margin
                [stat_mar, cg] = computeMargin(x_mass, x_f, x_r);
                    X_cg(end+1) = cg;
                    X_f(end+1) = x_f;
                    X_r(end+1) = x_r;
                    SM_values(end+1) = stat_mar;
            end
        end
    end

    % Finding plane with acceptable static margins

    % Flatten everything into column vectors
    X = X_cg(:); Y = X_f(:); Z = X_r(:); S = SM_values(:);
    
    % Define your target static margin and tolerance
    target_margin = 10;
    tolerance = 2;
    
    % Logical mask for values near the target
    mask = abs(S - target_margin) < tolerance;

    % test = a * 0.085 + b *0.085 + c_;         
    % fprintf("Expect rear wing placement: %.4f", test);
    
    % Plot the plane
    figure('Position', [100, 100, 1000, 700]);
    hold on;
    scatter3(X_cg, X_f, X_r, 36, SM_values, 'filled', MarkerFaceAlpha=0.2);
    scatter3(X_cg(mask), X_f(mask), X_r(mask), 60, 'red', 'x', 'LineWidth', 1.5);
    hold off;
    xlabel('x_{cg} (mm)');
    ylabel('x_{front wing} (mm)');
    zlabel('x_{rear wing} (mm)');
    title('Static Margin vs Component Positions');
    colormap parula;
    colorbar;
    grid on;
    view(135, 30);
end