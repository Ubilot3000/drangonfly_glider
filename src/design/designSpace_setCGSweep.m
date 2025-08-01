function designSpace_setCGSweep()
    x_cg_target = [0.080, 0.095, 0.110, 0.130];

    figure('Position',[200 200  1400 500]);

    for k = 1:length(x_cg_target)
       

        % Construct plane
        c = getConstructionVector();
    
        % Predefine sweep parameters
        num_stations = 25;
        x_f_vals = linspace(0.030, 0.120, num_stations);
        x_r_vals= linspace(0.080, 0.170, num_stations);

        % Preallocate arrays.
        results_count = 0;
        total_iterations = num_stations^3;
        X_f = zeros(total_iterations, 1);
        X_r = zeros(total_iterations, 1);
        L_D = zeros(total_iterations, 1);
    
        for i = 1:num_stations
            x_f = x_f_vals(i);
            for j = 1:num_stations
                x_r = x_r_vals(j);
    
                % Checking valid wing configuration
                if x_r < x_f + 0.015
                    continue;
                end
    
                % Updating construction
                c_mod = c;
                c_mod.x_f_wing = x_f;
                c_mod.x_r_wing = x_r;
                c_mod.x_mass = xMassForTargetCG(x_cg_target(k), c_mod);
                p = getParameterVector(c_mod);
    
                [L_D_best, ~] = stabilize_computeEfficiency(p);
    
                % Save results
                results_count = results_count + 1;
                X_f(results_count) = x_f;
                X_r(results_count) = x_r;
                L_D(results_count) = L_D_best;     
            end
        end
    
        % Trimming values
        Xf = X_f(1:results_count);
        Xr = X_r(1:results_count);
        Z = L_D(1:results_count);
    
        % ---- OPTION B: gridded surface ----
        % define a fine mesh covering the convex hull of (Xf,Xr)
        xi = linspace(min(Xf), max(Xf), 50);
        yi = linspace(min(Xr), max(Xr), 50);
        [Xi, Yi] = meshgrid(xi, yi);
    
        % interpolate the scattered data onto the mesh
        Zi = griddata(Xf, Xr, Z, Xi, Yi, 'natural');
    
        [maxZ, linIdx] = max(Zi(:));               % global max of the grid
        [row, col] = ind2sub(size(Zi), linIdx);    % convert to (row,col)
        xMax = Xi(row, col);
        yMax = Yi(row, col);
    
        subplot(1, length(x_cg_target), k);
        hold on;
        contourf(Xi, Yi, Zi, 10);
        plot(xMax, yMax, 'rx', 'MarkerSize',12, 'LineWidth',2);
        text(xMax, yMax, sprintf('  Max = %.2f', maxZ),'Color','red', 'FontWeight','bold', 'VerticalAlignment','bottom');
        hold off;
        shading flat;
        clim([min(Z(~isnan(Z))), max(Z(~isnan(Z)))]);
        colorbar;
        xlabel('$x_{\rm fore}$ (m)','Interpreter','latex');
        ylabel('$x_{\rm aft}$ (m)','Interpreter','latex');
        zlabel('$L/D$','Interpreter','latex');
        title('Interpolated $L/D$ Surface','Interpreter','latex');
        grid on;
    
        fprintf("Highest Efficiency: %.2f\n", maxZ);
    end
    
end