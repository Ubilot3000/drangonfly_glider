function designSpace_compareEfficiencyDampening()
%--- Targets & Figure Setup ---
    x_cg_targets = [0.110, 0.130, 0.150];
    num_cgs      = numel(x_cg_targets);
    figure('Position',[100 100 1600 800]);

    for k = 1:num_cgs
        x_cg_target = x_cg_targets(k);
        fprintf("Analyzing C.G. = %.3f m\n", x_cg_target);

        % Base geometry
        c = getConstructionVector();
        c.x_r_wing = 0.210;

        % Prealloc
        num_stations  = 25;
        x_f_vals      = linspace(0.030,0.120, num_stations);
        decalage_vals = deg2rad(linspace(-5,2, num_stations));
        L_Ds = zeros(num_stations);
        DRs  = zeros(num_stations);

        for i = 1:num_stations
            c.x_f_wing = x_f_vals(i);
            for j = 1:num_stations
                decalage     = decalage_vals(j);
                c.x_mass     = xMassForTargetCG(x_cg_target, c);
                c.i_r        = c.i_f + decalage;

                % Efficiency
                p = getParameterVector(c);
                p.i_r = p.i_f + decalage;
                [L_D_best, ~] = stabilize_computeEfficiency(p);
                if L_D_best > 0
                    L_Ds(j,i) = L_D_best;
                else
                    L_Ds(j,i) = NaN;
                end


                % Stability via linearization
                try
                    [~, ~, drs] = findEigenvalues(p);
                    if isempty(drs)
                        DRs(j,i) = NaN;
                    elseif min(drs) < 0
                        DRs(j,i) = NaN;
                    else
                        DRs(j, i) = min(drs);
                    end
                catch
                    DRs(j,i) = NaN;
                end
            end
        end

        L_rod = c.L_rod;

        % Top: L/D
        subplot(2,num_cgs,k);
        contourf(x_f_vals / L_rod * 100, rad2deg(decalage_vals), L_Ds,15);
        colorbar; title(sprintf('L/D @ CG = %.3f%%',x_cg_target/L_rod*100));
        xlabel('x_f (% of Fuselage)'); ylabel('Decalage (°)'); grid on;

        % Bottom: Damping Ratio
        subplot(2,num_cgs,num_cgs+k);
        contourf(x_f_vals / L_rod * 100, rad2deg(decalage_vals), DRs,15);
        colorbar; title(sprintf('Min Damping Ratio (%%) @ CG = %.3f%%',x_cg_target/L_rod*100));
        xlabel('x_f (% of Fuselage)'); ylabel('Decalage (°)'); grid on;
    end
end
