function linearize_massPosSweeper()
    % ===========================
    % --- MAIN FUNCTION START ---
    % ===========================
    graph1 = false;
    graph2 = true;

    x_cg_min = 0; 
    x_cg_max = 0.230; 
    num_stations = 50;
    all_x_cg = linspace(x_cg_min, x_cg_max, num_stations);
    
    
    
    % Preallocating arrays
    omega_short    = NaN(1,num_stations);
    omega_phugoid  = NaN(1,num_stations);
    damp_short     = NaN(1,num_stations);
    damp_phugoid   = NaN(1,num_stations);
    % eigens_short   = NaN(1,num_stations);
    % eigens_phugoid = NaN(1,num_stations);
    eig_all = NaN(4, num_stations);

    % Model construction
    c = getConstructionVector();
    c.x_r_wing = 0.200;
    v0 = 4;    % a realistic trim speed


    for j = 1:num_stations
        target_x_cg = all_x_cg(j);
        c.x_mass = xMassForTargetCG(target_x_cg, c);
        p = getParameterVector(c);

        %--- Call the trimmed‐based eigen solver --- 
        [eigens, omega_nats, damping_ratios] = findEigenvalues(p, v0);

        if numel(eigens)==4
            eig_all(:,j) = eigens;
        else
            fprintf("Skipping idx %d: got %d total modes\n", j, numel(eigens));
        end
    end


    % if graph1
    %     % Find instability boundary (where phugoid goes unstable)
    %     instability_idx = find(damp_phugoid < 0, 1, 'first');
    %     x_unstable_start = NaN;
    %     if ~isempty(instability_idx)
    %         x_unstable_start = all_x_cg(instability_idx);
    %     end
    % 
    %     % --- Plotting ---
    %     figure('Position', [100, 100, 1000, 400]);
    %     sgtitle("Flight Stability Analyses");
    % 
    %     % Frequency Plot
    %     subplot(1, 2, 1); 
    %     hold on;
    %     plot(all_x_cg, omega_short, 'LineWidth', 1.5);
    %     plot(all_x_cg, omega_phugoid, 'LineWidth', 1.5);
    %     hold off;
    %     title('Natural Frequency vs Mass Position');
    %     ylabel('\omega_n (rad/s)'); xlabel('x_{mass} (m)');
    %     legend("Short Mode", "Phugoid Mode", 'Location', 'northeast');
    %     grid on;
    % 
    %     % Damping Ratio Plot
    %     subplot(1, 2, 2); 
    %     hold on;
    % 
    %     % Shade unstable phugoid region
    %     if ~isnan(x_unstable_start)
    %         xfill = [x_unstable_start, max(all_x_cg)];
    %         fill([xfill(1), xfill(2), xfill(2), xfill(1)], [-1, -1, 1, 1], [1, 0.8, 0.8], 'EdgeColor', 'none', "FaceAlpha",0.5);
    %     end
    % 
    %     % Plot damping curves
    %     plot(all_x_cg, damp_short, 'LineWidth', 1.5);
    %     plot(all_x_cg, damp_phugoid, 'LineWidth', 1.5);
    % 
    %     title('Damping Ratio vs Mass Position');
    %     ylabel('\zeta'); xlabel('x_{mass} (m)');
    %     legend("Instability Region", "Short Mode", "Phugoid Mode", 'Location', 'southeast');
    %     grid on;
    %     hold off;
    % end

     if graph2
        % --- 1) figure & common setup ---
        figure('Position',[100 100 1000 400]);
        valid = all(~isnan(eig_all),1);      % columns where we got 4 modes
        xs   = all_x_cg(valid)*1000;         % CG positions [mm]
        
        % build a color array, two entries per station (for each conj pair)
        cmap_vals = repmat(xs, 2, 1);        % 2×K
        cmap_vals = cmap_vals(:);            % (2K)×1

        % --- 2) Phugoid subplot (rows 1&2) ---
        subplot(1,2,1); 
        ph = eig_all(1:2, valid);            % 2×K
        x_ph = real(ph(:));
        y_ph = imag(ph(:));
        scatter(x_ph, y_ph, 36, cmap_vals, 'filled');
        hold on;
        % plot(x_ph, y_ph, 'k-', 'LineWidth',1);
        hold off;
        axis equal; 
        % grid on;
        xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
        title('Phugoid Modes');
        sgrid;
        colormap(jet);
        cb = colorbar;
        cb.Label.String = 'CG position [% of Length]';
        clim([x_cg_min / c.L_rod * 100, x_cg_max / c.L_rod *100]);
        sgrid;

        % --- 3) Short‐period subplot (rows 3&4) ---
        subplot(1,2,2); 
        sp = eig_all(3:4, valid);            % 2×K
        x_sp = real(sp(:));
        y_sp = imag(sp(:));
        scatter(x_sp, y_sp, 36, cmap_vals, 'filled');
        hold on;
        % plot(x_sp, y_sp, 'k-', 'LineWidth',1);
        hold off;
        axis equal; 
        % grid on;
        xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
        title('Short‐period Modes');
        colormap(jet);
        cb = colorbar;
        cb.Label.String = 'CG position [% of Length]';
        clim([x_cg_min / c.L_rod * 100, x_cg_max / c.L_rod *100]);
        sgrid;
    end


        % subplot(1, 2, 1);
        % hold on;
        % cmap = jet(length(eigens_short));
        % for k = 1:length(eigens_short)-1
        %   plot(real(eigens_short(k:k+1)), ...
        %        imag(eigens_short(k:k+1)), ...
        %        'LineWidth',2, 'Color',cmap(k,:));
        % end
        % xline(0, '--k');
        % yline(0, '--k');
        % hold off;
        % axis equal;
        % xlabel('Real Part'); 
        % ylabel('Imaginary Part');
        % title('Short Mode Trajectories');
        % sgrid; 
        % 
        % subplot(1, 2, 2);
        % hold on;
        % cmap = jet(length(eigens_phugoid));
        % hold on
        % for k = 1:length(eigens_phugoid)-1
        %     plot(real(eigens_phugoid(k:k+1)), ...
        %         imag(eigens_phugoid(k:k+1)), ...
        %         'LineWidth',2, 'Color',cmap(k,:));
        % end
        % colormap(jet);       % set the colormap
        % c = colorbar;        % draw the colorbar
        % c.Label.String = 'C.G. Positions (mm)';  % optional label
        % clim([x_cg_min * 1000, x_cg_max * 1000]);
        % xline(0, '--k');
        % yline(0, '--k');
        % hold off;
        % axis equal;
        % xlabel('Real Part'); 
        % ylabel('Imaginary Part');
        % title('Phugoid Mode Trajectories');
end