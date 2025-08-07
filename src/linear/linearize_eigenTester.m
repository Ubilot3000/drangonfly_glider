function linearize_eigenTester()
    x_cg_min = 0; 
    x_cg_max = 0.230; 
    num_stations = 100;
    all_x_cg = linspace(x_cg_min, x_cg_max, num_stations);
    
    % Preallocating arrays
    eigs1 = NaN(4, num_stations);
    eigs2 = NaN(4, num_stations);
    eigs3 = NaN(4, num_stations);

    % Model construction
    c = getConstructionVector();
    c.x_r_wing = 0.200;
    v0 = 4;    % a realistic trim speed

    for i = 1:num_stations
        target_x_cg = all_x_cg(i);
        c.x_mass = xMassForTargetCG(target_x_cg, c);
        p = getParameterVector(c);

        %--- Call the trimmed‐based eigen solver --- 
        [e1, trim, ok] = findEigenvaluesSolve(p, v0);
        if ~ok
            fprintf("Not Okay, skipping: %.3fmm\n", target_x_cg * 100);
            continue;
        end
        v_trim = trim.V;

        fprintf("CG: %.3fmm, Equilibrium Velocity: %.1fm/s\n", target_x_cg * 100, v_trim)

        [e2, ~, ~] = findEigenvalues(p, v_trim);
        [e3, ~]    = findEigenvaluesRobust(p, v_trim);

        % --- SORT by |Im(lambda)| so phugoid (low-freq) are 1–2, SP are 3–4 ---
        [~, idx1] = sort(abs(imag(e1)));
        e1 = e1(idx1);
        [~, idx2] = sort(abs(imag(e2)));
        e2 = e2(idx2);
        [~, idx3] = sort(abs(imag(e3)));
        e3 = e3(idx3);

        % disp(e1);
        % disp(real(e1(1)));
        if real(e1(1)) >= 0
            fprintf("CG: %.3fmm unstable with: %d\n", target_x_cg * 100, real(e1(1)));
        end

        % store as columns
        eigs1(:,i) = e1;
        eigs2(:,i) = e2;
        eigs3(:,i) = e3;
    end


    figure("Position", [100, 100, 1200, 600]);

    % --- Phugoid Subplot ---
    subplot(1, 2, 1);
    hold on; % Turn hold on once at the beginning
    
    % -- Plot the real data --
    % Phugoid modes are the first two rows (complex conjugate pair)
    scatter(real(eigs1(1:2,:)), imag(eigs1(1:2,:)), 36, "filled", "MarkerFaceColor", "#0072BE");
    scatter(real(eigs2(1:2,:)), imag(eigs2(1:2,:)), 36, "filled", "MarkerFaceColor", "#DA5319");
    scatter(real(eigs3(1:2,:)), imag(eigs3(1:2,:)), 36, "filled", "MarkerFaceColor", "#EEB220");

    % -- Plot the trajectory arrows --
    % We only need to trace the path of one of the conjugate pairs (e.g., the one with positive imag)
    x1 = real(eigs1(1,:)); 
    y1 = imag(eigs1(1,:));
    x2 = real(eigs2(1,:)); 
    y2 = imag(eigs2(1,:));
    x3 = real(eigs3(1,:)); 
    y3 = imag(eigs3(1,:));
    quiver(x1(1:end-1), y1(1:end-1), diff(x1), diff(y1), 0, 'Color', "#0072BE", 'MaxHeadSize',0.25, 'AutoScale','off');
    quiver(x2(1:end-1), y2(1:end-1), diff(x2), diff(y2), 0, 'Color', "#DA5319", 'MaxHeadSize',0.25, 'AutoScale','off');
    quiver(x3(1:end-1), y3(1:end-1), diff(x3), diff(y3), 0, 'Color', "#EEB220", 'MaxHeadSize',0.25, 'AutoScale','off');
    
    % -- Create invisible proxy graphics for the legend --
    % We plot them at NaN to make them invisible on the chart, but give them a DisplayName.
    h1 = scatter(NaN, NaN, 36, "filled", "MarkerFaceColor", "#0072BE", "DisplayName", "fSolve");
    h2 = scatter(NaN, NaN, 36, "filled", "MarkerFaceColor", "#DA5319", "DisplayName", "2x2 Simplified");
    h3 = scatter(NaN, NaN, 36, "filled", "MarkerFaceColor", "#EEB220", "DisplayName", "4x4 Matrix Simplified");

    hold off;
    
    % -- Finalize the plot --
    axis equal;
    xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
    title('Phugoid vs CG');
    sgrid; % Call sgrid AFTER plotting all data
    legend([h1, h2, h3], "Location", "best"); % Create the legend from the proxies

    % --- Short-Period Subplot ---
    subplot(1, 2, 2);
    hold on;

    % -- Plot the real data --
    % Short-period modes are the last two rows (complex conjugate pair)
    scatter(real(eigs1(3:4,:)), imag(eigs1(3:4,:)), 36, "filled", "MarkerFaceColor", "#0072BE");
    scatter(real(eigs2(3:4,:)), imag(eigs2(3:4,:)), 36, "filled", "MarkerFaceColor", "#DA5319");
    scatter(real(eigs3(3:4,:)), imag(eigs3(3:4,:)), 36, "filled", "MarkerFaceColor", "#EEB220");
    
    % -- Plot the trajectory arrows --
    % Travel data
    a1 = real(eigs1(3,:)); 
    b1 = imag(eigs1(3,:));
    a2 = real(eigs2(3,:)); 
    b2 = imag(eigs2(3,:));
    a3 = real(eigs3(3,:)); 
    b3 = imag(eigs3(3,:));

    quiver(a1(1:end-1), b1(1:end-1), diff(a1), diff(b1), 0, 'Color', "#0072BE", 'MaxHeadSize',0.25, 'AutoScale','off');
    quiver(a2(1:end-1), b2(1:end-1), diff(a2), diff(b2), 0, 'Color', "#DA5319", 'MaxHeadSize',0.25, 'AutoScale','off');
    quiver(a3(1:end-1), b3(1:end-1), diff(a3), diff(b3), 0, 'Color', "#EEB220", 'MaxHeadSize',0.25, 'AutoScale','off');

    % --- Fake legend data ---
    g1 = scatter(NaN, NaN, 36, "filled", "MarkerFaceColor", "#0072BE", "DisplayName", "fSolve");
    g2 = scatter(NaN, NaN, 36, "filled", "MarkerFaceColor", "#DA5319", "DisplayName", "2x2 Simplified");
    g3 = scatter(NaN, NaN, 36, "filled", "MarkerFaceColor", "#EEB220", "DisplayName", "4x4 Matrix Simplified");
    hold off;
    % -- Finalize the plot --
    axis equal;
    xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
    title('Short Period vs CG');
    sgrid; 
    legend([g1, g2, g3], "Location", "best");


end
