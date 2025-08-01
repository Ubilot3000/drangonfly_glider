function [eigens, omega_nats, damping_ratios] = findEigenvalues(p, v0)
    %--- 1) Trim at speed v0 to get alpha_trim, gamma0, theta0, q0 ----
    alpha_trim = findTrimCondition(v0, p);
    if isnan(alpha_trim)
        eigens = []; omega_nats = []; damping_ratios = [];
        return;
    end
    [~, D_total] = calculateAeroForces(alpha_trim, v0, p);
    gamma0 = asin(-D_total/(p.m_total*p.g));
    theta0 = gamma0 + alpha_trim;
    q0     = 0;
    z_eq   = [v0; gamma0; theta0; q0];

    %--- 2) Build Jacobian A via central differences ---------------
    ode_func = @(z) linearize_ODE(0, z, p);
    n = numel(z_eq);
    A = zeros(n);
    dz = 1e-6;
    for i = 1:n
        zp = z_eq; zm = z_eq;
        zp(i) = zp(i) + dz;
        zm(i) = zm(i) - dz;
        A(:,i) = (ode_func(zp) - ode_func(zm))/(2*dz);
    end

    %--- 3) Full eigendecomposition ------------------------------
    [V_all, D_all] = eig(A);
    lam_all        = diag(D_all);           % 4×1 complex eigenvalues
    omega_all      = abs(lam_all);          % nat freqs
    zeta_all       = -real(lam_all)./omega_all;  % damping ratios

    %--- 4) Pick the two oscillatory, positive‐Im modes ----------
    osc = find(imag(lam_all)~=0 & imag(lam_all)>0);
    if numel(osc) < 2
        % not enough oscillatory modes
        eigens = []; omega_nats = []; damping_ratios = [];
        return;
    end

    % isolate their eigenvectors
    Vpos = V_all(:, osc);     % 4×2
    lam2 = lam_all(osc);      % 2×1
    w2   = omega_all(osc);
    z2   = zeta_all(osc);

    %--- 5) Score “short‐periodness” by |θ|+|q| participation -------
    % state ordering: z = [ V; gamma; theta; q ]
    sp_scores = abs(Vpos(3,:)) + abs(Vpos(4,:));
    % whichever score is larger is the short‐period mode
    [~, i_sp] = max(sp_scores);
    i_ph     = 3 - i_sp;  % the other index

    %--- 6) Assemble outputs in the order [short; phugoid] ------
    eigens       = [ lam2(i_sp);    lam2(i_ph)   ];
    omega_nats   = [ w2( i_sp);      w2( i_ph)    ];
    damping_ratios = [ z2(i_sp);     z2(i_ph)     ];
end
