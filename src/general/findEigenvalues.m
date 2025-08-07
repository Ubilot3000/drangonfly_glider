function [eigens, omega_nats, damping_ratios] = findEigenvalues(p, v0)
    % Finds eigenvalues based off simplified textbook example linearized
    % system. Uses the isolated and decoupled Phugoid and Short Period
    % matrixes.
    % 
    % First finds trim conditon at a given velocity and then uses the
    % construction information to determine all other quantities.
    % 
    % Output:
    %   eigens
    %   omega_nats
    %   damping_ratios
   
    % Finding alpha and forces at trim condition
    [alpha_trim, ~] = findTrimCondition(v0, p);
    
    % --- Finding Cl/Cd ---
    % Front Wing
    alpha_f = alpha_trim + p.i_f;

    Cl_f = p.Cl_a * alpha_f;
    Cd_f = p.Cd_0 + p.k_f * Cl_f^2;
    
    % Downwash
    epsilon = (2 * Cl_f) / (pi * p.AR_f) * (1 / (1 + p.k_wings * (p.delta_wings / p.span_f)^2));
    
    % Rear Wing
    alpha_r_eff = alpha_trim + p.i_r - epsilon;
    Cl_r = p.Cl_a * alpha_r_eff;

    Cd_r = p.Cd_0 + p.k_r * Cl_r^2;
    Cl_total  = (Cl_f * p.S_f + Cl_r * p.S_r) / p.S_tot;
    Cd_total  = (Cd_f * p.S_f + Cd_r * p.S_r) / p.S_tot;

    L_D = Cl_total / Cd_total;

    % --- Phugoid simplified matrix ---
    Au = -2 * p.g / (L_D * v0);
    Ag = -   p.g /   v0;
    Bu =  2 * p.g /   v0;
    A_ph = [ Au,  Ag;  Bu,   0 ];

    % --- Pitching coefficient ---
    % Set flight deviation
    alpha_delta = deg2rad(0.1);

    % Calculate Cm_alpha
    M_1 = calculateMoment(p, v0, alpha_trim);
    M_2 = calculateMoment(p, v0, alpha_trim + alpha_delta);

    % Moment coefficients
    Q = 0.5 * p.rho * v0^2; % Dynamic pressure
    Cm_1 = M_1 / (Q * p.S_tot * p.c_bar);
    Cm_2 = M_2 / (Q * p.S_tot * p.c_bar);

    % Moment alpha coefficient
    Cm_alpha = (Cm_2 - Cm_1) / alpha_delta;

    % Moment derivative coefficient
    Cm_q = - Cm_alpha * (p.c_bar / (2 * v0));

    % Short period simplified matrix
    Ca = - (0.5 * p.rho * p.S_tot * v0^2 * p.Cl_a) / (p.m_total * v0);
    Cq = 1;
    Da =   (0.5 * p.rho * p.S_tot * p.c_bar * v0^2 * Cm_alpha) / p.I_total;
    Dq =   (0.5 * p.rho * p.S_tot * p.c_bar^2 * v0   * Cm_q) / p.I_total;
    A_sp = [ Ca, Cq; Da, Dq ];

    % --- Finding eigenvalues ---
    phug_eigs = eig(A_ph);
    sp_eigs   = eig(A_sp);

    % Gathering all data
    all_eigs = [ phug_eigs; sp_eigs ];
    wn_all   = abs(all_eigs);
    zt_all   = -real(all_eigs) ./ wn_all;

    % sort so that index 1 = lower‐freq Phugoid, 2 = its pair,
    % and 3&4 = Short‐period pair (higher freq)
    [~,ord]      = sort(wn_all);
    eigens       = all_eigs(ord);
    omega_nats   = wn_all(ord);
    damping_ratios = zt_all(ord);
    
end
