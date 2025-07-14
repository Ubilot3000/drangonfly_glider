function analyzePitchAuthority()
    % --- Environmental Parameters ---
    p.g = 9.81;
    p.rho = 1.225;

    % --- Component Parameters ---
    % Nose
    c.m_nose = 0.004;
    c.x_nose = 0;

    % Mass
    c.m_mass = 0.0035;
    % c.x_mass = 0; % Not needed for changing mass

    % Front wing
    c.m_f_wing = 0.004;
    c.x_f_wing = 0.085;

    % Back wing
    c.m_r_wing = 0.006;            
    c.x_r_wing = 0.140;

    % Rudder
    c.m_rudder = 0.002;            
    c.x_rudder = 0.230;
    
    % CF rod
    c.m_rod = 0.005;
    c.L_rod = 0.230;
    c.x_rod = c.L_rod / 2;

    % --- Aerodynamic Parameters ---
    % Aerodynamic Coefficients
    p.Cl_a = 2 * pi;     % Lift curve slope (per radian)
    p.Cd_0 = 0.02;       % Parasitic drag
    oswald_e = 0.9;     % Oswald efficiency factor 

    % Front Wing
    p.S_f = 0.007;        % Wing area (m^2)
    p.x_f_wing = c.x_f_wing;
    p.i_f = deg2rad(6);  % Incidence angle (rad)
    p.AR_f = 11.2;
    p.k_f = 1 / (pi * oswald_e * p.AR_f);

    % Rear Wing
    p.S_r = 0.0084;        % Rear wing is larger
    p.x_r_wing = c.x_r_wing;
    p.i_r = deg2rad(3);
    p.AR_r = 9.33;
    p.k_r = 1 / (pi * oswald_e * p.AR_r);

    % Pitching Moment 
    p.c_ref = (p.S_f / p.AR_f^0.5 + p.S_r / p.AR_r^0.5) / 2; % A weighted average chord
    p.S_total = p.S_f + p.S_r;

    % --- Flight Parameters ---
    % Initial conditions
    V_test = 4;
    alpha_test = deg2rad(0);
    num_stations = 50;
    x_mass_min = 0.00;
    x_mass_max = 0.04;
    mass_positions = linspace(x_mass_min, x_mass_max, num_stations);

    % Data Arrays
    cg_positions = zeros(1, num_stations);
    moment = zeros(1, num_stations);
    moment_coefficients = zeros(1, num_stations);
    
    for i = 1:num_stations
        c.x_mass = mass_positions(i);

        [m_total, x_cg, ~] = tandemWings_Mech(c);


        p.m_total = m_total
        p.x_cg = x_cg;
        cg_positions(i) = x_cg;
        
        M_total = calculateMoment(p, V_test, alpha_test);
        moment(i) = M_total;

        Q = 0.5 * p.rho * V_test^2; % Dynamic pressure
        Cm = M_total / (Q * p.S_total * p.c_ref);
        moment_coefficients(i) = Cm;       
    end

    figure("Position", [100, 100, 900, 700]);
    sgtitle('Pitch Control Authority Analysis', 'FontSize', 16);
    subplot(2, 1, 1);
    hold on;
    plot(mass_positions * 1000, moment_coefficients, 'b-', 'LineWidth', 2, "Color", "#0072BE");
    yline(0, 'k--');
    hold off;
    grid on;
    title('Moment Coefficient caused by Mass Translation');
    xlabel('Movable Mass Position (mm from nose)');
    ylabel('Pitching Moment Coefficient (C_m)');
    


    % Cm / alpha slope
    plot_coefficients = polyfit(mass_positions, moment_coefficients, 1);
    slope = plot_coefficients(1);
    intercept = plot_coefficients(2);
    Cm_function = @(x_mass) slope*x_mass + intercept;

    x_mass_trim = fzero(Cm_function, (-intercept / slope));

    subplot(2, 1, 2);
    hold on;
    plot(mass_positions * 1000, moment, "r-", "LineWidth", 2, "Color", "#DA5319");
    yline(0, "k--");
    hold off;
    grid on;
    title("Moment caused by Mass Translation");
    xlabel("Movable Mass Position (mm from nose)");
    ylabel("Pitching Moment (Nm)")
    
    fprintf('Analysis Complete. Plot generated.\n');
    fprintf('C_m range: [%.3f, %.3f]\n', min(moment_coefficients), max(moment_coefficients));
    fprintf("Plotted slope: %.2f, intercept: %.2f\n", slope, intercept);
    fprintf("%.3fkg Mass Trim Position: %.2fmm\n", c.m_mass, x_mass_trim * 1000);
end