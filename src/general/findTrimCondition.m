function alpha_trim = findTrimCondition(v, p)
    % This function finds the true, flyable trim Angle of Attack for a given
    % design 'p' at a specific airspeed 'v'.
    %
    % DEEPER EXPLANATION: A true trim state requires both force balance (L & D
    % balance W) AND moment balance (M=0). For an uncontrolled, stable
    % aircraft, the moment equation acts as a CONSTRAINT. There is only one
    % angle of attack (alpha) where the moments will naturally balance for
    % a given speed. This function first finds that single, mandated alpha,
    % and then checks if the forces produced at that alpha allow for a
    % physically possible steady glide.

    % --- Step 1: Find the α where Moment is Zero (the Constraint) ---
    % Define the anonymous function for the solver. Its job is to find the
    % pitching moment for any given alpha.
    moment_function = @(alpha) calculateMoment(alpha, v, p);
    
    initial_alpha_guess = deg2rad(3);
    options = optimset('Display', 'off');

    try
        % Find the single alpha that makes the total moment zero.
        [alpha_moment_trim, ~, exitflag] = fzero(moment_function, initial_alpha_guess, options);

        if exitflag ~= 1
            alpha_trim = NaN; % Solver failed to find a moment trim.
            return;
        end
    catch
        alpha_trim = NaN; % Solver threw an error.
        return;
    end

    % --- Step 2: Verify this is a valid flight state (the Consequence) ---
    % Now that we have the ONLY alpha the plane can fly at, calculate the
    % drag force it produces at that alpha.
    D_total = calculateDrag(alpha_moment_trim, v, p);

    % If drag is greater than or equal to the glider's weight, it is
    % impossible to maintain this speed in a steady glide. This is not a
    % valid trim condition.
    if D_total >= p.m_total * p.g
        alpha_trim = NaN;
    else
        % If Drag < Weight, a valid glide angle exists. This is our true
        % and only trim alpha for this speed.
        alpha_trim = alpha_moment_trim;
    end
end

% =========================================================================
% === LOCAL HELPER FUNCTIONS (Visible only inside this file) ============
% =========================================================================

function M_total = calculateMoment(alpha, V, p)
    % This local function calculates the total pitching moment.
    
    % --- Aerodynamic Calculation Block ---
    alpha_f = alpha + p.i_f;
    Cl_f = p.Cl_a * alpha_f;
    L_f = 0.5 * p.rho * V^2 * p.S_f * Cl_f;
    D_f = 0.5 * p.rho * V^2 * p.S_f * (p.Cd_0 + p.k_f * Cl_f^2);
    
    epsilon = (2 * Cl_f) / (pi * p.AR_f);
    
    alpha_r_eff = alpha + p.i_r - epsilon;
    Cl_r = p.Cl_a * alpha_r_eff;
    L_r = 0.5 * p.rho * V^2 * p.S_r * Cl_r;
    D_r = 0.5 * p.rho * V^2 * p.S_r * (p.Cd_0 + p.k_r * Cl_r^2);
    
    % Forces normal to the fuselage chord are used for moment calculation
    N_f = L_f * cos(alpha) + D_f * sin(alpha);
    N_r = L_r * cos(alpha) + D_r * sin(alpha);
    
    % Sum moments about the center of gravity
    M_total = (N_f * (p.x_cg - p.x_f_wing)) + (N_r * (p.x_cg - p.x_r_wing));
end

function D_total = calculateDrag(alpha, V, p)
    % This local function calculates only the total drag force.
    
    % --- Aerodynamic Calculation Block ---
    alpha_f = alpha + p.i_f;
    Cl_f = p.Cl_a * alpha_f;
    D_f = 0.5 * p.rho * V^2 * p.S_f * (p.Cd_0 + p.k_f * Cl_f^2);
    
    epsilon = (2 * Cl_f) / (pi * p.AR_f);
    
    alpha_r_eff = alpha + p.i_r - epsilon;
    Cl_r = p.Cl_a * alpha_r_eff;
    D_r = 0.5 * p.rho * V^2 * p.S_r * (p.Cd_0 + p.k_r * Cl_r^2);
    
    D_total = D_f + D_r;
end