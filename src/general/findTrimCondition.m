function alpha_trim = findTrimCondition(v, p)
    % 

    % --- Step 1: Find the α where Moment is Zero (the Constraint) ---
    % Define the anonymous function for the solver. Its job is to find the
    % pitching moment for any given alpha.
    moment_function = @(alpha) calculateMoment(p, v, alpha);
    
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
    [~, D_total] = calculateAeroForces(alpha_moment_trim, v, p);

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