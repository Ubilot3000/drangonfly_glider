% function alpha_trim = findTrimCondition(v, p)
%     % 
% 
%     % --- Step 1: Find the α where Moment is Zero (the Constraint) ---
%     % Define the anonymous function for the solver. Its job is to find the
%     % pitching moment for any given alpha.
%     moment_function = @(alpha) calculateMoment(p, v, alpha);
% 
%     initial_alpha_guess = deg2rad(3);
%     options = optimset('Display', 'off');
% 
%     try
%         % Find the single alpha that makes the total moment zero.
%         [alpha_moment_trim, ~, exitflag] = fzero(moment_function, initial_alpha_guess, options);
% 
%         if exitflag ~= 1
%             alpha_trim = NaN; % Solver failed to find a moment trim.
%             return;
%         end
%     catch
%         alpha_trim = NaN; % Solver threw an error.
%         return;
%     end
% 
%     % --- Step 2: Verify this is a valid flight state (the Consequence) ---
%     % Now that we have the ONLY alpha the plane can fly at, calculate the
%     % drag force it produces at that alpha.
%     [~, D_total] = calculateAeroForces(alpha_moment_trim, v, p);
% 
%     % If drag is greater than or equal to the glider's weight, it is
%     % impossible to maintain this speed in a steady glide. This is not a
%     % valid trim condition.
%     if D_total >= p.m_total * p.g
%         alpha_trim = NaN;
%     else
%         % If Drag < Weight, a valid glide angle exists. This is our true
%         % and only trim alpha for this speed.
%         alpha_trim = alpha_moment_trim;
%     end
% end

function [alpha_trim, gamma_trim] = findTrimCondition(v, p)
    % Solve M(alpha)=0 and L(alpha)=W*cos(gamma) simultaneously
    %
    W = p.m_total * p.g;

    % fsolve options
    opts = optimoptions('fsolve', ...
        'Display','off', ...
        'FunctionTolerance',1e-8, ...
        'StepTolerance',1e-8);

    % initial guess [alpha; gamma]
    x0 = [deg2rad(3); 0];

    % call fsolve on the nested function
    [xsol, ~, exitflag] = fsolve(@trimEquations, x0, opts);

    if exitflag > 0
        alpha_trim = xsol(1);
        gamma_trim = xsol(2);
        success     = true;
    else
        alpha_trim = NaN;
        gamma_trim = NaN;
        success     = false;
    end

    return;


    %--------------------------------------------------------------------
    function F = trimEquations(x)
        % x(1) = alpha, x(2) = gamma
        % 1) moment eqn:
        M = calculateMoment(p, v, x(1));
        % 2) lift eqn from calculateAeroForces:
        [L, ~] = calculateAeroForces(x(1), v, p);
        % residuals:
        F = [ M; L - W*cos(x(2)) ];
    end
end
