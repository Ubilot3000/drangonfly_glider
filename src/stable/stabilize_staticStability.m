function [static_margin, neutral_point, x_cg] = stabilize_StaticStability(c)
    % --- Variable Setup ---
    % Parameters
    p = getParameterVector(c);
    x_cg = p.x_cg;

    % Initial Conditions
    v0  = 3;
    alpha0 = deg2rad(2);
    initial_guess_np = c.x_f_wing + 0.05; % neutral point behind front wing

    % Utility functions
    Cma_func = @(x_cg) stabilize_computeCma(x_cg, c, v0, alpha0);
    
    % Finding neutral point
    options = optimset('Display', 'off'); 
    [neutral_point, ~, exitflag] = fzero(Cma_func, initial_guess_np, options);
    
    if exitflag ~= 1
        error('Solver could not find the Neutral Point. Check configuration.');
    end

    % Finding static margin
    static_margin = (neutral_point - x_cg) / p.c_bar * 100;
end