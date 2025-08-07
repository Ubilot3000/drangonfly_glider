function [eigs, trim, ok] = findEigenvaluesSolve(p, v_guess)
    % Finds equilibrium condition of for anglular variables AND velocity.
    %
    % Inputs:
    %   p       struct of aircraft params
    %   v_guess initial guess for airspeed [m/s]
    % Outputs:
    %   eigs    4×1 vector of eigenvalues of A
    %   A       4×4 Jacobian at trim
    %   trim    struct with fields V, alpha, gamma, theta, q
    %   success true if fsolve converged
    
    %--- set up fsolve
    opts = optimoptions('fsolve', ...
        'Display','off', ...
        'FunctionTolerance',1e-8, ...
        'StepTolerance',1e-8...
        );
    
    % unknowns: x = [v_guess; alpha; gamma]
    x0    = [v_guess; deg2rad(3); 0];
    [xsol, ~, exitflag] = fsolve(@equations, x0, opts);
    
    if exitflag > 0
        V0     = xsol(1);
        alpha0 = xsol(2);
        gamma0 = xsol(3);
        ok = true;
    else
        eigs    = [];
        trim    = struct('V',NaN,'alpha',NaN,'gamma',NaN,'theta',NaN,'q',NaN);
        ok = false;
        return;
    end
    
    %--- build trim struct & state
    theta0 = alpha0 + gamma0;
    q0 = 0;
    trim = struct('V', V0, 'alpha', alpha0, 'gamma', gamma0, 'theta', theta0, 'q', q0);
    
    x_eq = [V0; gamma0; theta0; q0];
    
    %--- numeric Jacobian of the full 4×4 f
    A    = numericJacobian(@stateDeriv, x_eq);
    
    %--- eigenvalues
    eigs = eig(A);
    
    
    % nested: your three trim equations
    function F = equations(x)
        V     = x(1);
        alpha = x(2);
        gamma = x(3);
        M     = calculateMoment(p, V, alpha);
        [L,D] = calculateAeroForces(alpha, V, p);
        Vdot    = -D/p.m_total - p.g*sin(gamma);
        gammadot=  L/(p.m_total*V) - p.g*cos(gamma)/V;
        F = [ M; Vdot; gammadot ];
    end

    % nested: full state derivative
    function dx = stateDeriv(x)
        % unpacking variables
        V = x(1);
        gamma = x(2);
        theta = x(3);
        q = x(4);
        alpha = theta - gamma;
        
        % Calculating derivatives
        [L,D] = calculateAeroForces(alpha, V, p);
        M = calculateMoment(p, V, alpha);

        % Compiling results
        dx = [-D/p.m_total - p.g*sin(gamma);             % Vdot
            L/(p.m_total*V) - p.g*cos(gamma)/V;         % gammadot
            q;                                          % thetadot
            M/p.I_total                                % qdot
            ];
    end

% nested: central‐difference Jacobian
    function J = numericJacobian(fun, x0)
        eps0 = 1e-6;
        n = numel(x0);
        J = zeros(n);
        for i = 1:n
            h = eps0 * max(abs(x0(i)),1);
            xp = x0;  xm = x0;
            xp(i) = xp(i) + h;
            xm(i) = xm(i) - h;
            fp = fun(xp);
            fm = fun(xm);
            J(:,i) = (fp - fm) / (2*h);
        end
    end
end
