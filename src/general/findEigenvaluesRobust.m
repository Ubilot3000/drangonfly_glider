function [eigs, A] = findEigenvaluesRobust(p, v0)
  % 1) Trim for alpha & gamma at V0
  [alpha0, gamma0] = findTrimCondition(v0, p);

  % pack the equilibrium state (theta0 = alpha0 + gamma0, q0 = 0)
  x0 = [ v0; gamma0; alpha0 + gamma0; 0];

  % 2) Define the state‐derivative function
  f = @(x) stateDeriv(x, p);

  % 3) Numerical Jacobian (central difference)
  A = numericJacobian(f, x0);

  % 4) Eigen‐analysis
  eigs = eig(A);
end

%----------------------------------------------------------------------
function dx = stateDeriv(x, p)
  % x = [V; gamma; theta; q]
  V     = x(1);
  gamma = x(2);
  theta = x(3);
  q     = x(4);

  % angle of attack
  alpha = theta - gamma;

  % aerodynamic forces & moment
  [L, D] = calculateAeroForces(alpha, V, p);
  M      = calculateMoment(p, V, alpha);

  % eqns of motion
  Vdot      = -D/p.m_total - p.g * sin(gamma);
  gammadot  =  L/(p.m_total*V) - p.g * cos(gamma)/V;
  thetadot  =  q;
  qdot      =  M/p.I_total;

  dx = [Vdot; gammadot; thetadot; qdot];
end

%----------------------------------------------------------------------
function J = numericJacobian(fun, x0)
  % central‐difference Jacobian
  eps0 = 1e-6;
  n    = numel(x0);
  fx0  = fun(x0);
  J    = zeros(n);
  for i = 1:n
    h    = eps0 * max(abs(x0(i)),1);
    xp   = x0; xm = x0;
    xp(i)= x0(i) + h;
    xm(i)= x0(i) - h;
    fp   = fun(xp);
    fm   = fun(xm);
    J(:,i) = (fp - fm) / (2*h);
  end
end
