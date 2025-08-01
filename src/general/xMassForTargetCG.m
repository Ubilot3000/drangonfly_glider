function x_mass = xMassForTargetCG(target_cg, c)
    % This function calculates the required position of a movable mass to
    % achieve a desired overall Center of Gravity for the aircraft.
    % It solves the Center of Gravity equation for x_mass.

    % Sum of moments
    fixed_moments = c.m_nose * c.x_nose + c.m_f_wing * c.x_f_wing + c.m_r_wing * c.x_r_wing + c.m_rudder * c.x_rudder + c.m_rod * c.x_rod;

    % Sum of masses
    m_total = c.m_nose + c.m_f_wing + c.m_r_wing + ...
              c.m_rudder + c.m_rod + c.m_mass;

    % Solve for x_mass
    x_mass = (target_cg * m_total - fixed_moments) / c.m_mass;
end