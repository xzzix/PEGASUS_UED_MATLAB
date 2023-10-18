function [v1, v2, v3] = get_lattice_vectors(a, b, c, alpha, beta, gamma)
    % Input angles should be in degrees
    % Convert all angles to radians for computation
    alpha = deg2rad(alpha);
    beta = deg2rad(beta);
    gamma = deg2rad(gamma);
    % Calculate the components of each lattice vector
    v1 = [a, 0, 0];
    v2 = [b*cos(gamma), b*sin(gamma), 0];
    v3_x = c*cos(beta);
    v3_y = c*(cos(alpha) - cos(beta)*cos(gamma))/sin(gamma);
    v3_z = sqrt(c^2 - v3_x^2 - v3_y^2);
    v3 = [v3_x, v3_y, v3_z];
end