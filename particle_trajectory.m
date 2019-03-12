function [dpath_dt] = particle_trajectory(q, m, path)
    %PARTICLE_TRAJECTORY - Calculate change of particle trajectory
    %
    % Syntax: dpath_dt = particle_trajectory(q, m, path)
    %
    % Inputs:
    %   q    - Charge of particle in C
    %   m    - Mass of particle in kg
    %   path - Cartesian column vector of particle position and velocity
    %          [r_x; r_y; r_z; v_x; v_y; v_z]
    %
    % Outputs:
    %   d_path_dt - First derivate of particle trajectory
    %               [v_x; v_y; v_z; a_x; a_y; a_z]
    %
    % Example: 
    %   dpath_dt = particle_trajectory(10, 20, [10; 10; 0; 5; 0; 5])
    %
    % Other m-files required: b_earth.m f_lorentz.m
    % Subfunctions: none
    % MAT-files required: none
    %
    % Author: Liam Shepherd
    % Mar 2019; Last revision: 12-Mar-2019

    % Check input(s) are of correct type or bail otherwise
    if (not(isscalar(q)))
        error('q should be a scalar');
    end
    if (not(isscalar(m)))
        error('m should be a scalar');
    end
    if (not(isequal(size(path), [6, 1])))
        error('path should be a 6x1 column vector');
    end

    % Position column vector
    r = path(1:3);
    % Velocity column vector
    v = path(4:6);
    
    % Approximate there being no E field
    E = [0; 0; 0];
    % Calculate B field at this position
    B = b_earth(r);
    
    % Work out force on particle
    F = f_lorentz(E, B, q, v);

    % Return column vector of velocity and acceleration
    dpath_dt = [v; F/m];
end

