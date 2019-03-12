function [F] = f_lorentz(E, B, q, v)
    %F_LORENTZ - Calculate Lorentz force on a particle
    %
    % Syntax: F = f_lorentz(E, B, q, v)
    %
    % Inputs:
    %   E - Cartesian column vector of electric field in N C^-1
    %   B - Cartesian column vector of magnetic field in T
    %   q - Charge of particle in C
    %   v - Cartesian column vector of particle velocity in m s^-1
    %
    % Outputs:
    %   F - Cartesian column vector of force on particle in N
    %
    % Example: 
    %    F = f_lorentz([100; 0; 0], [0; 100; 0], 10, [100; 5; 100])
    %
    % Other m-files required: none
    % Subfunctions: none
    % MAT-files required: none
    %
    % Author: Liam Shepherd
    % Mar 2019; Last revision: 12-Mar-2019

    % Check input(s) are of correct type or bail otherwise
    if (not(isequal(size(E), [3, 1])))
        error('E should be a 3x1 column vector');
    end
    if (not(isequal(size(B), [3, 1])))
        error('B should be a 3x1 column vector');
    end
    if (not(isscalar(q)))
        error('q should be a scalar');
    end
    if (not(isequal(size(v), [3, 1])))
        error('v should be a 3x1 column vector');
    end

    % Calculate Lorentz force on particle
    F = q*(E + cross(v, B));
end

