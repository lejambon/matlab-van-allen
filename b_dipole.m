function [B] = b_dipole(m, r)
    %B_DIPOLE Calculate magnetic field at position from dipole at origin
    % Note that magnetic field at origin is undefined
    %
    % Syntax:  B = b_dipole(m, r)
    %
    % Inputs:
    %   m - Cartesian column vector of magnetic moment in A m^2
    %   r - Cartesian column vector of position to calculate at in m
    %
    % Outputs:
    %   B - Cartesian column vector of magnetic field in T
    %
    % Example: 
    %    B = b_dipole([0; 0; 100], [100; 100; 0]);
    %
    % Other m-files required: none
    % Subfunctions: none
    % MAT-files required: none
    %
    % Author: Liam Shepherd
    % Mar 2019; Last revision: 12-Mar-2019

    % Check input(s) are of correct type or bail otherwise
    if (not(isequal(size(m), [3, 1])))
        error('m should be a 3x1 column vector');
    end
    if (not(isequal(size(r), [3, 1])))
        error('r should be a 3x1 column vector');
    end
    
    % Error if position given is the origin
    if isequal(r, [0; 0; 0])
        error('Magnetic field undefined at origin');
    end

    % Vacuum permeability in T m A^-1
    u_0 = 1.257e-6;
    
    % Calculate dipole magnetic field at position using magnetic moment 
    B = (u_0)/(4*pi)*(3*r.*(dot(m, r))/(norm(r)^5) - m/(norm(r)^3));
end