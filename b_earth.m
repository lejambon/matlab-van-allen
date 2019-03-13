function [B] = b_earth(r)
    %B_EARTH Calculate Earth's magnetic field at a position
    % Approximates Earth's magnetic field as a dipole
    %
    % Syntax:  B = b_earth(r)
    %
    % Inputs:
    %   r - Cartesian column vector of position to calculate at in m
    %
    % Outputs:
    %   B - Cartesian column vector of magnetic field in T
    %
    % Example: 
    %   B = b_earth([100; 100; 0]);
    %
    % Other m-files required: b_dipole.m
    % Subfunctions: none
    % MAT-files required: none
    %
    % Author: Liam Shepherd
    % Mar 2019; Last revision: 12-Mar-2019
    
    % Check input(s) are of correct type or bail otherwise
    if (not(isequal(size(r), [3, 1])))
        error('r should be a 3x1 column vector');
    end

    % Magnetic dipole of Earth in A m^2
    % Ref: http://wdc.kugi.kyoto-u.ac.jp/poles/polesexp.html
    m_earth = [0; 0; -7.71e22];
    
    % Magnetic field of Earth at point r in T
    B = b_dipole(m_earth, r);
end