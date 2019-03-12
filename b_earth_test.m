%B_EARTH_TEST Unit tests for b_earth.m
% 
% Other m-files required: b_earth.m
% Subfunctions: none
% MAT-files required: none
%
% Author: Liam Shepherd
% Mar 2019; Last revision: 12-Mar-2019

%% Clear down
clear
clc

%% Test dipole direction at central axis

% Position to test at
r_axis = [0; 0; 1]; % Set 1m above axis as value is undef at origin

% Get B value at test point
B = b_earth(r_axis);

% Output B vector
disp('B at z-axis:');
disp(B);

% Throw error if B_x or B_y aren't zero
if or(not(isequal(B(1), 0)), not(isequal(B(2), 0)))
    warning('B_x & B_y on z-axis not zero');
end
% Check B_z is anti-parallel to z-axis
if not(B(3) < 0)
    warning('B_z on z_axis not in -ve z direction');
end
    
%% Test dipole direction in X/Y plane far from X-axis

% Test position vectors 
r_tests = horzcat([1e10; 0; 0], [-1e10; 1e10; 0],  [0; -1e10; 0]);

% Iterate over r_test columns
for n=1:size(r_tests, 2)
    r = r_tests(:, n);
    
    % Output r position
    disp('Checking B value at:');
    disp(r);
    
    % Get B value at test point
    B = b_earth(r);
    
    % Output B vector
    disp('B at point:');
    disp(B);
    
    % Throw error if B_x or B_y aren't zero
    if or(not(isequal(B(1), 0)), not(isequal(B(2), 0)))
        warning('B_x & B_y on X/Y plane far from z-axis not zero');
    end
    % Check B_z is parallel to z-axis
    if not(B(3) > 0)
        warning('B_z on X/Y plane far from z-axis not in +ve z direction');
    end
end