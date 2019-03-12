%F_LORENTZ_TEST Unit tests for f_lorentz.m
% 
% Other m-files required: f_lorentz.m
% Subfunctions: none
% MAT-files required: none
%
% Author: Liam Shepherd
% Mar 2019; Last revision: 12-Mar-2019

%% Initialise

% Clear down
clear
clc

% Test charge
q = 10;

%% Test stationary charge

% Test velocity
v = [0; 0; 0];

% Moving: no; Electric: yes; Magnetic: no
if not(isequal(f_lorentz([5; 5; 5], [0; 0; 0], q, v), [50; 50; 50]))
    warning('Unexpected force for stationary charge with electric only');
end

% Moving: no; Electric: no; Magnetic: yes
if not(isequal(f_lorentz([0; 0; 0], [5; 5; 5], q, v), [0; 0; 0]))
    warning('Unexpected force for stationary charge with magnetic only');
end

% Moving: no; Electric: yes; Magnetic: yes
if not(isequal(f_lorentz([5; 5; 5], [5; 5; 5], q, v), [50; 50; 50]))
    warning('Unexpected force for stationary charge with both fields');
end

%% Test moving charge

% Test velocity
v = [10; 10; 0];

% Moving: yes; Electric: yes; Magnetic: no
if not(isequal(f_lorentz([5; 5; 5], [0; 0; 0], q, v), [50; 50; 50]))
    warning('Unexpected force for moving charge with electric only');
end

% Moving: yes; Electric: no; Magnetic: yes
if not(isequal(f_lorentz([0; 0; 0], [0; 5; 5], q, v), [500; -500; 500]))
    f_lorentz([0; 0; 0], [0; 5; 5], q, v)
    warning('Unexpected force for moving charge with magnetic only');
end

% Moving: yes; Electric: yes; Magnetic: yes
if not(isequal(f_lorentz([5; 5; 5], [0; 5; 5], q, v), [550; -450; 550]))
    warning('Unexpected force for moving charge with both fields');
end
