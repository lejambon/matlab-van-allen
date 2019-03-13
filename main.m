%MAIN Plot particle trajectories in Earth's magnetic field
% 
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Liam Shepherd
% Mar 2019; Last revision: 13-Mar-2019

%% Initialise

% Clear down
clear
clc
clf

% Constants
R_e = 6.371e6;   % Radius of Earth in m
q_e = 1.602e-19; % Elementary charge in C
m_p = 1.673e-27; % Mass of proton in kg
m_e = 9.109e-31; % Mass of electron in kg
c   = 2.998e8;   % Speed of light in m s^-1

% Helper functions
% Convert eV to J
eVtoJ = @(ev) ev*1.60218e-19;
% Get velocity from relativistic kinetic energy
energytovelocity = @(e, m) c*sqrt(1 - 1/(e/(m*c^2) + 1)^2); %

%% Proton (Inner belt)

% 50 MeV proton
proton_eV = 50e6;
% Get the velocity from the energy
v_p = energytovelocity(eVtoJ(proton_eV), m_p);
% Initial conditions for proton: [r_x r_y r_z v_x v_y v_z]
%  Place proton 1 Earth radius away from equator
%  Put velocity at 45 degree angle with magnetic field in X/Z
path_p0 = [2*R_e; 0; 0; sind(45)*v_p; 0; cosd(45)*v_p];
% Setup differential equation to solve for proton
trajectory_p = @(t, s) particle_trajectory(+q_e, m_p, s);
% Solve for proton path over time 0s to 52s
[t_p, path_p] = ode45(trajectory_p, [0 52], path_p0);

%% Electron (Outer belt)

% 10 MeV electron
electron_eV = 10e6;
% Get the velocity from the energy
v_e = energytovelocity(eVtoJ(electron_eV), m_e);
% Initial conditions for electron: [r_x r_y r_z v_x v_y v_z]
%  Place electron 7 Earth radii away from equator
%  Put velocity at 30 degree angle with magnetic field in X/Z
path_e0 = [8*R_e; 0; 0; sind(30)*v_e; 0; cosd(30)*v_e];
% Setup differential equation to solve for electron
trajectory_e = @(t, s) particle_trajectory(-q_e, m_e, s);

% Below are three different unsuccessful methods for getting electron path

% Solve for electron path using ODE function for stiff eqn from 0 to 1s
%[t_e, path_e] = ode23s(trajectory_e, [0 1], path_e0);

% Euler method, large step. Produces jerky very broken motion
% Note that ode1.m came from https://uk.mathworks.com/matlabcentral/answers/98293-is-there-a-fixed-step-ordinary-differential-equation-ode-solver-in-matlab-8-0-r2012b
%path_e = ode1(trajectory_e, linspace(0, 1, 10000), electron_p0);

% Euler method, short step. Produces good motion for a short bit of the
%  trajectory. Gets more chaotic later if time range expanded. Plot heavy
%  due to huge number of points
% Note that ode1.m came from https://uk.mathworks.com/matlabcentral/answers/98293-is-there-a-fixed-step-ordinary-differential-equation-ode-solver-in-matlab-8-0-r2012b
%path_e = ode1(trajectory_e, linspace(0, 10, 100000000), electron_p0);

%% Plot

% Create new figure
figure(1);

% Plot the proton path in red using units of Earth's radius
plot3(path_p(:,1)./R_e , path_p(:,2)./R_e , path_p(:,3)./R_e ,...
    'LineWidth',1.0,'Color','Red');

% If we have an electron path to plot ...
if exist('path_e', 'var')
    scale = 9; % ... then set plot size as 9 Earth radii
else
	scale = 2.5; % Otherwise set plot size as 2.5 Earth radii
end
% Set the axis limits in all directions based on above
axis([-scale scale -scale scale -scale scale]);
% Axis labels
xlabel('x / R_e');
ylabel('y / R_e');
zlabel('z / R_e');
% Set title
title("Particle trajectory in Earth's magnetic field");

% If we have an electron path to plot ...
if exist('path_e', 'var')
    hold on; % Same plot
    % Plot the electron path in blue using units of Earth's radius
    plot3(path_e(:,1)./R_e, path_e(:,2)./R_e, path_e(:,3)./R_e,...
        'LineWidth',1.0,'Color','Blue');
end

% Plot the Earth
hold on; % Same axis
colormap winter; % Blue/Green colours
[x,y,z] = sphere; % Unit sphere
surf(x, y, z); % Surface plot unit sphere