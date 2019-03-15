%MAIN Plot particle trajectories in Earth's magnetic field
% 
% Other m-files required: particle_trajectory.m
% Subfunctions: none
% MAT-files required: none
%
% Author: Liam Shepherd
% Mar 2019; Last revision: 14-Mar-2019

%% Initialise

% Clear down
clear
clc
close all

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
energytovelocity = @(e, m) c*sqrt(1 - 1/(e/(m*c^2) + 1)^2);

%% Init plot

% Create new figure
figure(1);

% Same figure for everything
hold all;

% Set viewport as 45 degrees azimuth and 15 degress elevation
view(45, 15);

% Plot axis size in Earth radii
scale = 6;
% Set the axis limits in all directions based on above
axis([-scale scale -scale scale -scale scale]);

% Add ticks at every 2nd Earth radii
xticks(-floor(scale):2:floor(scale));
yticks(-floor(scale):2:floor(scale));
zticks(-floor(scale):2:floor(scale));

% Draw grid
grid on;

% Axis labels
xlabel('x / R_e');
ylabel('y / R_e');
zlabel('z / R_e');

% Set title
title("Particle trajectory in Earth's magnetic field");

% Plot the Earth
colormap winter; % Blue/Green colours
[x,y,z] = sphere; % Unit sphere
surf(x, y, z); % Surface plot unit sphere

% Flush the plot buffer
drawnow;

%% Proton (Inner belt)

% 50 MeV proton
proton_eV = 50e6;
% Get the velocity from the energy
v_p = energytovelocity(eVtoJ(proton_eV), m_p);
% Initial conditions for proton: [r_x r_y r_z v_x v_y v_z]
%  Place proton 1 Earth radius away from equator
%  Put velocity at 45 degree angle with magnetic field in X/Z
path_p0 = [2*R_e; 0; 0; sind(45)*v_p; 0; cosd(45)*v_p];
% Timespan to solve for (roughly one Earth revolution)
t_max = 32.3;

% Setup differential equation to solve for proton
trajectory_p = @(t, s) particle_trajectory(+q_e, m_p, s);

% Decrease relative tolerance when ODE solving from default of 1e-3
%  to stop proton bounce having different height after one revolution
opts = odeset('RelTol',1e-4);

% Solve the proton path in steps to show motion
% Step size for each loop
dt = 0.5;
% For each [t, t + dt] time step
for t=0:dt:t_max-dt
    % Output where we've got up to
    fprintf('Solving proton path in inner belt for %.1f to %.1f s\n',...
        t, t+dt);
    % Solve for this step
    [~, path_p] = ode45(trajectory_p, [t t+dt], path_p0, opts);

    % Set the next loop initial conditions to the last point in this solve
    path_p0 = path_p(end, :);

    % Plot the proton path in red using units of Earth's radius
    plot3(path_p(:,1)./R_e, path_p(:,2)./R_e, path_p(:,3)./R_e,...
        'Color', 'Red');

    % Flush the plot buffer
    drawnow;
end

%% Electron (Outer belt)

% 10 MeV electron
electron_eV = 10e6;
% Get the velocity from the energy
v_e = energytovelocity(eVtoJ(electron_eV), m_e);
% Initial conditions for electron: [r_x r_y r_z v_x v_y v_z]
%  Place electron 4.5 Earth radii away from equator
%  Put velocity at 45 degree angle with magnetic field in X/Z
path_e0 = [5.5*R_e; 0; 0; sind(45)*v_e; 0; cosd(45)*v_e];
% Timespan to solve for
t_max = 1200;

% Setup differential equation to solve for electron
trajectory_e = @(t, s) particle_trajectory(-q_e, m_e, s);

% Decrease ODE solving relative tolerance to reduce bounce height drifting
%  It's still not perfect but decreasing further increases solve time
opts = odeset('RelTol', 5e-7);

% Solving the electron path in one go takes a huge amount of memory
%  and gives no feedback whilst solving. Instead break into small chunks
% Step size for each loop
dt = 0.1;
% For each [t, t + dt] time step
for t=0:dt:t_max-dt
    % Output where we've got up to
    fprintf('Solving electron path in outer belt for %.1f to %.1f s\n',...
        t, t+dt);
    % Solve for this step
    [~, path_e] = ode45(trajectory_e, [t t+dt], path_e0, opts);

    % Set the next loop initial conditions to the last point in this solve
    path_e0 = path_e(end, :);

    % Plotting all points becomes extremely slow. Instead strip only
    %  every 100th point. This breaks the smooth spiral nature of the
    %  electron path but makes the plot dragging more usable
    path_e=path_e(1:100:end, :);

    % Plot the electron path in blue using units of Earth's radius
    plot3(path_e(:,1)./R_e, path_e(:,2)./R_e, path_e(:,3)./R_e,...
        'Color', 'Blue');

    % Flush the plot buffer
    drawnow;
end