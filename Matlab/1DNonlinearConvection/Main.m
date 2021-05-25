% This script runs the "nonlinearconv.m" script for solving the nonlinear
% convection equation in 1D with different finite difference schems.

clear
close all

% Defining the domain of solution
l = 2*pi;

% Number of spatial grid points
nx = 101;

% Number of time steps
nt = 200;
dt = 0.0016;

% Discritize the domain
x = linspace(0, l, nx);

% Create the initial condition
vis = 0.1;
u_initial = Initial(x, nx, vis);

% Solving Non-linear convection equation and store the final time result!
[x1, u1] = nonlinearconv('LaxWendroffTwoStep', l, nx, nt, dt, u_initial);
[x2, u2] = nonlinearconv('LaxWendroff', l, nx, nt, dt, u_initial);
[x3, u3] = nonlinearconv('LaxFriedrichs', l, nx, nt, dt, u_initial);
[x4, u4] = nonlinearconv('MacCormak', l, nx, nt, dt, u_initial);

% Initial profile
x_init = x1;
u_init = Initial(x_init, nx, 0.1); % vis = 0.1;

% Setting up the plotting envirnoment
xmax = max(max([x1; x2; x3; x4]));
xmin = min(min([x1; x2; x3; x4]));
ymax = max(max([u1; u2; u3; u4; u_init]));
ymin = min(min([u1; u2; u3; u4; u_init]));

% Plotting
plot(x_init, u_init, '--k', x1, u1, '-o', x2, u2, '-s', x3, u3, '-^', x4, u4, '-*')
title('1D Non-linear convection');
xlabel('Spatial coordinate: $x$', 'Interpreter', 'latex');  % Using latex intepreteer
ylabel('Transport property: $u$', 'Interpreter', 'latex');  % Using latex intepreteer


axis([xmin xmax ymin-0.5 ymax+0.5])

legend('Initial', 'LaxWendroffTwoStep', 'LaxWendroff', 'LaxFriedrichs', 'MacCormak', 'Location', 'northwest')

