% Solving the Linear Convection problem in 1D:
%
%       du/dt + c du/dx = 0
%
% for 0 <= x <= L. And with periodic BC:
%       u(0) = u(L)
%
% with IC @ t = 0 with u(x, 0) = u_0(x)
% where u_0(x): a simple step function.

clear 
close all

% Input parameter
L = 1.0;        % Length of the domain
N = 100;        % Number of grid points

t_final = 10;   % Final time of simulation!
nt = 100;      % Number of time steps

% Wave speed (The speed of wave in the given domain)
c = 0.1;

% Discritized domain spatial and temporal
x = linspace(0, L, N);
dx = L/(N - 1);     % Alternatively dx = x(2) - x(1)
dt = t_final/nt;

% Generate the shape of the intial condition
u_0 = zeros(1, N);
u_0(4: floor(N/5)) = 1.0;

% Setup the 1D stencil accounting for periodic BC
ip = zeros(1, N);
im = zeros(1, N);
for i = 1:N
    ip(i) = i+1;
    im(i) = i-1;
end
ip(N) = 1;      % The final grid point to the first
im(1) = N;      % The first grid point to the final

% Setup the visualization
figure


% Preallocating the solution u_new
u1_new = zeros(1, N);
u2_new = zeros(1, N);

u1 = u_0;  % Apply the IC to the system
u2 = u_0;  % Apply the IC to the system

t = 0;    % Initialize time

% Marching in time!
while (t < t_final)
    
    % Update solution for the new time step
    for i = 1:N   
        
        % Simple explicit and forward in time mrthod
        u1_new(i) = u1(i) - c*(dt/dx)*(u1(i) - u1(im(i)));
        
        % Lax–Friedrichs method 
        u2_new(i) = 0.5*(u2(ip(i)) + u2(im(i))) - c*0.5*(dt/dx)*(u2(ip(i)) - u2(im(i)));
        
    end
    
    % Store new solutions
    u1 = u1_new;
    u2 = u2_new;
    
    % Step in time!
	t = t + dt;
    
    % Plotting live!
    plot(x, u1, '-o', x, u2, '-ro');
    
    % Setup the limits and labels
    xlim([0 L]) 
    ylim([0 max(u_0)])
    xlabel('x')
    ylabel('u(x)')
	
    drawnow;
end