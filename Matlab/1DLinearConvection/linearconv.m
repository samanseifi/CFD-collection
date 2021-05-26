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
N = 200;        % Number of grid points

t_final = 10;   % Final time of simulation!
nt = 100;      % Number of time steps

% Wave speed (The speed of wave in the given domain)
c = 0.1;

% Discritized domain spatial and temporal
x = linspace(0, L, N);
dx = L/(N - 1);     % Alternatively dx = x(2) - x(1)
dt = t_final/nt;

% Frist check the CFL codition (This is the same for all explicit schemes)
if c <= dx/dt
    fprintf("NOTE: The explicit method is STABLE!\n")
else
    fprintf("WARNING: The explicit method is UNSTABLE!\n")
end

% Generate the shape of the intial condition
u_0 = zeros(N, 1);
u_0(4: floor(N/5)) = 1.0;

% Setup the visualization
figure

% Preallocating the solution u_new
u1_new = zeros(N, 1);
u2_new = zeros(N, 1);

u1 = sparse(u_0);  % Apply the IC to the system
u2 = sparse(u_0);  % Apply the IC to the system

t = 0;    % Initialize time

% Solving the problem with matrices: [u]_k+1 = {A}[u]_k
% For explicit schemes matrix {A} is constant
% 1) Upwind Scheme
% **<@
dm = c * (dt/dx);           % ***> Off diagonal -1 values
d = 1 - c * (dt/dx);        % ***> Diagonal values

% Constructing the matrix!
A_upwind = diag(d * ones(1, N)) + diag(dm * ones(1, N-1), -1);
A_upwind(1, N) = dm;
A_upwind = sparse(A_upwind);
% **@>

% 2) Lax–Friedrichs method
% **<@
dm = 0.5*(1 + c*(dt/dx));   % ***> Off diagonal -1 values
dp = 0.5*(1 - c*(dt/dx));   % ***> Off diagonal +1 values

% Constructing the matrix!
A_LaxFr = diag(dm * ones(1, N-1), -1) + diag(dp * ones(1, N-1), 1);
A_LaxFr(1, N) = dm;
A_LaxFr(N, 1) = dp;
A_LaxFr = sparse(A_LaxFr);
% **@>


% Marching in time!
while (t < t_final)
    
    % Upwind scheme!
    u1_new = A_upwind * u1;
    
    % Lax–Friedrichs scheme1
    u2_new = A_LaxFr * u2;
    
    % Store new solutions
    u1 = u1_new;
    u2 = u2_new;
    
    % Step in time!
	t = t + dt;
    
    % Plotting live!
    plot(x, u1, '-bo', x, u2, '--ks');
    
    % Setup the limits and labels
    xlim([0 L]) 
    ylim([0 max(u_0)])
    xlabel('x')
    ylabel('u(x)')
	
    drawnow;
end