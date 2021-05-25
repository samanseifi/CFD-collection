function [x, u] = nonlinearconv(scheme, L, nx, nt, dt, initial)
% Solving the nonlinear convection equation in 1D using different numerical
% schemes: 
%
%   1) Lax-Wendroff, 
%   2) Lax-Friedrichs,
%   3) MacCormak
%
% choices for 'scheme': 'LaxWendroffTwoStep', 'LaxWendroff', 'LaxFriedrichs', 'MacCormak'
%
% nonlinear convection:     0 <= x <= L
%
%   du/dt = u * du/dx
%
% with periodic BC: u(0) = u(L) and initial condition of 
%
% Example run:
%
%   nonlinearconv('LaxWendroff', 2*pi, 201, 200, 0.0016)
%

% Discritized domain spatial
x = linspace(0, L, nx);
dx = L/(nx - 1);

% Initial Condition
u_0 = initial;

% Setup the 1D stencil accounting for periodic BC
ip = zeros(1, nx);
im = zeros(1, nx);
for i = 1:nx
    ip(i) = i+1;
    im(i) = i-1;
end
ip(nx) = 1;
im(1) = nx;


u = u_0;

u_new = zeros(size(u_0));

if strcmp(scheme, 'LaxWendroffTwoStep')
    for t = 1:nt  
        for i = 1:nx
            % Lax-Wendroff
            % First Step
            u1 = 0.5*(u(i) + u(ip(i))) - 0.25*(dt/dx)*(u(ip(i))^2 - u(i)^2);
            
            % Second Step
            u2 = 0.5*(u(im(i)) + u(i)) - 0.25*(dt/dx)*(u(i)^2 - u(im(i))^2);                              
            
            % Update the solution!
            u_new(i) = u(i) - 0.5*(dt/dx)*(u1^2 - u2^2); 
        end
        u = u_new;
    end

elseif strcmp(scheme, 'LaxWendroff')
    for t = 1:nt  
        for i = 1:nx
            % Lax-Wendroff
            u1 = (u(i) + (u(i) + u(ip(i)))*0.5)*0.5;
            u2 = (u(i) + (u(im(i)) + u(i))*0.5)*0.5;
            u_new(i) = u(i) - 0.25*(dt/dx)*(u(ip(i))^2 - u(im(i))^2) + 0.5*(dt^2/dx^2)*(u1*(u(ip(i))^2 - u(i)^2) - u2*(u(i)^2 - u(im(i))^2));
        end
        u = u_new;
    end

elseif strcmp(scheme, 'LaxFriedrichs')
    for t = 1:nt  
        for i = 1:nx
            % Lax–Friedrichs
            u_new(i) = 0.5*(u(ip(i)) + u(im(i))) - 0.25*(dt/dx)*(u(ip(i))^2 - u(im(i))^2);
        end
        u = u_new;
    end

elseif strcmp(scheme, 'MacCormak')
    u_star = u;
    for t = 1:nt      
        for i = 1:nx
            % Mac Cormak
            % Predictor
            u_star(i) = u(i) - (dt/dx)*0.5*(u(ip(i))^2 - u(i)^2); 
            
            % Corrector
            u_new(i) = 0.5*(u(i) + u_star(i)) - 0.5*(dt/dx)*0.5*(u_star(i)^2 - u_star(im(i))^2); 
        end
        u = u_new;
	end
end   
