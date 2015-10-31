% Solving the nonlinear convection equation using different numerical
% schemes: Lax-Wendroff, Lax-Friedrichs and MacCormak
%
% nonlinear convection: du/dt = u * du/dx
%

clear all
close all

Input;

% Initial Condition
vis = 0.1;
u = DoubleGaussian(x, nx, vis);

% Choice of scheme
%scheme = 'LaxWendroffTwoStep';
%scheme = 'LaxWendroff'; 
%scheme = 'LaxFriedrichs'; 
scheme = 'MacCormak';

for i = 1:nx
    ip(i) = i+1;
    im(i) = i-1;
end
ip(nx) = 1;
im(1) = nx;


for t = 1:nt  
    for i = 1:nx
        if strcmp(scheme, 'LaxWendroffTwoStep')
            % Lax-Wendroff two step
            u1 = 0.5*(u(i) + u(ip(i))) - 0.25*(dt/dx)*(u(ip(i))^2 - u(i)^2);
            u2 = 0.5*(u(im(i)) + u(i)) - 0.25*(dt/dx)*(u(i)^2 - u(im(i))^2);                              
            u_new(i) = u(i) - 0.5*(dt/dx)*(u1^2 - u2^2); 
        elseif strcmp(scheme, 'LaxWendroff')
            % Lax-Wendroff
            u1 = (u(i) + (u(i) + u(ip(i)))*0.5)*0.5;
            u2 = (u(i) + (u(im(i)) + u(i))*0.5)*0.5;
            u_new(i) = u(i) - 0.25*(dt/dx)*(u(ip(i))^2 - u(im(i))^2) + 0.5*(dt^2/dx^2)*(u1*(u(ip(i))^2 - u(i)^2) - u2*(u(i)^2 - u(im(i))^2));
        elseif strcmp(scheme, 'LaxFriedrichs')
            % Lax–Friedrichs
            u_new(i) = 0.5*(u(ip(i)) + u(im(i))) - 0.25*(dt/dx)*(u(ip(i))^2 - u(im(i))^2);
        elseif strcmp(scheme, 'MacCormak')
            % Mac Cormak
            us = u(i) - (dt/dx)*0.5*(u(ip(i))^2 - u(i)^2);
            u_new(i) = 0.5*(u(i) + us) - 0.5*(dt/dx)*0.5*(u(i)^2 - u(im(i))^2);
        end
    end

    u = u_new;
    
end
plot(x, u)
