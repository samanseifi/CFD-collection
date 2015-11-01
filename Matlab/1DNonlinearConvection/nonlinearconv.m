function [x, u] = nonlinearconv(scheme, l, nx, nt, dt)
% Solving the nonlinear convection equation using different numerical
% schemes: Lax-Wendroff, Lax-Friedrichs and MacCormak
%
% nonlinear convection: du/dt = u * du/dx
%
% sample: nonlinearconv('LaxWendroff', 2*pi, 201, 200, 0.0016)
%
% choices for 'scheme': 'LaxWendroffTwoStep', 'LaxWendroff', 'LaxFriedrichs', 'MacCormak'

dx = l/(nx - 1);
x = 0:dx:l;
%dt = 0.0016;

% Initial Condition
vis = 0.1;
u_init = DoubleGaussian(x, nx, vis);
%u_init = Step(x,nx);
u = u_init;


for i = 1:nx
    ip(i) = i+1;
    im(i) = i-1;
end
ip(nx) = 1;
im(1) = nx;


if strcmp(scheme, 'LaxWendroffTwoStep')
    for t = 1:nt  
        for i = 1:nx
            % Lax-Wendroff two step
            u1 = 0.5*(u(i) + u(ip(i))) - 0.25*(dt/dx)*(u(ip(i))^2 - u(i)^2);
            u2 = 0.5*(u(im(i)) + u(i)) - 0.25*(dt/dx)*(u(i)^2 - u(im(i))^2);                              
            u_new(i) = u(i) - 0.5*(dt/dx)*(u1^2 - u2^2); 
        end
        u = u_new;
    end
end

if strcmp(scheme, 'LaxWendroff')
    for t = 1:nt  
        for i = 1:nx
            % Lax-Wendroff
            u1 = (u(i) + (u(i) + u(ip(i)))*0.5)*0.5;
            u2 = (u(i) + (u(im(i)) + u(i))*0.5)*0.5;
            u_new(i) = u(i) - 0.25*(dt/dx)*(u(ip(i))^2 - u(im(i))^2) + 0.5*(dt^2/dx^2)*(u1*(u(ip(i))^2 - u(i)^2) - u2*(u(i)^2 - u(im(i))^2));
        end
        u = u_new;
    end
end

if strcmp(scheme, 'LaxFriedrichs')
    for t = 1:nt  
        for i = 1:nx
            % Lax–Friedrichs
            u_new(i) = 0.5*(u(ip(i)) + u(im(i))) - 0.25*(dt/dx)*(u(ip(i))^2 - u(im(i))^2);
        end
        u = u_new;
    end
end

if strcmp(scheme, 'MacCormak')
    u_star = u;
    for t = 1:nt  
        for i = 1:nx
            % Mac Cormak
            u_star(i) = u(i) - (dt/dx)*0.5*(u(ip(i))^2 - u(i)^2); % Predictor
            u_new(i) = 0.5*(u(i) + u_star(i)) - 0.5*(dt/dx)*0.5*(u_star(i)^2 - u_star(im(i))^2); % Corrector
        end
        u = u_new;
	end
end

end

   
