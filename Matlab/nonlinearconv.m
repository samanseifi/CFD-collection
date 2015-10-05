%clear all
%close all
hold on

nx = 201;
nt = 200;
l = 2*pi;
dx = l/(nx - 1);
x = 0:dx:l;
c = 1.0;
dt = 0.0016;
u_init = zeros(1,nx);
u = zeros(1,nx);
u_new = zeros(1,nx);

vis = 0.1;

ip = zeros(nx);
im = zeros(nx);

for i = 1:nx
	ip(i) = i + 1;
	im(i) = i - 1;
	phi(i) = exp(-0.25*(x(i)*x(i))/vis) + exp(-0.25*((2*pi - x(i))^2)/vis);
	dphi(i) = (-0.5*x(i)/vis)*exp(-0.25*(x(i)*x(i))/vis) + (0.5*((2*pi) - x(i))/vis)*exp(-0.25*((2*pi - x(i))^2)/vis);
end
ip(nx) = 1;
im(1) = nx;

for i = 1:nx
	u(i) = (-2.0*vis*(dphi(i)/phi(i))) + 4;
end


u_init = u;
for t = 1:nt
    
    for i = 1:nx
        
        % Lax-Wendroff two step
        u1 = 0.5*(u(i) + u(ip(i))) - 0.25*(dt/dx)*(u(ip(i))^2 - u(i)^2);
        u2 = 0.5*(u(im(i)) + u(i)) - 0.25*(dt/dx)*(u(i)^2 - u(im(i))^2);                              
        u_new(i) = u(i) - 0.5*(dt/dx)*(u1^2 - u2^2); 
        
        % Lax-Wendroff
        %u1 = (u(i) + (u(i) + u(ip(i)))*0.5)*0.5;
        %u2 = (u(i) + (u(im(i)) + u(i))*0.5)*0.5;
        %u_new(i) = u(i) - 0.25*(dt/dx)*(u(ip(i))^2 - u(im(i))^2) + 0.5*(dt^2/dx^2)*(u1*(u(ip(i))^2 - u(i)^2) - u2*(u(i)^2 - u(im(i))^2));
        
        % Lax–Friedrichs
        %u_new(i) = 0.5*(u(ip(i)) + u(im(i))) - 0.25*(dt/dx)*(u(ip(i))^2 - u(im(i))^2);
        
        % Mac Cormak
        us = u(i) - (dt/dx)*0.5*(u(ip(i))^2 - u(i)^2);
        u(i) = 0.5*(u(i) + us) - 0.5*(dt/dx)*0.5*(u(i)^2 - u(im(i))^2);
      
    end
	
    %u_new(1) = u(nx);
    %u_new(nx) = u(1);

    u = u_new;
	%plot(x, u);drawnow;
    
end
plot(x,u)
%for i = 1:nx
%    if (x(i) >= 0) && (x(i) <= 2+nt*dt)
%        u_final(i) = 1.0;
%    else
%        u_final(i) = 0;
%    end
%end

%plot(x, u_final)