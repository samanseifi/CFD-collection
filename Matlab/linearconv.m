clear all
close all

nx = 40;
nt = 200;
l = 4.0;
dx = l/(nx - 1);
x = 0:dx:l;
c = 1.0;
sigma = 0.1;
dt = sigma*(dx/c);
u_init = zeros(1,nx);
u = zeros(1,nx);
u_new = zeros(1,nx);

% --------- Initial Profile --------- %
for i = 1:nx
	if (x(i) >= 0) && (x(i) <= 1.0)
        u_0(i) = sin(4.0*pi*x(i));
    else
        u_0(i) = 0.0;
    end
end
% ---------- First Time Step --------- %

for i = 2:nx-1
   u_1(i) = u_0(i) - c*dt*(u_0(i) - u_0(i-1))/dx;
end
u_1(1) = 0;
u_1(nx) = 0;
u = u_1;
for t = 1:nt
    for i = 2:nx-1
        u_new(i) = u_1(i) - c*dt*(u_0(i+1) - u_0(i-1))/(dx);
    end
    u_new(1) = 0.0;
    u_new(nx) = 0.0;
    u_1 = u_new;
    u_0 = u_1;
    
	plot(x, u_new);
    drawnow;
end