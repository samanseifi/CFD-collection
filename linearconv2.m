%clear all
%close all
nx = 201;
nt = 50;
l = 4.0;
dx = l/(nx - 1);
x = 0:dx:l;
c = 1.0;

sigma = 0.8;
dt = sigma*(dx/c);
u_init = zeros(1,nx);
u = zeros(1,nx);
u_new = zeros(1,nx);
for i = 1:nx
	if (x(i) >= 0) && (x(i) <= 1.0)
        u_init(i) = sin(16.0*pi*x(i));
    else
        u_init(i) = 0.0;
    end
end
for i = 1:nx
    if (x(i) >= 0) && (x(i) <= 2)
        u_init(i) = 1.0;
    else
        u_init(i) = 0;
    end
end

u = u_init;
for t = 1:nt
    for i = 2:nx-1
        % --- Step One --- %
        u1(i) = 0.5*(u(i) + u(i+1)) - 0.5*c*(dt/dx)*(u(i+1) - u(i));
        u2(i) = 0.5*(u(i-1) + u(i)) - 0.5*c*(dt/dx)*(u(i) - u(i-1));
        u_new(i) = u(i) - 0.5*c*(dt/dx)*(u(i+1) - u(i-1)) + 0.5*c*c*(dt^2/dx^2)*(u(i+1) - 2*u(i) + u(i-1)); % Lax-Wendroff
        %u_new(i) =0.5*(u(i+1) + u(i-1))  - c*dt*(u(i+1) - u(i-1))/(2*dx);
        % --- second step --- %
        u_new(i) = u(i) - c*(dt/dx)*(u1(i) - u2(i));
    end
	u_new(1) = 1.0;
    u_new(nx) = 0.0;

    u = u_new;
	plot(x, u);
	drawnow;
    
end
hold on
for i = 1:nx
    if (x(i) >= 0) && (x(i) <= 2+nt*dt)
        u_final(i) = 1.0;
    else
        u_final(i) = 0;
    end
end

plot(x, u)