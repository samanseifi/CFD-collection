nx = 81;
nt = 80;
l = 4.0;
dx = l/(nx - 1);
x = 0:dx:l;
c = 1.0;
sigma = 0.8;
dt = sigma*(dx/c);
u_init = zeros(nx);
u = zeros(nx);
u_new = zeros(nx);

for i = 1:nx
	u_init(i) = 0.0;
end
u_init(1) = 1.0;

u = u_init;
for t = 1:nt
    for i = 2:nx-1
        u_new(i) = u(i) - c*dt*(u(i) - u(i-1))/dx
    end
	u_new(1) = 1.0;

    u = u_new;
	plot(x, u);
	drawnow;

end