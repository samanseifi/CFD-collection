using PyPlot

l = 2.0
nx, ny, nt = 40, 40, 10;
dt = 0.01;
dx = l/(nx - 1);
dy = l/(ny - 1);
x = 0:dx:l
y = 0:dy:l

vis = 0.1

u = zeros(nx, ny);
u_inital = zeros(nx, ny);
u_new = zeros(nx, ny);

for i = 1:nx
	for j = 1:ny
		if x[i] >= 0.5 && x[i] <= 1.0 && y[j] >= 0.5 && y[j] <= 1.0
			u_init[i, j] = 2.0;
		else
			u_init[i, j] = 1.0;
		end
	end
end

u_ = u_init;
contourf(x, y , u)
for t = 1:nt
	for i = 2:nx-1
		for j = 2:ny-1
			u_new[i, j] = u[i, j] - (u[i, j]*dt*(u[i, j] - u[i-1, j])/dx) + (vis*dt*(u[i+1, j] - 2*u[i, j] + u[i-1, j])/(dx*dx)) - (u[i, j]*dt*(u[i, j] - u[i, j-1])/dy) + (vis*dt*(u[i, j+1] - 2*u[i, j] + u[i, j-1])/(dy*dy));
		end	
	end
	u_new[1,:] = 1.0;
	u_new[nx,:] = 1.0;
	u_new[:,ny] = 1.0;
	u_new[:,1] = 1.0;
	u = u_new;
	contourf(x, y, u)
	#plot(x,u(:,7))
	draw()
	sleep(0.5)
end
