clear all
close all

Input;

%nx = 201;
%nt = 50;
%l = 4.0;
%dx = l/(nx - 1);
%x = 0:dx:l;
c = 1.0;

sigma = 0.8;
dt = sigma*(dx/c);

u = Step(x,nx);
u = WavePacket(x, nx, 16.0);

for i = 1:nx
    ip(i) = i+1;
    im(i) = i-1;
end
ip(nx) = 1;
im(1) = nx;


for t = 1:nt
    for i = 1:nx
        
        u_new(i) = 0.5*(u(ip(i)) + u(im(i))) - (dt/dx)*0.5*(u(ip(i)) - u(im(i)));
    
    
    end
    
    
    u_new(1) = u(1);
    u_new(end) = u(end);
    u = u_new;
	plot(x, u);
	drawnow;
    
end
%hold on
%for i = 1:nx
%    if (x(i) >= 0) && (x(i) <= 2+nt*dt)
%        u_final(i) = 1.0;
%    else
%        u_final(i) = 0;
%    end
%end

%plot(x, u)