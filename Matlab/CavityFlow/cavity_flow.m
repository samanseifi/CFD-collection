function [ u,v,p ] = cavity_flow( nx, ny, nt, nit, c, rho, nu, dt )

    x = linspace(0,2,nx);
    y = linspace(0,2,ny);
    [X,Y] = meshgrid(x,y);
    
    u = zeros(nx,ny); v = zeros(nx,ny);
    p = zeros(nx,ny); b = zeros(nx,ny);
    
    un = 0;
    vn = 0;
    b = zeros(nx, ny);
    
    for i = 1:nt
        un = u;
        vn = v;
        
        b = buildUpB(b, rho, dt, u, v, dx ,dy);
        p = pressPoisson(p, dx, dy, b);
        
        u(2:end-1, 2:end-1) = un(2:end-1, 2:end-1) - ...
            un(2:end-1, 2:end-1)*(dt/dx).*(un(2:end-1, 2:end-1) - un(2:end-1, 1:end-2)) - ...
            vn(2:end-1, 2:end-1)*(dt/dy).*(un(2:end-1, 2:end-1) - un(1:end-2, 2:end-1)) - ...
            dt/(2*rho*dx)*(p(2:end-1, 3:end) - p(2:end-1, 1:end-2)) + ...
            nu*(dt/dx^2*(un(2:end-1, 3:end) - 2*un(2:end-1, 2:end-1) + un(1:end-2, 2:end-1)));
        
        v(2:end-1, 2:end-1) = vn(2:end-1, 2:end-1) - ...
            un(2:end-1, 2:end-1)*(dt/dx).*(vn(2:end-1, 2:end-1) - vn(2:end-1, 1:end-2)) - ...
            vn(2:end-1, 2:end-2)*(dt/dy).*(vn(2:end-1, 2:end-1) - vn(1:end-2, 2:end-1)) - ...
            dt/(2*rho*dy)*(p(3:end, 2:end-1) - p(1:end-2, 2:end-1)) + ...
            nu*(dt/dx^2*(vn(2:end-1, 3:end) - 2*vn(2:end-1, 2:end-1) + vn(2:end-1, 1:end-2)) + ...
            (dt/dy^2*(vn(3:end, 2:end-1) - 2*vn(2:end-1, 2:end-1) + vn(1:end-2, 2:end-1)));
        
        u(1,:) = 0;
        u(:,1) = 0;
        u(:, end-1) = 0;
        u(end-1, :) = 1;
                
        v(1,:) = 0;
        v(:,1) = 0;
        v(:, end-1) = 0;
        v(end-1, :) = 1;

    end
end

