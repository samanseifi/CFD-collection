function [ output_args ] = cavity_flow( nx, ny, nt, nit, c, rho, nu, dt )

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
            

end

