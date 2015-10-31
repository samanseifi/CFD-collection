l = 4.0; 
% l = 2*pi %nonlinear  convec

nx = 201;
ny = 201;

nt = 100;

dx = l/(nx - 1);
dy = l/(ny - 1);

x = 0:dx:l;

c = 1.0;
dt = 0.0016;
