% Osher scheme (O-variant) with MUSCL for one-dimensional Euler equations
% Runge-Kutta time stepping

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory in Section 10.8 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program makes Figs. 10.32 -- 10.35 in the book

% Functions called: extrap,  osheulerstep,  problem_specification, Riemann

clear all
global  PRL  CRL MACHLEFT  gamma  pleft  pright  rholeft  rhoright  uleft...
	uright  tend  epsi  invariants... 
	lambda		% lambda = dt/dx

		% .....................Input............................
gamma = 1.4; 	% Ratio of specific heats
J = 48;		% Number of grid cells
flops = 0;
invariants = 0;	% Enter 0 for extrapolation of primitive variables or
		%   something else for extrapolation of Riemann invariants 
limtype = 1;	% Enter 1 for minmod or 2 for van Albada or something else
		%	if MUSCL not wanted

gammab = 1/(gamma - 1); gam1 = gamma-1; gamgam = gamma^gamma;
problem_specification	

h = 1/J;  				% Cell size  
dt = lambda*h;				% Time step
n = floor(tend/dt);			% Number of time-steps

% 		Definition of grid numbering 
%       x=0    					 x=1
% grid   |---o---|---o---|---o---  ...  --|---o---|
%            1   1   2   2   3           J-1  J 

xcenter = h*[1:J] - h/2;		% Location of cell centers

press = zeros(size(xcenter));		% Preallocation of pressure, 
rhoold = press; uold = press;		%       density, velocity,
rhonew = press; mnew = press;		%	momentum, 
totenew = press;			%	total energy

for j = 1:length(xcenter)		% Initial conditions
  if xcenter(j) < 0.5, press(j) = pleft; rhoold(j) = rholeft;  uold(j) = uleft;
  else,		      press(j) = pright; rhoold(j) = rhoright; uold(j) = uright;
  end
end

	% toten = rho*(total energy)
totenold = rhoold.*(0.5*uold.*uold + gammab*press./rhoold);
totenleft = totenold(1); totenright = totenold(J);
mold = rhoold.*uold;			% m is momentum
c = sqrt(gamma*press./rhoold);		% Sound speed

		% Preallocation of auxiliary variables
aalpha = zeros(J-1,1); c13 = aalpha; uH = aalpha; p13 = aalpha;
% 
count = 0;	% Counter of flux evaluations

% Preallocation of extrapolated states: Riemann invariants and primitive var.
Z1L = zeros(J-1,1); Z2L = Z1L; Z3L = Z1L;
Z1R = Z1L; Z2R = Z1L; Z3R = Z1L;
U1L = Z1L; U1R = Z1L; U2L = Z1L; U2R= Z1L;
U3R = Z1L; U3L = Z1L;

Z1 = zeros(J,1); Z2 = Z1; Z3 = Z1;	% Preallocation of Riemann invariants
flopcount = flops;			% Operations counter

t = 0;
for i = 1:n,  t = t + dt;
  rhostar = rhoold; mstar = mold; totenstar = totenold;
  rkalpha = 0.25;  osheulerstep
  rkalpha = 1/3;   osheulerstep
  rkalpha = 0.5;   osheulerstep
  rkalpha = 1;     osheulerstep
  
  rhonew = rhostar; mnew = mstar; totenew = totenstar;
  uold = mnew./rhonew; press = gam1*(totenew - 0.5*mnew.*uold);
  rhoold = rhonew; totenold = totenew; mold = mnew;
  c = sqrt(gamma*press./rhoold); 
end

flopcount = flops - flopcount;
mach = uold./c;  entropy = log(press./rhoold.^gamma);

figure(1), clf
subplot(2,3,1),hold on,title('DENSITY','fontsize',14),plot(xcenter,rhonew,'o')
subplot(2,3,2),hold on,title('VELOCITY','fontsize',14),plot(xcenter,uold,'o')
subplot(2,3,3),hold on,title('PRESSURE','fontsize',14),plot(xcenter,press,'o')
subplot(2,3,4),hold on,title('MACHNUMBER','fontsize',14),plot(xcenter,mach,'o')
subplot(2,3,5),hold on,title('ENTROPY','fontsize',14),plot(xcenter,entropy,'o')
subplot(2,3,6), axis('off'), hold on
	title('Osher scheme, O-variant','fontsize',14)
text(0,0.9,['lambda = ', num2str(lambda),'  t = ',num2str(n*dt)])
if invariants == 0,  s1 = 'Primitive extrapolation';
else,  		     s1 = 'Riemann extrapolation';
end
text(0,0.75,s1),	text(0,0.6,'SHK Runge-Kutta')
if limtype == 0,  	s1 = 'No MUSCL';
elseif limtype == 1,  	s1 = 'MUSCL, minmod limiter';
else,  			s1 = 'MUSCL, van Albada limiter';
end
text(0,0.45,s1)
text(0,0.3,['flopcount = ',num2str(flopcount)])

Riemann		% Plot exact solution

