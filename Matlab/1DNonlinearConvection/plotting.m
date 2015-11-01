clear all
close all

l = 2*pi;
nx = 101;
nt = 200;
dt = 0.0016;

% Plotting Non-linear convection
[x1, u1] = nonlinearconv('LaxWendroffTwoStep', l, nx, nt, dt);
[x2, u2] = nonlinearconv('LaxWendroff', l, nx, nt, dt);
[x3, u3] = nonlinearconv('LaxFriedrichs', l, nx, nt, dt);
[x4, u4] = nonlinearconv('MacCormak', l, nx, nt, dt);

x_init = x1;
u_init = DoubleGaussian(x_init, nx, 0.1); % vis = 0.1;
%u_init = Step(x_init, nx);
%u_init = WavePacket(x_init, 201, 16);

xmax = max(max([x1; x2; x3; x4]));
xmin = min(min([x1; x2; x3; x4]));
ymax = max(max([u1; u2; u3; u4; u_init]));
ymin = min(min([u1; u2; u3; u4; u_init]));

plot(x_init, u_init, '--k', x1, u1, '-o', x2, u2, '-s', x3, u3, '-^', x4, u4, '-*')
title('1D Non-linear convection');
xlabel('Spatial co-ordinate: $x$', 'Interpreter', 'latex');
ylabel('Transport property: $u$', 'Interpreter', 'latex');


axis([xmin xmax ymin-0.5 ymax+0.5])

legend('Initial', 'LaxWendroffTwoStep', 'LaxWendroff', 'LaxFriedrichs', 'MacCormak', 'Location', 'northwest')

