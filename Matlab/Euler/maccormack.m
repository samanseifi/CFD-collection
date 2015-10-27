hold on
clear
Input; % Input Values
Sod; % Sod shock initial condition
 
u = m./rho;
p = (gamma-1)*(e - 0.5*rho.*u.^2);
a = sqrt((gamma*p)./rho);
dt = CFL * dx/max(a);
 
%Time integration
for t = dt:dt:t_final
    Q = [rho; rho.*u; e]; 
    E = [rho.*u; rho.*u.^2+p; u.*(e+p)]; 
    
    %-------Predictor Step------%
    Q_bar(1:3,1:n-1) = 0.5*(Q(1:3,1:n-1) + Q(1:3,2:n)) - dt/(2*dx)*(E(1:3,2:n) - E(1:3,1:n-1));
    
    %Q_bar(1:3, 1:n-1) = Q(1:3, 1:n-1) - (dt/dx)*(E(1:3, 2:n) - E(1:3, 1:n-1)); 
    
    %update the answers for Q_bar
    rho(1:n-1) = Q_bar(1,1:n-1); 
    u = Q_bar(2,1:n-1)./rho(1:n-1); 
    e = Q_bar(3,1:n-1); 
    p = (gamma-1)*(e(1:n-1) - 0.5*rho(1:n-1).*u(1:n-1).^2);
    
    % calculate E_bar
    E_bar(1:3,1:n-1)=[rho(1:n-1).*u(1:n-1); rho(1:n-1).*u(1:n-1).^2+p(1:n-1); u(1:n-1).*(e(1:n-1)+p(1:n-1))];
    
        
    % Calculate viscous official parameter
    alpha = 0.4;
    A = [zeros(1,n-1);ones(1,n-1);u(1:n-1)];
    
    for i=1:3
        Artificial(i,2:n-1) = alpha*dx^2*rho(2:n-1).*abs((u(2:n-1)-u(1:n-2))/dx).*((u(2:n-1)-u(1:n-2))/dx).*A(i,2:n-1);
    end
    
    %--------Corrector Step------%
    %E_bar = E_bar-Artificial;
    Q(1:3,2:n-1) = Q(1:3,2:n-1) - (dt/dx)*(E_bar(1:3,2:n-1) - E_bar(1:3,1:n-2));
    
    %Q(1:3, 2:n-1) = 0.5.*(Q(1:3, 2:n-1) + Q_bar(1:3, 2:n-1) - (dt/dx)*(E_bar(1:3, 2:n-1) - E_bar(1:3, 1:n-2)));
    
    % Update the answers
    rho=Q(1,1:n); 
    u=Q(2,1:n)./rho(1:n); 
    e=Q(3,1:n); 
    p=(gamma-1)*(e-0.5*rho.*u.^2);
end
% hold on
% plot(x,rho,'g');
% plot(x,rho.*u,'g');
% plot(x,e,'g');