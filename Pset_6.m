%AERSP 458 
clear, clc

format long
LU = 3.844e05; 
rho = 5000/LU; %[LU]
lambda = 0.012151;
C = -1.59;

x0 = @(beta) lambda-1+rho*cos(beta);
y0 = @(beta) rho*sin(beta);
xdot0 = @(alpha,v0) v0*cos(alpha);
ydot0 = @(alpha,v0) v0*sin(alpha);
r1 = @(beta) sqrt((x0(beta)-lambda)^2+y0(beta)^2);
r2 = @(beta) sqrt((x0(beta)+(1-lambda))^2+y0(beta)^2);

tspan = 0:0.1:8*pi;
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
%Moon Location
Moon = zeros(length(tspan),2);
Moon(:,1) = lambda-1;
%Earth Location
Earth = zeros(length(tspan),2);
Earth(:,1) = lambda; 
%%
%part (a)
beta = 0; 
alpha = beta + pi/2;
%IC 
v0 = sqrt(2*(lambda/r2(beta)+(1-lambda)/r1(beta)+0.5*(x0(beta)^2+y0(beta)^2)+C));
K0 = [x0(beta); xdot0(alpha, v0); y0(beta); ydot0(alpha,v0)]; %x,xdot,y,ydot

[t,k] = ode45('StateSpace',tspan, K0, options);

figure(1);
plot(t,k(:,1),t,k(:,3),t,Earth(:,1),t,Earth(:,2),t,Moon(:,1),t,Moon(:,2),'linewidth',2);
legend('Satellite x','Satellite y','Earth x','Earth y','Moon x','Moon y');
%%
%part (b)
beta = pi/2; 
alpha = beta + pi/2;
%IC 
v0 = sqrt(2*(lambda/r2(beta)+(1-lambda)/r1(beta)+0.5*(x0(beta)^2+y0(beta)^2)+C));
K0 = [x0(beta); xdot0(alpha, v0); y0(beta); ydot0(alpha,v0)]; %x,xdot,y,ydot

[t,k] = ode45('StateSpaceb',tspan, K0, options);
figure(2);
plot(t,k(:,1),t,k(:,3),t,Earth(:,1),t,Earth(:,2),t,Moon(:,1),t,Moon(:,2),'linewidth',2);
legend('Satellite x','Satellite y','Earth x','Earth y','Moon x','Moon y');
%%
%part (c)
beta = 0.625*pi; 
alpha = beta + pi/2;
%IC 
v0 = sqrt(2*(lambda/r2(beta)+(1-lambda)/r1(beta)+0.5*(x0(beta)^2+y0(beta)^2)+C));
K0 = [x0(beta); xdot0(alpha, v0); y0(beta); ydot0(alpha,v0)]; %x,xdot,y,ydot

[t,k] = ode45('StateSpaceC',tspan, K0, options);
figure(3);
plot(t,k(:,1),t,k(:,3),t,Earth(:,1),t,Earth(:,2),t,Moon(:,1),t,Moon(:,2),'linewidth',2);
legend('Satellite x','Satellite y','Earth x','Earth y','Moon x','Moon y');

%% Problem 2 
clear, clc
lambda = 0.012151;
%Earth-Moon L2
xStar = -1.1556;
yStar = 0;
r1Star = @(lam) 1.1556 + lam;
r2Star = @(lam) 1.1556 - (1-lam);
%dxCoef = @(lam) (1-lam)*(1/r1Star^3-3*(xStar-lam)^2/(r1Star)^5)+lam*(1/r2Star^3-3*(xStar+1-lam)^2/r2Star^5);
%dyCoef = @(lam) (1-lam)*(1/r1Star^3)+lam*(1/r2Star^3);

A = [0 1 0 0; 7.3888 0 0 2; 0 0 0 1; 0 -2 -2.1944 0];
[V,E] = eig(A);











