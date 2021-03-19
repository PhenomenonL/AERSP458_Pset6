function Kdot = StateSpace(t,K)
%STATESPACE output state space representation of CR3BP EOMs
%state space variable K
Kdot = zeros(4,1);
lambda = 0.012151;
LD = 3.844e05; 
rho = 5000/LD; %[LD]
%part (a)
beta = 0; 

x0 = lambda-1+rho*cos(beta);
y0 = rho*sin(beta);

r1 = sqrt((x0-lambda)^2+y0^2);
r2 = sqrt((x0+(1-lambda))^2+y0^2);

Kdot(1) = K(2);
Kdot(2) = 2*K(4)+K(1) - (1-lambda)*(K(1)-lambda)/r1^3 - lambda*(K(1)+1-lambda)/r2^3;
Kdot(3) = K(4);
Kdot(4) = -2*K(2) + K(3) - (1-lambda)*K(3)/r1^3 - lambda*K(3)/r2^2;
return

