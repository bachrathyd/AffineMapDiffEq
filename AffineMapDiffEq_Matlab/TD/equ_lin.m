function yprime = equ_lin(t,y,par)
% Single wheel (no rider) with brush tire model - LINEARIZED
% Right hand side with finite difference discretization of the PDE

Omega=par.Omega;
k=par.k;
m=par.m;
RE=par.RE;
g=par.g;
x=par.x;
dx=par.dx;
n=par.n;

sigma1 = y(1);
sigma2 = y(2);
sigma3 = y(3);
theta = y(4);
q = y(5:n+5);

% =======================================
% discretized PDE / tire deformation
qd = zeros(n+1,1);

xd = - RE*Omega;
qd(1) = -sigma1 - x(1)*sigma3 - xd*(q(2)-q(1))/dx;
for i=2:n
    qd(i) = -sigma1 - x(i)*sigma3 - xd*(q(i+1)-q(i-1))/(2*dx);
end
qd(n+1) = 0; % BC for brush model / no deformation at the leading edge

Fy = 0;
Mz = 0;
for i=2:n+1
    Fy = Fy + k * (q(i)+q(i-1))/2 * dx;
    Mz = Mz + k * (x(i)+x(i-1))/2 * (q(i)+q(i-1))/2*dx;
end

% =======================================
% ODEs / wheel dynamics

f1 = (5*Fy)/m + Omega*RE*sigma3 + 4*g*theta;
f2 = (4*Fy + 2*m*Omega*RE*sigma3 + 4*g*m*theta)/(m*RE);
f3 = (4*Mz)/(m*RE^2) - 2*Omega*sigma2;

% =======================================
% Composing RHS
 
yprime = [f1; f2; f3; sigma2; qd];


