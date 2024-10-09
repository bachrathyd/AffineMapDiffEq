function yprime = equ(t,y,par)
% Single wheel (no rider) with brush tire model
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
xd = zeros(n+1,1);
qd = zeros(n+1,1);

xd = - RE*Omega + q./cos(theta)^2*sigma3;
qd(1) = -sigma1 - q(1)*tan(theta)*sigma2 - x(1)*sigma3 - xd(1)*(q(2)-q(1))/dx;
for i=2:n
    qd(i) = -sigma1 - q(i)*tan(theta)*sigma2 - x(i)*sigma3 - xd(i)*(q(i+1)-q(i-1))/(2*dx);
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

f1 = (5*Fy + m*Omega*RE*sigma3 + 4*g*m*sin(theta) + m*sigma1*sigma2*1/cos(theta)^3*(-4*sin(theta) + sin(3*theta)) - 5*m*RE*sigma2^2*tan(theta) - 5*m*Omega*RE*sigma3*tan(theta)^2 - 5*m*RE*sigma3^2*tan(theta)^3)/(m*(-4 + 5*1/cos(theta)^2));
f2 = -((2*Fy + 2*m*Omega*RE*sigma3 + 2*Fy*cos(2*theta) + 4*g*m*sin(theta) - 2*m*RE*sigma2^2*sin(2*theta) + m*RE*sigma3^2*tan(theta))/(m*RE*(-3 + 2*cos(2*theta))));
f3 = (-(m*RE^2*sigma2*(8*Omega + sigma3*1/cos(theta)^3*(29*sin(theta) + 5*sin(3*theta)))) + 16*(Mz + m*RE*sigma1*sigma3*1/cos(theta)^2*tan(theta)))/(4*RE^2*(m + 6*m*tan(theta)^2));

% =======================================
% Composing RHS
 
yprime = [f1; f2; f3; sigma2; qd];


