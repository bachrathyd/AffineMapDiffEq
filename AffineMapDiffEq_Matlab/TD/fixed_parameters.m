function [A,invMM,K0,C1]=fixed_parameters()
% fixed parameters for e-scooter (Xiaomi M365)
mch = 9.5;
mfork = 2.797;
mwheel1 = 2.894;
mwheel2 = 1.136;
p = 0.829;
b = 0.2418;
e = 0.0281;
epsilon = deg2rad(14.2);% ez nagyon lassú egyébként
R = 0.111;
xCGch = 0.4314;
xCGfork = 0.01; % optimum: 0.027 m
zCGch = 0.0782;
zCGfork = 0.271586;
Aw1 = 0.01134;
Bw1 = 0.00628;
Aw2 = 0.00536;
Bw2 = 0.0029;
Ach = 0.1014;
Bch = 0.4374;
Cch = 0.3561;
Dch = -0.1074;
Afork = 0.341171;
Bfork = 0.4018;
Cfork = 0.0840292;
Dfork = 0.151596;
v = 0;
g = 9.81;

% mass matrix, checked
M11 = (8*Ach+4*Afork+8*Aw1+8*Aw2+4*Cfork+4*b^2*mfork+e^2*mfork+8*mch*R^2+3*mfork*R^2+8*mwheel1*R^2+8*mwheel2*R^2+4*mfork*xCGfork^2+16*mch*R*zCGch+8*mch*zCGch^2+8*b*mfork*zCGfork+4*mfork*zCGfork^2+4*mfork*(3*b*R+e*xCGfork+3*R*zCGfork)*cos(epsilon)+4*(Afork-Cfork+mfork*(b^2+R^2-xCGfork^2+2*b*zCGfork+zCGfork^2))*cos(2*epsilon)+4*b*mfork*R*cos(3*epsilon)-4*e*mfork*xCGfork*cos(3*epsilon)+4*mfork*R*zCGfork*cos(3*epsilon)-e^2*mfork*cos(4*epsilon)+mfork*R^2*cos(4*epsilon)+4*b*e*mfork*sin(epsilon)+4*mfork*R*xCGfork*sin(epsilon)+4*e*mfork*zCGfork*sin(epsilon)-8*Dfork*sin(2*epsilon)+4*e*mfork*R*sin(2*epsilon)+8*b*mfork*xCGfork*sin(2*epsilon)+8*mfork*xCGfork*zCGfork*sin(2*epsilon)+4*b*e*mfork*sin(3*epsilon)+4*mfork*R*xCGfork*sin(3*epsilon)+4*e*mfork*zCGfork*sin(3*epsilon)+2*e*mfork*R*sin(4*epsilon))/8;
M12 = (cos(epsilon)^3*(2*Dfork*(e+p)+2*Dch*e*1/cos(epsilon)^2-2*e*mch*xCGch*(R+zCGch)*1/cos(epsilon)^2-2*Afork*p*tan(epsilon)-2*Aw1*p*tan(epsilon)+2*Afork*(e+p)*tan(epsilon)-2*Cfork*(e+p)*tan(epsilon)-2*mwheel1*p*R^2*1/cos(epsilon)^2*tan(epsilon)+4*Dfork*p*tan(epsilon)^2-2*Dfork*(e+p)*tan(epsilon)^2-2*Aw1*p*tan(epsilon)^3-2*Cfork*p*tan(epsilon)^3-mfork*(2*(e+p)*xCGfork+e*(e+2*p+e*cos(2*epsilon))*1/cos(epsilon)-2*e*R*sin(epsilon)-2*e*(b+zCGfork)*tan(epsilon)+2*p*xCGfork*tan(epsilon)^2)*(b+zCGfork+R*1/cos(epsilon)+xCGfork*tan(epsilon)+sin(epsilon)*(e-R*tan(epsilon)))))/(2*p);
M21 = M12;
M22 = (4*Aw2*e^2*cos(epsilon)^2+4*Cch*e^2*cos(epsilon)^2+4*e^2*mch*xCGch^2*cos(epsilon)^2+4*Aw1*(e+p)*cos(epsilon)^2*(p+e*cos(epsilon)^2)+4*Cfork*(e+p)*cos(epsilon)^2*(p+e*cos(epsilon)^2)+4*Dfork*e*(e+p)*cos(epsilon)^3*sin(epsilon)+4*mwheel1*p^2*R^2*sin(epsilon)^2+2*Aw1*p*(e+2*p+e*cos(2*epsilon))*sin(epsilon)^2+2*Cfork*p*(e+2*p+e*cos(2*epsilon))*sin(epsilon)^2+4*Dfork*e*p*cos(epsilon)*sin(epsilon)^3-Dfork*p*(e+2*p+e*cos(2*epsilon))*sin(2*epsilon)+Dfork*(e+p)*(e+2*p+e*cos(2*epsilon))*sin(2*epsilon)-Afork*e*p*sin(2*epsilon)^2-Aw1*e*p*sin(2*epsilon)^2+Afork*e*(e+p)*sin(2*epsilon)^2+Aw1*e*(e+p)*sin(2*epsilon)^2+mfork*(2*(e+p)*xCGfork*cos(epsilon)^2+2*p*xCGfork*sin(epsilon)^2+e*cos(epsilon)*(e+2*p+e*cos(2*epsilon)-2*(b+zCGfork)*sin(epsilon)-R*sin(2*epsilon)))^2)/(4*p^2);
MM = [M11 M12; M21 M22];
invMM = inv(MM);

% damping matrix, checked
C1 = [0 0;
    0 0];

% stiffness matrix, checked
K11 = -(g*(2*mch*R+mfork*R+2*mwheel1*R+2*mwheel2*R+2*mch*zCGch+2*mfork*(b+zCGfork)*cos(epsilon)+mfork*R*cos(2*epsilon)+2*mfork*xCGfork*sin(epsilon)+e*mfork*sin(2*epsilon)))/2;
K12 = -((-3*e^2*g*mfork*R-4*e*g*R*(mfork*p+mch*xCGch))*cos(epsilon)+R*(-2*e*g*mfork*xCGfork-4*g*mfork*p*xCGfork-2*e*g*mfork*xCGfork*cos(2*epsilon)-e^2*g*mfork*cos(3*epsilon)+e*g*mfork*R*sin(epsilon)-4*g*mwheel1*p*R*sin(epsilon)+2*b*e*g*mfork*sin(2*epsilon)+2*e*g*mfork*zCGfork*sin(2*epsilon)+e*g*mfork*R*sin(3*epsilon)))/(4*p*R);
K21 = K12;
K22 = -(-(e*g*mfork*p*R^2)+4*g*mwheel1*p^2*R^2-2*mfork*R*(b*e*g*p+e*g*p*zCGfork)*cos(epsilon)-4*g*mwheel1*p^2*R^2*cos(2*epsilon)+2*b*e*g*mfork*p*R*cos(3*epsilon)+2*e*g*mfork*p*R*zCGfork*cos(3*epsilon)+e*g*mfork*p*R^2*cos(4*epsilon)+2*e*g*mfork*p*R*xCGfork*sin(epsilon)+8*g*mfork*p^2*R*xCGfork*sin(epsilon)+2*e^2*g*mfork*p*R*sin(2*epsilon)+4*e*g*mfork*p^2*R*sin(2*epsilon)+4*e*g*mch*p*R*xCGch*sin(2*epsilon)+2*e*g*mfork*p*R*xCGfork*sin(3*epsilon)+e^2*g*mfork*p*R*sin(4*epsilon))/(8*p^2*R);
K0 = [K11 K12; K21 K22];

% set of equations of motion
A = zeros(4,4);
A(1:2,3:4)= eye(2,2);
A(3:4,1:2) = -invMM*K0;
A(3:4,3:4) = -invMM*C1;
