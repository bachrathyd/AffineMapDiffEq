
function dydt = motor_rhs_control_linear_hierarchical_multipledelays(t,y,Z,par)

Kpphi = par.Kpphi;
Kdphi = par.Kdphi;


% history of y
ylag1 = Z(:,1); % history of delta
ylag2 = Z(:,2); % history of phi

Kpdelta = 10;
Kddelta = -5;

invMM=par.invMM;
%K0=par.K0;
%C1=par.C1;
A=par.A;

% control matrices
Pphi = [0 0; -Kpphi*Kpdelta 0];
Pdelta = [0 0; 0 -Kpdelta];
Dphi = [0 0; -Kdphi*Kpdelta 0];
Ddelta = [0 0; 0 -Kddelta];


B1 = zeros(4,4); % not the same B as in semi-discretization!
B1(3:4,1:2) = invMM*Pdelta;
B1(3:4,3:4) = invMM*Ddelta;

B2 = zeros(4,4); % not the same B as in semi-discretization!
B2(3:4,1:2) = invMM*Pphi;
B2(3:4,3:4) = invMM*Dphi;

dydt = A*y + B1*[0;ylag1(2);0;ylag1(4)] + B2*[ylag2(1);0;ylag2(3);0]; % time delay in delta and phi!

end