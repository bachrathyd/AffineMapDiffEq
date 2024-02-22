function dydt = DelayMathieu(t,y,Z,par)
% Differential equations function for DDEX1
dydt=zeros(2,1);
dydt(1) =  y(2);
dydt(2) = 1.0+-(par.delta + par.eps * cos(t))* y(1)  -2 * par.zeta * y(2) + par.b * Z(1);
end