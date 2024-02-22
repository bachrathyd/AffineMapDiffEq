function v=LinMap(s,par)
%TODO: It is redundant, I use an interpolant to find v, then, I redu the
%interpolant again based on v - however it is safe in this way

hist = griddedInterpolant(par.StateSmapled,s','cubic');%faster%Note, it works with transposed
%hist = griddedInterpolant(par.StateSmapled,s','linear');%slower%Note, it works with transposed
par.h=hist;

%options = ddeset('AbsTol',1e-8,'RelTol',1e-8);
options = ddeset();
sol = dde23(@DelayMathieu,[par.tau],@historyInterp,[0, par.T],options,par);

histFinal = griddedInterpolant(sol.x,sol.y','spline');%faster
%histFinal = griddedInterpolant(sol.x,sol.y','linear');%slower

newStateSmapled=par.StateSmapled+par.T;
overlap=newStateSmapled<0.0;%necessary if T<tau
v=[historyInterp(newStateSmapled(overlap),par),histFinal(newStateSmapled(~overlap))'];


if par.doplot
    plot(par.StateSmapled,s,'--',LineWidth=2),hold on
    plot(sol.x,sol.y)
    plot(par.StateSmapled+par.T,v,'--',LineWidth=2),hold on
    %pause
end

end


