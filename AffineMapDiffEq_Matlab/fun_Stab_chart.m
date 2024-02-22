function H=fun_Stab_chart(ax,par)
H=zeros(1,size(ax,2));
for k=1:size(ax,2)
    par.delta=ax(1,k);
    par.eps=ax(2,k);


    s0start=ones(par.Ndim,par.Nstep);
    s0flat=s0start(:);
    v0flat=LinMapsqueeze(s0flat,par);


    eigsN=3;
    d = eigs(@(s) LinMapsqueeze(s+s0flat,par)-v0flat,numel(s0start),eigsN,'largestabs',...
        StartVector=s0start(:),SubspaceDimension=eigsN+4,Tolerance=1e-5);%SubspaceDimension=eigsN+5
    H(1,k)=max(abs(d)-1);
end


end