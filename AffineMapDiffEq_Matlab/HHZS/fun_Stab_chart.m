function H=fun_Stab_chart(ax,par)
H=zeros(1,size(ax,2));
for k=1:size(ax,2)

    par.Kpphi=ax(1,k);
    par.Kdphi=ax(2,k);

    s0start=ones(par.Ndim,par.Nstep);
    s0flat=s0start(:);
    v0flat=LinMapsqueeze(s0flat,par);


    eigsN=1;
    d = eigs(@(s) LinMapsqueeze(s+s0flat,par)-v0flat,numel(s0start),eigsN,'largestabs',...
        StartVector=s0start(:),Tolerance=1e-4);%SubspaceDimension=eigsN+5
    H(1,k)=max(log(abs(d)));
end


end