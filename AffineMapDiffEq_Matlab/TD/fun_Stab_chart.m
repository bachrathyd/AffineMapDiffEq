function H=fun_Stab_chart(ax,par)
H=zeros(1,size(ax,2));
for k=1:size(ax,2)
display(k/size(ax,2))
    par.Omega=ax(1,k);
    par.a=ax(2,k);

    s0start=ones(par.Ndim,par.Nstep);
    s0flat=s0start(:);
    v0flat=LinMapsqueeze(s0flat,par);


    eigsN=6;

opts=[];
opts.p=30;%30;
opts.tol = (1e-10);%15;
opts.disp=5;
opts.issym=false;
opts.v0=ones(par.Ndim*par.Nstep,1);

opts.maxit=30;
opts.v0=par.s0flat;
opts.fail='keep';
    mus = eigs(@(s) LinMapsqueeze(s+s0flat,par)-v0flat,numel(s0start),eigsN,'largestabs',...
        StartVector=s0start(:),Tolerance=1e-4);%SubspaceDimension=eigsN+5
    %H(1,k)=max(log(abs(d)));
    H(1,k)=max(real(log(mus)/par.T));
end


end