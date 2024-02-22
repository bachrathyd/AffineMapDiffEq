
%% Parameters for the diff.eq.
par=[];
par.zeta = 0.03;
par.delta  = 1.1;%0.8;1.1;
par.eps = 0.9;
par.tau = 2*pi;
par.b = 0.05;

%% paramteres for the spactrum calculation
par.T = 2.0*pi; % Pariod of Mapping
par.diffun=@DelayMathieu; % The function handle for the diff.eq.
par.Ndim=2; % Dimension of the system

%defining the sampling of the delay (state)
taumax=2*pi; %max delay - it must be larger or equal to the maximal delyay (important if it is not constant)

par.StateSmapled=linspace(-taumax,0,100);%no-need for uniform sampling
par.Nstep=length(par.StateSmapled);

par.doplot=false;%true;%false % plot the mapped resilts - for debuging

%% testing the mapping
s0start=ones(par.Ndim,par.Nstep);
v0start=LinMap(s0start,par);
%Xout=LinMap(Xout,par)



s0flat=s0start(:);
v0flat=LinMapsqueeze(s0flat,par);

%%testing the spectrum & the fixed point
eigsN=20;
tic
%d = eigs(@(s) LinMapsqueeze(s,par),numel(s0start),eigsN,'largestabs');%Linear
[V,D] = eigs(@(s) LinMapsqueeze(s+s0flat,par)-v0flat,numel(s0start),eigsN,'largestabs');%Affine
toc

d=diag(D)
for k=1:40%fixpoint iteration s0 will (should) converge to the fixpoint
    %convergence speed depends on the abs. val. of largest mu, which is not computed
    figure(44),clf
    par.doplot=true;
    s0flat =real( v0flat - V * ( (V'*V) \ (V' * (v0flat-s0flat)) .* ((d) ./ (d - 1.0))) ) ;% imaginary part is created only by numerical errors
    v0flat=LinMapsqueeze(s0flat,par);
    errorestimation=norm(v0flat-s0flat)/numel(s0flat)
    if (errorestimation)<1e-10
        a0flat=s0flat;
        break
    end
    drawnow
end
a0flat=s0flat;
a0=reshape(a0flat,par.Ndim,[]);


figure(2345)
plot(par.StateSmapled,a0)

figure(1),clf
subplot(2,2,1)
plot(log(abs(d))),grid on

fi=linspace(0,2*pi,100);
subplot(2,2,2)
plot(sin(fi),cos(fi),'k-')
hold on
plot(d,'ro'),grid on

subplot(2,2,3)
plot(par.StateSmapled,a0)

%long simulation from the fix point s0
subplot(2,2,4)
par.T = 2.0*pi*20; % Pariod of Mapping
par.doplot=true;
a0_long=LinMap(a0,par);
%


%% -------- Stability Chart with MDBM-------------

addpath('C:\Users\Bacharthy\Documents\GitHub\MBDM\code_folder')
par.T = 2.0*pi; % Pariod of Mapping
par.b = 0.005;
par.zeta=0.0001;
par.doplot=false;
par.Nstep=100;
ax=[];
mdbm_options=mdbmset('timelimit',Inf);

ax(1).val=linspace(-1,6,8);%delata
ax(2).val=linspace(-0.1,8,8);%eps
mdbm_sol=mdbm(ax,'fun_Stab_chart',3,mdbm_options,par);


figure(3)
plot_mdbm(mdbm_sol,'k');
view(2)
