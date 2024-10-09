
%% Parameters for the diff.eq.
par = [];
par.Kpphi = -2000;
par.Kdphi = -100;
[par.A,par.invMM,par.K0,par.C1]=fixed_parameters();


tau1 = 0.01; % for steering (delta)
tau2 = 0.01; % for lean (phi)

%% paramteres for the spactrum calculation
par.T = max(tau1,tau2); % Pariod of Mapping
par.diffun=@motor_rhs_control_linear_hierarchical_multipledelays; % The function handle for the diff.eq.
par.Ndim=4; % Dimension of the system
par.options = ddeset('MaxStep',0.0001,'RelTol',1e-5,'AbsTol',1e-6); %options of the dde solver
par.lags=[tau1,tau2];%lags used in the dde solver


%defining the sampling of the delay (states)
taumax=max(tau1,tau2); %max delay - it must be larger or equal to the maximal delyay (important if it is not constant)

par.StateSmapled=linspace(-taumax,0,100);%no-need for uniform sampling
par.Nstep=length(par.StateSmapled);

par.doplot=false;%true;%false % plot the mapped results - for debuging only

%% testing the mapping- for the matrix form
s0start=zeros(par.Ndim,par.Nstep); %necessary initalization for the fixpoint and the eigen value calculations
% s0start=ones(par.Ndim,par.Nstep);
% s0start=rand(par.Ndim,par.Nstep);
v0start=LinMap(s0start,par);
%Xout=LinMap(Xout,par)


% testing the mapping- for the vectorized form
s0flat=s0start(:);
v0flat=LinMapsqueeze(s0flat,par);

%%testing the spectrum & the fixed point
eigsN=20;
Norm_of_perturbation=1e-5;
tic
[V,D] = eigs(@(s) (LinMapsqueeze(s*Norm_of_perturbation+s0flat,par)-v0flat)/Norm_of_perturbation,numel(s0start),eigsN,'largestabs');%Affine
toc
% plotting the specttral properties
d=diag(D);


figure(1)%,clf
subplot(2,2,1)
plot(d,'gs','LineWidth',2,'MarkerSize',7), hold on
fi=linspace(0,2*pi,1000);
plot(sin(fi),cos(fi)),grid on
title('mu')
subplot(2,2,2)
lam_sim=log(d)/par.T;

plot(lam_sim,'gs','LineWidth',2,'MarkerSize',7), hold on,grid on
for k=-1:1
 plot(lam_sim+k*2i*pi/par.T,'gs','LineWidth',2,'MarkerSize',7), hold on,grid on
end
title('lam')


subplot(2,2,3:4)
semilogy(abs(d),'gs','LineWidth',2,'MarkerSize',7),grid on
title('abs(mu)')

drawnow


do_fixpoint_iteration=false;
if do_fixpoint_iteration
    for k=1:40%fixpoint iteration s0 will (should) converge to the fixpoint
        %convergence speed depends on the abs. val. of largest mu, which is not computed
        figure(44),clf
        par.doplot=true;
        s0flat =real( v0flat - V * ( (V'*V) \ (V' * (v0flat-s0flat)) .* ((d) ./ (d - 1.0))) ) ;% imaginary part is created only by numerical errors
        v0flat=LinMapsqueeze(s0flat,par);
        errorestimation=norm(v0flat-s0flat)/numel(s0flat);
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

    %
    % %long simulation from the fix point a0
    % %it stay closed for "long" time even in case of unstable fixpoint
    % subplot(2,2,4)
    % par.T = par.T*20; % Update the time period of a longer simulation Pariod of Mapping
    % par.doplot=true;
    % a0_long=LinMap(a0,par); %start a long simulation form the detected fixed point
    %
end


%% -------- Stability Chart with MDBM-------------
warning off
%Very slow solution! ~120s
addpath('C:\Users\Bachrathy\Documents\git\MDBM-Matlab\code_folder')
ax=[];
mdbm_options=mdbmset('timelimit',Inf);

ax(1).val=linspace(-5000,100,10);%Kpphi
ax(2).val=linspace(-150,0,10);%Kdphi
mdbm_sol=mdbm(ax,'fun_Stab_chart',2,mdbm_options,par);

figure(3)
plot_mdbm(mdbm_sol,'k');
view(2)
