% clear all
% close all
% clc

%% initial time value


% Parameters
Omega = 15; % Rotational speed of the wheel[rad/s]
k = 6e4; % Distributed lateral stiffness [N/m^2]
m = 10; % Mass of the wheel [kg]
RE = 0.2; % (Effective rolling) Radius of the wheel [m]
a = 0.05; % Half-length of the contact patch [m]
g = 9.81; % Gravitational acceleration [m/s^2]

Omega_cr = sqrt(g/3/RE); % Rigid wheel critical speed

% Parameters for the discretization
for n = 160:-5:30%10:1:100 % Number of discretization elements in the contact patch
disp(n)
    dx = 2*a/n;
x = linspace(-a,a,n+1);
% timestep size
dt = dx/(RE*Omega)/2;



% IC - Initial Conditions
%       y = [ sigma1; sigma2; sigma3; theta; q(:) ];
%           sigma1 - lateral velocity of the contact patch center point
%           sigma2 - tilting rate / rotational speed along the x-axis
%           sigma3 - yaw rate / rotational speed along z-axis
%           theta - tilting angle
%           q(:) - discretized lateral tire deformation q(x,t)
y0 = [0; 0.1; 0; 0; zeros(n+1,1)];
%y0 = [0; 2.5; 0; 0; zeros(n+1,1)];

% y0 =ones(5+n,1);

par = [];
par.Omega=Omega;
par.k=k;
par.m=m;
par.RE=RE;
par.g=g;
par.x=x;
par.dx=dx;
par.n=n;


% initial time value
t0=0;
% final time value
tf = 2.0;

% [t,y] = ode23(@equ,[t0,tf],y0,odeset('RelTol',1e-7,'AbsTol',1e-7),par);
[t,y] = ode23(@equ_lin,[t0,tf],y0,odeset('RelTol',1e-7,'AbsTol',1e-7),par);

% Animate the deformation
if false
    for i=1:length(t)
        figure(3)
        q=y(i,5:n+5);
        plot(x(:),q,'.-')
        grid on;
        ylim([-0.002,0.002]);
        drawnow;
    end
end

figure(1);%clf
set(gcf,'position',[200,100,500,800]);
plot_the_solution(t,y)

%% ------------------ Affine mapping -----------------
% Parameters for the diff.eq.


% paramteres for the spactrum calculation
par.T = 0.001; % Pariod of Mapping
par.Ndim=n+5; % Dimension of the system

par.diffun=@equ; % The function handle for the diff.eq. par.diffun=@(t,y,Z,par) equ_lin(t,y,par); % The function handle for the diff.eq.
par.options = odeset('RelTol',1e-2,'AbsTol',1e-2); %options of the dde solver
par.lags=[];%lags used in the dde solver
%defining the sampling of the delay (states)
taumax=0; %max delay - it must be larger or equal to the maximal delyay (important if it is not constant)
Nstep=1;

% par.diffun=@(t,y,Z,par) equ(t,y,par); % The function handle for the diff.eq.
% % par.diffun=@(t,y,Z,par) equ_lin(t,y,par); % The function handle for the diff.eq.
% par.options = ddeset('RelTol',1e-4,'AbsTol',1e-4)%,'MaxStep',par.T/100); %options of the dde solver
% par.lags=[tau];%lags used in the dde solver
% %defining the sampling of the delay (states)
% tau = 0.001; %  %it is not delayed!!! ->dummy parameter
% %defining the sampling of the delay (states)
% taumax=tau; %max delay - it must be larger or equal to the maximal delyay (important if it is not constant)
% Nstep=1000


par.StateSmapled=linspace(-taumax,0,Nstep);%no-need for uniform sampling


par.doplot=false;%true;%false % plot the mapped results - for debuging only

% testing the mapping- for the matrix form
%s0start=zeros(par.Ndim,Nstep); %necessary initalization for the fixpoint and the eigen value calculations
% s0start=ones(par.Ndim,Nstep)*0.001;
% s0start=rand(par.Ndim,Nstep);
s0start=repmat(y0,1,Nstep);
% tic
v0start=LinMap(s0start,par);
% toc

figure(1)
%plotting the test initial condition-----------
plot_the_solution(par.StateSmapled,s0start')
%plotting the test mapped result-----------
plot_the_solution(par.StateSmapled+par.T,v0start')


% testing the mapping- for the vectorized form
s0flat=s0start(:);
v0flat=LinMapsqueeze(s0flat,par);


%plotting the test initial condition-----------
plot_the_solution(par.StateSmapled,reshape(s0flat,par.Ndim,[])')
%plotting the test mapped result-----------
plot_the_solution(par.StateSmapled+par.T,reshape(v0flat,par.Ndim,[])')

% <<<<<<<<<<<<<<<<<<<<< testing the spectrum & the fixed point >>>>>>>>>>>>>>>>

Norm_of_perturbation=10^-2;


eigsN=10;%20;
tic
opts=[];
opts.p=30;%30;
% opts.tol = (1e-20)/n^2;%15;
opts.tol = (1e-20);%15;
opts.disp=0;
opts.issym=false;
% opts.v0=ones(par.Ndim*Nstep,1)
% opts.v0=zeros(par.Ndim*Nstep,1)
% opts.v0(2)=1
opts.maxit=30;
opts.v0=s0flat;
opts.fail='keep';
par.s0flat=s0flat;
% Norm_of_perturbation=10e-5
[V,D] = eigs(@(s) (LinMapsqueeze(s .* Norm_of_perturbation,par)) ./ Norm_of_perturbation,numel(s0start),eigsN,'largestabs',opts);%LinearMap
% [V,D] = eigs(@(s) (LinMapsqueeze(s*Norm_of_perturbation+s0flat,par)-v0flat)/Norm_of_perturbation,numel(s0start),eigsN,'largestabs','Tolerance',1e-9);%Affine
toc
% plotting the specttral properties
mus=diag(D);



figure(11)
%clf
% subplot(2,2,1)
% plot(mus,'s','LineWidth',2,'MarkerSize',7), hold on
% fi=linspace(0,2*pi,1000);
% plot(sin(fi),cos(fi)),grid on
% title('mu')
% xlabel('real(mu)')
% ylabel('imag(mu)')

subplot(2,2,2)
lam_sim=log(mus)/par.T;

plot(lam_sim,'+','LineWidth',2,'MarkerSize',7), hold on,grid on
% for k=-1:1
% plot(lam_sim+k*2i*pi/par.T,'gs','LineWidth',2,'MarkerSize',7), hold on,grid on
% end
title('lam')
xlabel('real(lam)')
ylabel('imag(lam)')

subplot(2,2,3)
plot(n+ 0*real(lam_sim),real(lam_sim)*n^2,'+','LineWidth',2,'MarkerSize',7),hold on, grid on
% plot(real(lam_sim)*n^2,'o','LineWidth',2,'MarkerSize',7),hold on, grid on
title('real(lam)*n^2')
xlabel('n')
ylabel('real(lam)*n2')

subplot(2,2,4)
plot(n+ 0*real(lam_sim),real(lam_sim),'+','LineWidth',2,'MarkerSize',7),hold on, grid on
% plot(real(lam_sim)*n^2,'o','LineWidth',2,'MarkerSize',7),hold on, grid on
title('real(lam)')
xlabel('n')
ylabel('real(lam)')
drawnow

subplot(2,2,1)
 plot(real(lam_sim),'+','LineWidth',2,'MarkerSize',7),hold on, grid on
 title('real(lam_i)')
xlabel('i')
ylabel('real(lam)_i')
drawnow
end

% %plot the mode shapes
% figure(2),clf
% for kmode=1:eigsN
%         figure(2),clf
%     ai_modeshape=reshape(V(:,kmode),par.Ndim,[]);
%     norm((D(kmode,kmode)*ai_modeshape)-LinMap(ai_modeshape,par))/length(ai_modeshape(:))
%     plot_the_solution(par.StateSmapled,imag(ai_modeshape.'))
%     plot_the_solution(par.StateSmapled,real(ai_modeshape.'))
%     plot_the_solution(par.StateSmapled+par.T,imag((D(kmode,kmode)*ai_modeshape).'))
%     plot_the_solution(par.StateSmapled+par.T,real((D(kmode,kmode)*ai_modeshape).'))
%     plot_the_solution(par.StateSmapled+par.T,imag(LinMap(ai_modeshape,par).'))
%     plot_the_solution(par.StateSmapled+par.T,real(LinMap(ai_modeshape,par).'))
%     title({"Modeshape i, \mu_i",num2str(kmode),num2str(mus(kmode))})
%     drawnow
%     pause(1)
% end


% % do_fixpoint_iteration=false;
% % if do_fixpoint_iteration
% %     for k=1:40%fixpoint iteration s0 will (should) converge to the fixpoint
% %         %convergence speed depends on the abs. val. of largest mu, which is not computed
% %         figure(44),clf
% %         par.doplot=true;
% %         s0flat =real( v0flat - V * ( (V'*V) \ (V' * (v0flat-s0flat)) .* ((d) ./ (d - 1.0))) ) ;% imaginary part is created only by numerical errors
% %         v0flat=LinMapsqueeze(s0flat,par);
% %         errorestimation=norm(v0flat-s0flat)/numel(s0flat);
% %         if (errorestimation)<1e-10
% %             a0flat=s0flat;
% %             break
% %         end
% %         drawnow
% %     end
% %     a0flat=s0flat;
% %     a0=reshape(a0flat,par.Ndim,[]);
% %
% %     figure(2345)
% %     plot(par.StateSmapled,a0)
% %
% %     %
% %     % %long simulation from the fix point a0
% %     % %it stay closed for "long" time even in case of unstable fixpoint
% %     % subplot(2,2,4)
% %     % par.T = par.T*20; % Update the time period of a longer simulation Pariod of Mapping
% %     % par.doplot=true;
% %     % a0_long=LinMap(a0,par); %start a long simulation form the detected fixed point
% %     %
% % end

% % % % % % % % 
% % % % % % % % Omega = 15; % Rotational speed of the wheel[rad/s]
% % % % % % % % k = 6e4; % Distributed lateral stiffness [N/m^2]
% % % % % % % % m = 10; % Mass of the wheel [kg]
% % % % % % % % RE = 0.2; % (Effective rolling) Radius of the wheel [m]
% % % % % % % % a = 0.05; % Half-length of the contact patch [m]
% % % % % % % % g = 9.81; % Gravitational acceleration [m/s^2]
% % % % % % % % %% -------- Stability Chart with MDBM-------------
% % % % % % % % %Ezt lefut, csak épp nincs stabil terület
% % % % % % % % warning off
% % % % % % % % %Very slow solution! ~120s
% % % % % % % % addpath('C:\Users\Bacharthy\Documents\GitHub\MBDM\code_folder')
% % % % % % % % ax=[];
% % % % % % % % ax(1).val=linspace(5,15,5);%Omega
% % % % % % % % ax(2).val=linspace(0.01,0.2,5);%a
% % % % % % % % mdbm_options=mdbmset('timelimit',Inf);
% % % % % % % % mdbm_sol=mdbm(ax,'fun_Stab_chart',2,mdbm_options,par);
% % % % % % % % 
% % % % % % % % figure(33)
% % % % % % % % plot_mdbm(mdbm_sol,'k');
% % % % % % % % view(2)
