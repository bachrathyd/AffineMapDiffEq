%% Simulation for balancing a motorcycle 
% Control with steering, hierarchical control
% no longitudinal speed, v=0
% phi' and delta' as a pseudovelocities
% multiple delays (for steering and for lean)

clear all;
% close all;
clc

global Kpphi Kdphi Kpdelta Kddelta
tau1 = 0.01; % for steering (delta)
tau2 = 0.2*1e-2; % for lean (phi)
lags = [tau1,tau2];
t_end = 1.5;

color = 'k';
options = ddeset('InitialY',[0; 0; 0.01; 0],'MaxStep',0.001,'RelTol',1e-5,'AbsTol',1e-7); % ,'InitialStep',0.001
sol = dde23(@motor_rhs_control_linear_hierarchical_multipledelays, lags, @(t) zeros(4,1) ,[0, t_end], options); 
 
do_plot=0
if do_plot
figure %(2); hold on; %clf;
%set(gcf,'Units','centimeters','Position',[2 6 10 5.5])
line_x = [0 t_end];
line_y = [0 0];
% phi(t)
subplot(5,1,1)
plot(sol.x,sol.y(1,:),color,'Linewidth', 1)
ax = gca;
ax.FontName = 'TimesNewRoman';
ax.FontSize = 8;
xlabel('$t [\mathrm{s}]$','Interpreter','latex','FontSize',8,'FontName','TimesNewRoman');
ylabel('$\varphi [\mathrm{rad}]$','Interpreter','latex','FontSize',8,'FontName','TimesNewRoman');
grid off
title(['$K_{p \varphi} = \,$',num2str(round(Kpphi,4)),' \,, $K_{p \delta} = \,$',num2str(round(Kpdelta,4)), ' {\rm Nm} \,, $K_{d \varphi} = \, $',num2str(round(Kdphi,4)),'$\, {\rm s} \,, K_{d \delta} = \, $',num2str(round(Kddelta,4)), '$ \, {\rm Nms} \,, \tau_{\delta} = \,$', num2str(round(tau1,5)), '$ \, {\rm s} \,, \tau_{\varphi} = $', num2str(round(tau2,5)), '$ \, {\rm s}$'],'Interpreter','latex')
hold on
line(line_x,line_y,'Color','black','LineStyle','-','Linewidth', 0.2)
% ylim([-10*1e-4 6*1e-4])
% delta(t)
subplot(5,1,2)
plot(sol.x,sol.y(2,:),color,'Linewidth', 1)
ax = gca;
ax.FontName = 'TimesNewRoman';
ax.FontSize = 8;
xlabel('$t [\mathrm{s}]$','Interpreter','latex','FontSize',8,'FontName','TimesNewRoman');
ylabel('$\delta [\mathrm{rad}]$','Interpreter','latex','FontSize',8,'FontName','TimesNewRoman');
grid off
hold on
line(line_x,line_y,'Color','black','LineStyle','-','Linewidth', 0.2)
% ylim([-0.01 0.08])
% phi'(t)
subplot(5,1,3)
plot(sol.x,sol.y(3,:),color,'Linewidth', 1)
ax = gca;
ax.FontName = 'TimesNewRoman';
ax.FontSize = 8;
xlabel('$t [\mathrm{s}]$','Interpreter','latex','FontSize',8,'FontName','TimesNewRoman');
ylabel('$\dot{\varphi} [\mathrm{rad/s}]$','Interpreter','latex','FontSize',8,'FontName','TimesNewRoman');
grid off
hold on
line(line_x,line_y,'Color','black','LineStyle','-','Linewidth', 0.2)
% ylim([-0.025 0.025])
% delta'(t)
subplot(5,1,4)
plot(sol.x,sol.y(4,:),color,'Linewidth', 1)
ax = gca;
ax.FontName = 'TimesNewRoman';
ax.FontSize = 8;
xlabel('$t [\mathrm{s}]$','Interpreter','latex','FontSize',8,'FontName','TimesNewRoman');
ylabel('$\dot{\delta} [\mathrm{rad/s}]$','Interpreter','latex','FontSize',8,'FontName','TimesNewRoman');
grid off
hold on
line(line_x,line_y,'Color','black','LineStyle','-','Linewidth', 0.2)
% ylim([-1 2])
% M(t) =
% -Kpphi*Kpdelta*phi(t-tau2)-Kpdelta*delta(t-tau1)-Kdphi*Kpdelta*phi'(t-tau2)-Kddelta*delta'(t-tau1)
subplot(5,1,5)
plot(sol.x,-Kpphi*Kpdelta*sol.y(1,:)-Kpdelta*sol.y(2,:)-Kdphi*Kpdelta*sol.y(3,:)-Kddelta*sol.y(4,:),color,'Linewidth', 1)
ax = gca;
ax.FontName = 'TimesNewRoman';
ax.FontSize = 8;
xlabel('$t [\mathrm{s}]$','Interpreter','latex','FontSize',8,'FontName','TimesNewRoman');
ylabel('$M [\mathrm{Nm}]$','Interpreter','latex','FontSize',8,'FontName','TimesNewRoman');
grid off
hold on
line(line_x,line_y,'Color','black','LineStyle','-','Linewidth', 0.2)
% ylim([-3 4])
end
