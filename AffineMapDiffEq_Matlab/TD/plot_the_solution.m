function plot_the_solution(t,y)
 subplot(4,1,1)
 plot(t,y(:,1),'-');xlabel('Time [sec]');ylabel('Lat.speed [m/s]');grid on;hold on
 subplot(4,1,2)
 plot(t,y(:,2),'-');xlabel('Time [sec]');ylabel('Tilting rate [rad/s]');grid on;hold on
 subplot(4,1,3)
 plot(t,y(:,3),'-');xlabel('Time [sec]');ylabel('Yaw rate [rad/s]');grid on;hold on
 subplot(4,1,4)
 plot(t,y(:,4),'-');xlabel('Time [sec]');ylabel('Tilting angle [rad]');grid on;hold on
