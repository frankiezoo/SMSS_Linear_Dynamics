% function SMSS_lin_dynamics(dr,dm)
%% set up geometry of system
mu_0=4*pi*10^-7;    % permeability of free space
mu=0.8815;          % magnetic field strength of .75 in D N42 spherical NdFeB magnet

r_FC=[0 0 0.01]';   % field cooled position of magnet, separation distance of 1 cm above superconductor
mu_FC=mu*[0 0 1]';     % orientation of field cooled magnet, pointing upward
r_e=[0 0 0.01]'; % position of magnet, not at equilibrium
m_e=mu*[0 0 1]';    % orientation of magnet, pointing upward

a_s=[0 0 1]';       % superconductor normal
r_fe=(eye(3)-2*(a_s*a_s'))*r_FC;  % position of frozen image
r_me=(eye(3)-2*(a_s*a_s'))*r_e;   % position of mobile image
m_fe=mu_FC;          % orientation of frozen image
m_me=(eye(3)-2*(a_s*a_s'))*m_e;  % orientation of mobile image
rho_fe=r_e-r_fe;    % equilibrium distance from frozen image to magnet
rho_me=r_e-r_me;    % equilibrium distance from mobile image to magnet
m=0.0272;           % mass of spherical magnet
M=m*eye(3);         % mass matrix of spherical magnet
I=2/5*M*(0.75*2.54/100)^2;  % inertia of spherical magnet
%% find analytical solution with linearized plant
dt=1e-3;        % timestep for output dynamics
tf=1;          % total simulated time
t=0:dt:tf;      % simulation time vector
dr=[0 0 0 ]'; % initial translational displacement
theta=0;     % initial angular displacement
Rx=[1 0 0;
    0 cos(theta) sin(theta);
    0 -sin(theta) cos(theta)];
Ry=[cos(theta) 0 -sin(theta);
    0 1 0;
    sin(theta) 0 cos(theta)];
dm=m_e-Ry*m_e;    % initial orientation displacement
dq=[dm; 0];
dv=[0 0 0]';    % initial translational velocity
domega=[0 0 0]';% initial angular velocity
z0=[dr; dv; dq; domega];    % initial state of system

[A, b]=SMSS_state_space(r_e,m_e,a_s);

z=zeros(13,length(t));
Ftau=zeros(6,length(t));
for i=1:length(t)
    z(:,i)=expm(A*t(i))*z0;
    dz=A*z(:,i)*dt;
    Ftau(:,i)=[M*dz(4:6);I*dz(11:13)];
end

dx=z(1,:);
dy=z(2,:);
dz=z(3,:);
dvx=z(4,:);
dvy=z(5,:);
dvz=z(6,:);
dmx=z(7,:);
dmy=z(8,:);
dmz=z(9,:);
dms=z(10,:);
dox=z(11,:);
doy=z(12,:);
doz=z(13,:);

%% find simulation solution
r_mag_0=r_e+dr;
mu_mag_0=m_e+dm;
sim_time=tf;
damping=0;

tic
sim('SMSS_nonlinear')
toc

r_mag_sx = spline(r_mag.time,r_mag.data(:,1),t);
r_mag_sy = spline(r_mag.time,r_mag.data(:,2),t);
r_mag_sz = spline(r_mag.time,r_mag.data(:,3),t);
dx_s = r_mag_sx'-ones(size(t'))*r_FC(1);
dy_s = r_mag_sy'-ones(size(t'))*r_FC(2);
dz_s = r_mag_sz'-ones(size(t'))*r_FC(3);
dvx_s = spline(vel_N.time,vel_N.data(:,1),t);
dvy_s = spline(vel_N.time,vel_N.data(:,2),t);
dvz_s = spline(vel_N.time,vel_N.data(:,3),t);
err_x = sqrt((dx_s'-dx).^2)/(max(dx)-min(dx));
err_y = sqrt((dy_s'-dy).^2)/(max(dy)-min(dy));
err_z = sqrt((dz_s'-dz).^2)/(max(dz)-min(dz));
err_vx = sqrt((dvx_s-dvx).^2)/(max(dvx)-min(dvx));
err_vy = sqrt((dvy_s-dvy).^2)/(max(dvy)-min(dvy));
err_vz = sqrt((dvz_s-dvz).^2)/(max(dvz)-min(dvz));

mu_mag_sx = spline(mu_mag.time,mu_mag.data(:,1),t);
mu_mag_sy = spline(mu_mag.time,mu_mag.data(:,2),t);
mu_mag_sz = spline(mu_mag.time,mu_mag.data(:,3),t);
dmx_s = mu_mag_sx'-ones(size(t'))*mu_FC(1);
dmy_s = mu_mag_sy'-ones(size(t'))*mu_FC(2);
dmz_s = mu_mag_sz'-ones(size(t'))*mu_FC(3);
dox_s = spline(omega_b.time,omega_b.data(:,1),t);
doy_s = spline(omega_b.time,omega_b.data(:,2),t);
doz_s = spline(omega_b.time,omega_b.data(:,3),t);
err_phix = sqrt((dmx_s'-dmx).^2)/(max(dmx)-min(dmx));
err_phiy = sqrt((dmy_s'-dmy).^2)/(max(dmy)-min(dmy));
err_phiz = sqrt((dmz_s'-dmz).^2)/(max(dmz)-min(dmz));
err_omegax = sqrt((dox_s-dox).^2)/(max(dox)-min(dox));
err_omegay = sqrt((doy_s-doy).^2)/(max(doy)-min(doy));
err_omegaz = sqrt((doz_s-doz).^2)/(max(doz)-min(doz));
%% graphics
close all

figure
subplot(2,2,1)
plot(t,dx,'k',t,dy,'--k',t,dz,':k')
xlabel('time [sec]')
ylabel('position displacement [m]')
title('linearized dr')
subplot(2,2,2)
plot(t,dx_s,'k',t,dy_s,'--k',t,dz_s,':k')
xlabel('time [sec]')
ylabel('position displacement [m]')
title('simulation dr')
subplot(2,2,3)
plot(t,dx_s'-dx,'k',t,dy_s'-dy,'--k',t,dz_s'-dz,':k')
xlabel('time [sec]')
ylabel('error position displacement [m]')
title('error of dr')
subplot(2,2,4)
plot(t,err_x,'k',t,err_y,'--k',t,err_z,':k')
xlabel('time [sec]')
ylabel('error position displacement [m]')
title('NMRSD error dr')
legend('x','y','z')

figure
subplot(2,2,1)
plot(t,dvx,'k',t,dvy,'--k',t,dvz,':k')
xlabel('time [sec]')
ylabel('translational velocity [m/s]')
title('linearized dv')
subplot(2,2,2)
plot(t,dvx_s,'k',t,dvy_s,'--k',t,dvz_s,':k')
xlabel('time [sec]')
ylabel('translational velocity [m/s]')
title('simulation dv')
subplot(2,2,3)
plot(t,dvx_s-dvx,'k',t,dvy_s-dvy,'--k',t,dvz_s-dvz,':k')
xlabel('time [sec]')
ylabel('error translational velocity [m/s]')
title('error dv')
subplot(2,2,4)
plot(t,err_vx,'k',t,err_vy,'--k',t,err_vz,':k')
xlabel('time [sec]')
ylabel('error trans velocity [m/s]')
title('NMRSD error dv')
legend('x','y','z')

figure
subplot(2,2,1)
plot(t,dmx,'k',t,dmy,'--k',t,dmz,':k')
xlabel('time [sec]')
ylabel('angular displacement [rad]')
title('linearized dm')
subplot(2,2,2)
plot(t,dmx_s,'k',t,dmy_s,'--k',t,dmz_s,':k')
xlabel('time [sec]')
ylabel('angular displacement [rad]')
title('simulation dm')
subplot(2,2,3)
plot(t,dmx_s-dmx','k',t,dmy_s-dmy','--k',t,dmz_s-dmz',':k')
xlabel('time [sec]')
ylabel('angular displacement [rad]')
title('error dm')
subplot(2,2,4)
plot(t,err_phix,'k',t,err_phiy,'--k',t,err_phiz,':k')
xlabel('time [sec]')
ylabel('angular displacement [rad]')
title('NMRSD error dm')
legend('x','y','z')

figure
subplot(2,2,1)
plot(t,dox,'k',t,doy,'--k',t,doz,':k')
xlabel('time [sec]')
ylabel('angular velocity [rad/s]')
title('linearized \omega')
subplot(2,2,2)
plot(t,dox_s,'k',t,doy_s,'--k',t,doz_s,':k')
xlabel('time [sec]')
ylabel('angular velocity [rad/s]')
title('simulation \omega')
subplot(2,2,3)
plot(t,dox_s-dox,'k',t,doy_s-doy,'--k',t,doz_s-doz,':k')
xlabel('time [sec]')
ylabel('angular velocity [rad/s]')
title('error\omega')
subplot(2,2,4)
plot(t,err_omegax,'k',t,err_omegay,'--k',t,err_omegaz,':k')
xlabel('time [sec]')
ylabel('angular velocity [rad/s]')
title('NMRSD error \omega')
legend('x','y','z')