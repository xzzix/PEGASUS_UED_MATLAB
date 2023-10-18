clear all
close all
ek=7.75;
% xis=linspace(0.8,0.85,5);
xis=0.7;
dbfps=linspace(-0.02,0.02,10);
Gg0=[2.2197    3.6599    1.4900]; % guess gradient green quads.
s0=[0.0118    0.0218    0.0380]; %initial guess for PMQ peak field positions
[Xp, Yp]=meshgrid(xis,dbfps);


I1s=zeros(size(Xp));
I2s=zeros(size(Xp));
I3s=zeros(size(Xp));
s12s=zeros(size(Xp));
s23s=zeros(size(Xp));
R12s=zeros(size(Xp));
R34s=zeros(size(Xp));
R11s=zeros(size(Xp));
R33s=zeros(size(Xp));

[S T]=size(Xp)
for u=1:S
    for v=1:T
    
    dbfp=Yp(u,v); %beam divergence approximately
    xi=Xp(u,v); %for gaussian dist
    Distance_sample_to_detector=1.5736;
    Diffraction_optics_optimization
    
    I1s(u,v)=Gg(1)/0.45;
    I2s(u,v)=Gg(2)/0.45;
    I3s(u,v)=Gg(3)/0.45;
    s12s(u,v)=s12;
    s23s(u,v)=s23;
    R11s(u,v)=Cpointsx(end);
    R33s(u,v)=Cpointsy(end);
    R12s(u,v)=Spointsx(end);
    R34s(u,v)=Spointsy(end);
    end
end


%%
figure(11)

ax1=subplot(1,3,1);
plot(Yp*1e+3,I1s,'-o')
xlabel('\delta (mm)')
ylabel('I_1 (A)')
title('I_1 (A)')
set(gca,'xscale','linear','FontSize',14)
box on
ax = gca;
ax.LineWidth = 2;

ax2=subplot(1,3,2);
plot(Yp*1e+3,I2s,'-o')
xlabel('\delta (mm)')
ylabel('I_2 (A)')
set(gca,'xscale','linear','FontSize',14)
box on
ax = gca;
ax.LineWidth = 2;
hold on

ax3=subplot(1,3,3);
plot(Yp*1e+3,I3s,'-o')
xlabel('\delta (mm)')
ylabel('I_3 (A)')
set(gca,'xscale','linear','FontSize',14)
box on
ax = gca;
ax.LineWidth = 2;
hold on


figure(2)

ax1=subplot(1,2,1);
plot(Yp*1e+3,s12s*1e+3,'-o')
xlabel('\delta (mm)')
ylabel('s_{12} (mm)')
set(gca,'xscale','linear','FontSize',14)
box on
ax = gca;
ax.LineWidth = 2;
hold on

ax2=subplot(1,2,2);
plot(Yp*1e+3,s23s*1e+3,'-o')
xlabel('\delta (mm)')
ylabel('s_{23} (mm)')
set(gca,'xscale','linear','FontSize',14)
box on
ax = gca;
ax.LineWidth = 2;
hold on


figure(3)
ax1=subplot(2,2,1);
plot(Yp*1e+3,R11s,'-o')
xlabel('\delta (mm)')
ylabel('R_{11}')
set(gca,'xscale','linear','FontSize',14)
box on
ax = gca;
ax.LineWidth = 2;
hold on

ax2=subplot(2,2,2);
plot(Yp*1e+3,R12s,'-o')
xlabel('\delta (mm)')
ylabel('R_{12} (m)')
set(gca,'xscale','linear','FontSize',14)
box on
ax = gca;
ax.LineWidth = 2;
hold on

ax3=subplot(2,2,3);
plot(Yp*1e+3,R33s,'-o')
xlabel('\delta (mm)')
ylabel('R_{33}')
set(gca,'xscale','linear','FontSize',14)
box on
ax = gca;
ax.LineWidth = 2;
hold on

ax4=subplot(2,2,4);
plot(Yp*1e+3,R34s,'-o')
xlabel('\delta (mm)')
ylabel('R_{34} (m)')
set(gca,'xscale','linear','FontSize',14)
box on
ax = gca;
ax.LineWidth = 2;
hold on

% %%
% figure(1)
% 
% ax1=subplot(1,3,1);
% contourf(Xp,Yp*1e+3,I1s)
% xlabel('\xi (m)')
% ylabel('\delta (mm)')
% title('I_1 (A)')
% set(gca,'xscale','linear','FontSize',14)
% box on
% colorbar
% ax = gca;
% ax.LineWidth = 2;
% 
% ax2=subplot(1,3,2);
% contourf(Xp,Yp*1e+3,I2s)
% xlabel('\xi (m)')
% ylabel('\delta (mm)')
% title('I_2 (A)')
% set(gca,'xscale','linear','FontSize',14)
% box on
% colorbar
% ax = gca;
% ax.LineWidth = 2;
% hold on
% 
% ax3=subplot(1,3,3);
% contourf(Xp,Yp*1e+3,I3s)
% xlabel('\xi (m)')
% ylabel('\delta (mm)')
% title('I_3 (A)')
% set(gca,'xscale','linear','FontSize',14)
% box on
% colorbar
% ax = gca;
% ax.LineWidth = 2;
% hold on
% 
% figure(2)
% 
% ax1=subplot(1,2,1);
% contourf(Xp,Yp*1e+3,s12s*1e+3)
% xlabel('\xi (m)')
% ylabel('\delta (mm)')
% title('s_{12} (mm)')
% set(gca,'xscale','linear','FontSize',14)
% box on
% colorbar
% ax = gca;
% ax.LineWidth = 2;
% hold on
% 
% ax2=subplot(1,2,2);
% contourf(Xp,Yp*1e+3,s23s*1e+3)
% xlabel('\xi (m)')
% ylabel('\delta (mm)')
% title('s_{23} (mm)')
% set(gca,'xscale','linear','FontSize',14)
% box on
% colorbar
% ax = gca;
% ax.LineWidth = 2;
% hold on
% 
% 
% figure(3)
% ax1=subplot(2,2,1);
% contourf(Xp,Yp*1e+3,R11s)
% xlabel('\xi (m)')
% ylabel('\delta (mm)')
% title('R_{11}')
% set(gca,'xscale','linear','FontSize',14)
% box on
% colorbar
% ax = gca;
% ax.LineWidth = 2;
% hold on
% 
% ax2=subplot(2,2,2);
% contourf(Xp,Yp*1e+3,R12s)
% xlabel('\xi (m)')
% ylabel('\delta (mm)')
% title('R_{12} (m)')
% set(gca,'xscale','linear','FontSize',14)
% box on
% colorbar
% ax = gca;
% ax.LineWidth = 2;
% hold on
% 
% ax3=subplot(2,2,3);
% contourf(Xp,Yp*1e+3,R33s)
% xlabel('\xi (m)')
% ylabel('\delta (mm)')
% title('R_{33}')
% set(gca,'xscale','linear','FontSize',14)
% box on
% colorbar
% ax = gca;
% ax.LineWidth = 2;
% hold on
% 
% ax4=subplot(2,2,4);
% contourf(Xp,Yp*1e+3,R34s)
% xlabel('\xi (m)')
% ylabel('\delta (mm)')
% title('R_{34} (m)')
% colorbar
% set(gca,'xscale','linear','FontSize',14)
% box on
% ax = gca;
% ax.LineWidth = 2;
% hold on