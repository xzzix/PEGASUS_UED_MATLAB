
%% Matrix Simulation Parameters
%physical constants
m0 = 9.10938215.*10^(-31);
e0 = 1.60217646.*10.^(-19);
c0 = 2.99792458.*10.^8;
eps0=8.85418782e-12; %m^-3 kg^-1 s^4 A^2
%beam parameters
% ek= 8.3671; %kinetic energy
dG=1e-5; %rms energy spread
gamma = 1+ek/0.511;
beta = (1-1./(gamma.^(2))).^(1/2);
bbrho = gamma.*beta.*m0.*c0./e0;
beamdiv=0.1e-3; %beam divergence approximately
Spotsize=100.0e-6; %for gaussian dist
Qbeam=-1.6e-14*0 ;%30e-15;   %beam chargeR   microscope has 1pc potentially
ztime=10e-12;  %temporal length
Lbeam=c0.*beta.*ztime;
lambda=Qbeam./Lbeam; %lin charge density
perveance=-e0*lambda./(gamma.^3.*m0.*(c0*beta).^2*2*pi*eps0);
geoemitx=beamdiv*Spotsize;
geoemity=beamdiv*Spotsize;
CSalphax=0;
CSalphay=0;
CSbetax=Spotsize^2/geoemitx;
CSbetay=Spotsize^2/geoemity;
CSgammax=(1+CSalphax^2)/CSbetax;
CSgammay=(1+CSalphay^2)/CSbetay;
sigma0=eye(4);
sigma0(1,1)=geoemitx*CSbetax;
sigma0(1,2)=-geoemitx*CSalphax;
sigma0(2,1)=-geoemitx*CSalphax;
sigma0(2,2)=geoemitx*CSgammax;
sigma0(3,3)=geoemity*CSbetay;
sigma0(3,4)=-geoemity*CSalphay;
sigma0(4,3)=-geoemity*CSalphay;
sigma0(4,4)=geoemity*CSgammay;

%Detector plane and back focal planes.
Distance_sample_to_detector=1.5736;


LBFP=xi+dbfp/2;
bfpshift=dbfp;



%% Green triplet

lg=[0.0768 0.0768 0.0768]; %effective lengths -- green quads.
bg=[135 135 135]; %field curvature parameter -- green quads.
delta = 0.085; %spacing between green quads.
sg=[0.104 0.104+delta 0.104+2*delta]; %positions of green quad centers.
% Gg0=[2.2197    3.6599    1.4900]; % guess gradient green quads.

%% matrix solution -- fringe field calculation + linear space charge
fun1=@(Gg)greentriplet_sc(sg,lg,Gg,bg,gamma,LBFP,0,sigma0(1,1),sigma0(3,3),sigma0(2,2),sigma0(4,4),bfpshift);
options = optimset('Display','iter');
options.MaxIter = 500000 ;
options.MaxFunEvals = 5000000 ;
options.FunctionTolerance = 1e-21;
options.OptimalityTolerance = 1e-21;
options.StepTolerance = 1e-21;
lb = [0.35 0.35 0.35];
ub = 0.47*[8 8 8];
A = [];
b = [];
Aeq = [];
beq = [];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','PlotFcns',@optimplotfval);
nonlcon = [];
[Gg,fval1]=fmincon(fun1,Gg0,A,b,Aeq,beq,lb,ub,nonlcon,options)
Gg0=Gg;
%% Green triplet matrix solution up to BFP-- fringe field calculation + linear space charge
[M_all_green,zsteps1]=greentripletRmatall_sc(sg,lg,Gg,bg,gamma,LBFP,0,sigma0(1,1),sigma0(3,3),sigma0(2,2),sigma0(4,4));
disp('linear transport matrix:');
R0=M_all_green(:,:,end);
%% PMQ matrix optimization-- fringe field calculation + linear space charge
%PMQ lattice design specs
L0=Distance_sample_to_detector-LBFP;%distance from BFP to detector
l0=[6.0e-03 6.0e-03 3.0e-03]; %PMQ effective lengths
G0=[510.4 517.7 416.7]; %PMQ gradients Veronica's measurements
bs=[1.675493408057606   1.675493408057606   1.639343181086262]*1e+3; %PMQ field curvature parameter
% s0=[0.0118    0.0218    0.0380]; %initial guess for PMQ peak field positions

fun1=@(s)quaddiffractionfringe_sc(s,l0,G0,bs,gamma,L0,0,sigma0(1,1),sigma0(3,3),sigma0(2,2),sigma0(4,4),R0,bfpshift);
options = optimset('Display','iter','PlotFcns',@optimplotfval);
options.MaxIter = 500000 ;
options.MaxFunEvals = 5000000 ;
options.FunctionTolerance = 1e-21;
options.OptimalityTolerance = 1e-21;
options.StepTolerance = 1e-21;
[s,fval1]=fminsearch(fun1,s0,options);
s0=s
%% PMQ matrix solution-- fringe field calculation + linear space charge
[M_all,zsteps]=quaddiffractionfringeRmatall_sc(s,l0,G0,bs,gamma,L0,0,sigma0(1,1),sigma0(3,3),sigma0(2,2),sigma0(4,4));
disp('linear transport matrix:');
Mf=M_all(:,:,end)*R0

%% check the quad spacings
disp('check the quad spacings:');
z12=s(2)-s(1)-l0(1)
z23=s(3)-s(2)-l0(1)/2-l0(3)/2

disp('center-to-center spacings:');
s12=s(2)-s(1)
s23=s(3)-s(2)
zsteps=zsteps+LBFP;

%% Principle ray's
%normalized Bfield, to overlay fields with principle rays.
Gz_PMQ=@(z)(zquadrupole(z,G0(1),bs(1),l0(1),s(1))-zquadrupole(z,G0(2),bs(2),l0(2),s(2))+zquadrupole(z,G0(3),bs(3),l0(3),s(3)));
Gz_green=@(z)(zquadrupole(z,Gg(1),bg(1),lg(1),sg(1))-zquadrupole(z,Gg(2),bg(2),lg(2),sg(2))+zquadrupole(z,Gg(3),bg(3),lg(3),sg(3)));

%C and S matrix elements Green triplet transport from sample to BFP
Cpointsx1=zeros(1,length(M_all_green));
Spointsx1=zeros(1,length(M_all_green));
Cpointsy1=zeros(1,length(M_all_green));
Spointsy1=zeros(1,length(M_all_green));

for j=1:length(Spointsx1)
Cpointsx1(j)=M_all_green(1,1,j);
Spointsx1(j)=M_all_green(1,2,j);
Cpointsy1(j)=M_all_green(3,3,j);
Spointsy1(j)=M_all_green(3,4,j);
end

%compose green triplet transport with PMQ, to plot entire transport.
Rmat_green=M_all_green;
Rmat_PMQ=M_all;
Rmat=zeros(4,4,length(Rmat_PMQ));
for j=1:length(Rmat)
    Rmat(:,:,j)=Rmat_PMQ(:,:,j)*Rmat_green(:,:,end);
end

%C and S matrix elements PMQ triplet transport from BFP to detector
Cpointsx=zeros(1,length(Rmat));
Spointsx=zeros(1,length(Rmat));
Cpointsy=zeros(1,length(Rmat));
Spointsy=zeros(1,length(Rmat));

for j=1:length(Spointsx)
Cpointsx(j)=Rmat(1,1,j);
Spointsx(j)=Rmat(1,2,j);
Cpointsy(j)=Rmat(3,3,j);
Spointsy(j)=Rmat(3,4,j);
end

% plot principle ray's
figure(2)
yyaxis left
plot(zsteps1,Cpointsx1,'-','Color',[0,1,0.75],'LineWidth',2)
hold on
plot(zsteps1,Cpointsy1,'-','Color',[0.72,0.27,1.00],'LineWidth',2)
hold on
plot(zsteps,Cpointsx,'-','Color',[0,1,0.75],'LineWidth',2)
hold on
plot(zsteps,Cpointsy,'-','Color',[0.72,0.27,1.00],'LineWidth',2)
hold on
xlabel('z(m)')
ylabel('C')
yyaxis right
plot(zsteps1,Gz_green(zsteps1)/max(Gz_green(zsteps1)),'k-')
hold on
plot(zsteps,Gz_PMQ(zsteps-LBFP)/max(Gz_PMQ(zsteps-LBFP)),'k-')
hold on
xline(0.43815) %box entrance
hold on
xline(1.1176) %box exit
set(gca,'FontSize',15)

figure(3)
yyaxis left
plot(zsteps1,Spointsx1,'-','Color',[0,1,0.75],'LineWidth',2)
hold on
plot(zsteps1,Spointsy1,'-','Color',[0.72,0.27,1.00],'LineWidth',2)
hold on
plot(zsteps,Spointsx,'-','Color',[0,1,0.75],'LineWidth',2)
hold on
plot(zsteps,Spointsy,'-','Color',[0.72,0.27,1.00],'LineWidth',2)
xlabel('z(m)')
ylabel('S(m)')
yyaxis right
hold on
plot(zsteps1,Gz_green(zsteps1)/max(Gz_green(zsteps1)),'k-')
hold on
plot(zsteps,Gz_PMQ(zsteps-LBFP)/max(Gz_PMQ(zsteps-LBFP)),'k-')
hold on
xline(0.43815) %box entrance
hold on
xline(1.1176) %box exit
set(gca,'FontSize',15)


