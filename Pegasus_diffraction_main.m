clear all
close all

path_to_gpt='"C:\Program Files\General Particle Tracer\bin\"';
path_to_in_out_file=' C:\Users\Paul\Desktop\Pegasus_diffraction_simulation_7_18_2023\' ;
license=' GPTlicense=1476385047';


%% Injector to sample GPT simulation
%injector parameteres
E0 = 0.4;
nps = 20000;
sigmaX = 150e-6;
sigmaY = 150e-6;
sigmaT = 30e-15;
qTotal = -100e-15;
gunLoopmV = 43;
gunPhaseDeg = 25;
solFac1 = 1.15;
S2val=0;
XPhaseDeg=21.11;
XFieldScaleFac=1.01*0;
linacPhaseDeg=218;
linacFieldScaleFac=18.75e+6;

infile = 'injection.in';
gdffile = 'injection_solution.gdf';
statsfile1 = 'injection_stats.gdf';
callGPT = strcat(path_to_gpt,'gpt -v -o ',path_to_in_out_file,gdffile,path_to_in_out_file,infile,...
    ' E0=',num2str(E0),...
    ' nps=',num2str(nps),' sigmaX=',num2str(sigmaX),' sigmaY=',num2str(sigmaY),...
    ' sigmaT=',num2str(sigmaT),' qTotal=',num2str(qTotal),' gunLoopmV=',num2str(gunLoopmV),...
    ' gunPhaseDeg=',num2str(gunPhaseDeg),' solFac1=',num2str(solFac1),' S2val=',num2str(S2val), ...
    ' XPhaseDeg=',num2str(XPhaseDeg), ' XFieldScaleFac=',num2str(XFieldScaleFac),...
    ' linacPhaseDeg=',num2str(linacPhaseDeg),' linacFieldScaleFac=',num2str(linacFieldScaleFac), ...
    license)
system(callGPT,'-echo') ;
callGPT2=strcat(path_to_gpt,'gdfa -o',path_to_in_out_file,statsfile1,path_to_in_out_file,gdffile,' time ',' stdBy ',' stdBx ',' stdx ',' stdy',' CSgammax',' CSgammay',' CSalphax',' CSalphay', ' CSbetax',' CSbetay',' nemixrms',' nemiyrms',' nemirrms',' avgz',' stdt',' stdz',' avgG',' stdG');
system(callGPT2,'-echo') ;

%% Impart diffraction onto momentum space
Diffract_beam
%% optimize diffraction optics
xi=0.882; %average back focal plane (BFP).
dbfp=0; %distance between BFP's.
Diffraction_optics_optimization
Distance_sample_to_detector=1.5736;
%% Diffraction optics GPT simulation
infile = 'diffraction_optics_setup.in';
gdffile = 'diffraction_GPT_solution.gdf';
statsfile2 = 'diffraction_GPT_stats.gdf'
callGPT = strcat(path_to_gpt,'gpt -v -o ',path_to_in_out_file,gdffile,path_to_in_out_file,infile,...
    ' d1=',num2str(s(1)),...
    ' Spotsize=',num2str(Spotsize),' beamdiv=',num2str(beamdiv),...
    ' d2=',num2str(s(2)),' d3=',num2str(s(3)),' G1=',num2str(G0(1)),' G2=',num2str(G0(2)),' G3=',num2str(G0(3)),...
    ' l1=',num2str(l0(1)),' l2=',num2str(l0(2)),' l3=',num2str(l0(3)),' b1=',num2str(bs(1)),' b2=',num2str(bs(2)),' b3=',num2str(bs(3)),...
    ' Lg=',num2str(xi),...
    ' sg1=',num2str(sg(1)),' sg2=',num2str(sg(2)),' sg3=',num2str(sg(3)),...
    ' Gg1=',num2str(Gg(1)),' Gg2=',num2str(Gg(2)),' Gg3=',num2str(Gg(3)),' Ltot=',num2str(Distance_sample_to_detector),license)
system(callGPT,'-echo') ;
callGPT2=strcat(path_to_gpt,'gdfa -o',path_to_in_out_file,statsfile2,path_to_in_out_file,gdffile,' position ',' stdBy ',' stdBx ',' stdx ',' stdy',' CSgammax',' CSgammay',' CSalphax',' CSalphay', ' CSbetax',' CSbetay',' nemixrms',' nemiyrms',' nemirrms',' avgt',' stdt',' stdz',' avgG',' stdG');
system(callGPT2,'-echo') ;

%% Show diffraction pattern at pimax screen
data=load_gdf('C:\Users\Paul\Desktop\Pegasus_diffraction_simulation_7_18_2023\diffraction_GPT_solution.gdf') ;
% convert to reciprocal space
Rdiff=Mf(:,:,end);
hbar=1.05e-34;
sx=gb*m0*c0*data(end).d.x./hbar./Rdiff(1,2);
sy=gb*m0*c0*data(end).d.y./hbar./Rdiff(3,4);
figure(11)
dscatter(sx*1e-10,sy*1e-10,'BINS',[1500,1500],'PLOTTYPE','surf')
view(2)
