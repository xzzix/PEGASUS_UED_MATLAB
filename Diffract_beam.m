fileloc='C:\Users\Paul\Desktop\Pegasus_diffraction_simulation_7_18_2023\';
sample_screen_gdf=gdffile;

particles=load_gdf('C:\Users\Paul\Desktop\Pegasus_diffraction_simulation_7_18_2023\injection_solution.gdf') ;

particles1=particles(end);
Npar=length(particles1.d.q);
gb=sqrt(mean(particles1.d.G).^2-1);
hplank=6.63e-34;
m0=9.11e-31;
c0=3e+8;
qe=1.6e-19;
rc=2.8179403262e-15;


a=3e-10; %atomic spacing along axis 1
b=3e-10; %atomic spacing along axis 2
c=13e-10; %atomic spacing along axis 3

angleac=90;
angleab=90;
anglebc=120;


[a1, a2, a3]=get_lattice_vectors(a,b,c,angleac,angleab,anglebc);
% rotate sample
theta=0;
phi=0;
psi=0;
X1=[cos(theta) -sin(theta) 0;
    sin(theta) cos(theta) 0;
    0 0 1];
Z=[1 0 0;
    0 cos(phi) -sin(phi);
    0 sin(phi) cos(phi)]
X2=[cos(psi) -sin(psi) 0;
    sin(psi) cos(psi) 0;
    0 0 1];
Rot=X1*Z*X2;
A=Rot*[a1; a2; a3]';
a1=A(:,1)
a2=A(:,2)
a3=A(:,3) %lattice vectors

V=dot(a3,cross(a1,a2)); %V unit cell

b1=cross(a2,a3)/V;
b2=cross(a3,a1)/V;
b3=cross(a1,a2)/V; %reciprocal lattice vectors

ind=randi([-3 3],[3 Npar]); %uniform rand assigned indiced
xkick=(ind(1,:)*b1(1)+ind(2,:)*b2(1))'*hplank./gb/m0/c0;
ykick=(ind(1,:)*b1(2)+ind(2,:)*b2(2))'*hplank./gb/m0/c0;

particles1.d.Bx=xkick+particles1.d.Bx;
particles1.d.By=ykick+particles1.d.By;
save_struct_to_gdf_file(strcat(fileloc,'matlabOut.gdf'), particles1);

ek=(mean(particles1.d.G)-1)*0.511; %?
figure(101)
dscatter(particles1.d.Bx,particles1.d.By,'BINS',[1500,1500],'PLOTTYPE','surf')
view(2)
% xlim([-1e-3,1e-3])
% ylim([-1e-3,1e-3])