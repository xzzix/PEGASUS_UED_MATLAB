function [M_all, zsteps]=quaddiffractionfringeRmatall_sc(s,l,G,b,gamma,L,perveance,sigxx,sigyy,sigxpxp,sigypyp)
   
    
    lq1=l(1);
    lq2=l(2);
    lq3=l(3); %effective lengths
    
    d1 = s(1);
    d2 = s(2);
    d3 = s(3); %quad centers
    
    
    G1=G(1);
    G2=G(2);
    G3=G(3); %gradients
    
    b1=b(1);
    b2=b(2);
    b3=b(3); %field curvature
    
    
    beta = (1-gamma^(-2))^(1/2);
    m0 = 9.10953.*10^(-31);
    e0 = 1.6021892.*10.^(-19);
    c0 = 299792458;
    brho = gamma*beta*m0*c0/e0;
    
    %PMQ field
    Gz=@(z)(zquadrupole(z,G1,b1,lq1,d1)-zquadrupole(z,G2,b2,lq2,d2)+zquadrupole(z,G3,b3,lq3,d3));
    zsteps=linspace(0,L,100000);    
    M_all = zeros(4,4,length(zsteps)); i=1;
    %Transfer Matrices%
    function x=drift(L)
        x=[[1 L 0 0]; [0 1 0 0]; [0 0 1 L]; [0 0 0 1]];
    end
    function x=Quadlensh(kqx,kqy,zstep)
        x=[[cosh(kqx*zstep) 1/kqx*sinh(kqx*zstep) 0 0];
                [kqx*sinh(kqx*zstep) cosh(kqx*zstep) 0 0];
                [0 0 cos(kqy*zstep) 1/kqy*sin(kqy*zstep)];
                [0 0 -kqy*sin(kqy*zstep) cos(kqy*zstep)]];
    end
    function x=Quadlensv(kqx,kqy,zstep)
        x=[[cos(kqx*zstep) 1/kqx*sin(kqx*zstep) 0 0];
                [-kqx*sin(kqx*zstep) cos(kqx*zstep) 0 0];
                [0 0 cosh(kqy*zstep) 1/kqy*sinh(kqy*zstep)];
                [0 0 kqy*sinh(kqy*zstep) cosh(kqy*zstep)]];
    end
    M_all = zeros(4,4,length(zsteps)); i=1;
    Mtot=eye(4);
    for j=1:size(zsteps,2)-1
        dz=zsteps(j+1)-zsteps(j); %small steps through
        spotx=sqrt(Mtot(1,1)^2*sigxx+Mtot(1,2)^2*sigxpxp);
        spoty=sqrt(Mtot(3,3)^2*sigyy+Mtot(3,4)^2*sigypyp);
        SCx=perveance/(spotx*(spotx+spoty));
        SCy=perveance/(spoty*(spotx+spoty));
        if Gz(zsteps(j)) == 0
            Mtot=drift(dz)*Mtot;
        
        elseif Gz(zsteps(j)) < 0
            kqx = (-Gz(zsteps(j))/brho-SCx)^(1/2);
            kqy = (-Gz(zsteps(j))/brho-SCy)^(1/2);
            Mtot=Quadlensh(kqx,kqy,dz)*Mtot;
           
            
        elseif Gz(zsteps(j)) > 0
            kqx = (Gz(zsteps(j))/brho-SCx)^(1/2);
            kqy = (Gz(zsteps(j))/brho-SCy)^(1/2);
            Mtot=Quadlensv(kqx,kqy,dz)*Mtot;
        end
        M_all(:,:,i)=Mtot;
        i=i+1;
    end
    M_all(:,:,i)=Mtot;
    end