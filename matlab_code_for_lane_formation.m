%% Description
% This program provides the velocity field for a rotating particle (one
% revolution). 

%% Initialize program

clear all;   %clears variables
warning off


%% sedimentation velocity magnitude

a = 2.0e-6;   %particle radius
g = 9.81;   %gravity
visc = 8.9e-4;   %visocity
densityliq = 997;   %density liquid
densitysol = 1220;   %density solid
deldens = (densitysol-densityliq);   %density difference

%F = (1/6)*a^3*g*deldens/visc;   %sedimentation velocity

F = (4/3)*pi*a^3*g*deldens;   %sedimentation velocity

%% Orbital variables

n = 20; %rotational speed (rpm)

%rsed = vsed*60/(2*pi*n);   %sedimentation radius

%% simulation parameters

rlength =  180e-6;   %simulation space

rotstep1 = -6;   %rotation step, degrees
rotend = -(360+rotstep1);   %end rotation 
rotstart = 0;   %start rotation
timestep = 60*abs(rotstep1)/(n*360);
stoptime = 0.2;

%% Locate center position of random circles that don't overlap
centrx1=dlmread('centrx.txt');
centry1=dlmread('centry.txt');
centrz1=dlmread('centrz.txt');


centrx=centrx1(:,end);
centry=centry1(:,end);
centrz=centrz1(:,end);

%r = sqrt(centrx.^2+centry.^2);
%        nn = find(r >(rlength));   
%        centrx(nn)=[];
%        centry(nn)=[];
%        centrz(nn)=[];
%clear nn

%nn=find(centrz>1.3*a);
%        centrx(nn)=[];
%        centry(nn)=[];
%        centrz(nn)=[];
%clear nn  

centrx_plot(:,1)= centrx;
centry_plot(:,1)= centry;
centrz_plot(:,1)= centrz;
%% Vector calculations

clear centrx1 centry1 centrz1

for run=1:1:300
tic;    
k=0;

		%%4 different rotation in one cycles

for rot = 1:1:4 			
 
            if rot==1
                rotstart = 0;   %start rotation
                rotstep = rotstep1;   %rotation step, degrees
                rotend = -(360+rotstep);   %end rotation 
            elseif rot==2
                rotstart = -360;   %start rotation
                rotstep = -rotstep1;   %rotation step, degrees
                rotend = -(180+rotstep);   %end rotation 
            elseif rot==3
                rotstart = -180;   %start rotation
                rotstep = rotstep1;   %rotation step, degrees
                rotend = -(540+rotstep);   %end rotation 
            else 
                rotstart = -540;   %start rotation
                rotstep = -rotstep1;   %rotation step, degrees
                rotend = -(360+rotstep);   %end rotation 
            end
            t=rotstart;

for u = 1:1:inf
    k=k+1;
    vsedtheta=degtorad(t-90);
    
    vxgrid=zeros(1,size(centrx,1))';
    vygrid=zeros(1,size(centrx,1))';
    vzgrid=zeros(1,size(centrx,1))';
    
    
    %%hardsphere model
        
    centrx1=zeros(1,size(centrx,1))';
    centry1=zeros(1,size(centrx,1))';
    centrz1=zeros(1,size(centrx,1))';
    
    parfor p = 1:1:size(centrx,1)
        
        centrx1(p)=centrx(p);
        centry1(p)=centry(p);
        centrz1(p)=centrz(p);
        
            x = centrx-centrx1(p);     
            y = centry-centry1(p);   
            z = centrz-centrz1(p);
             
            r = sqrt(x.^2+y.^2+z.^2);
             [B,I] = sort(r,'ascend'); % sort distances of all particles from particle p in ascending order
                                       % I is the original index location 

             mmm=find(r(I)<3*a); %
             mmm(1)=[]; % first value is zero, discard that index
    
             if isempty(mmm)
                 ;
             else
                for(j=1:1:size(mmm))
                mm=I(mmm(j));
                
                if((r(mm))<2*a)

                  centrx1(p) = centrx1(p) + (0.5-a/r(mm))*x(mm);
                  centry1(p) = centry1(p) + (0.5-a/r(mm))*y(mm);
                  centrz1(p) = centrz1(p) + (0.5-a/r(mm))*z(mm);
                  
                end
                
                if (centrz1(p)<a)
                centrz1(p)=a;
                end
                
                x = centrx-centrx1(p);     
                y = centry-centry1(p);   
                z = centrz-centrz1(p);
             
                r = sqrt(x.^2+y.^2+z.^2);
                end
                %clear mm
             end
             %clear B I mmm
    
    end
    centrx = centrx1;
    centry = centry1;
    centrz = centrz1;
    clear centrx1 centry1 centrz1
    
    
    clear x y z

%% hard-sphere model for side-wall boundary

r1 = sqrt(centrx.^2+centry.^2);
                    ee=find(r1>(rlength-a));
                    centrx(ee)=(rlength-a).*centrx(ee)./r1(ee);
                    centry(ee)=(rlength-a).*centry(ee)./r1(ee);
                    clear ee  

%% pair-wise colloidal interaction counted for all particles
   
     parfor p = 1:1:size(centrx,1)

            x = centrx(p)-centrx;     
            y = centry(p)-centry;   
            z = centrz(p)-centrz;
            
    [vx,vy,vz] = vsedfield(x,y,z,centrx,centry,centrz,rlength,a,F,vsedtheta);  
            
            vxgrid(p) = sum(vx)+vxgrid(p);
            vygrid(p) = sum(vy)+vygrid(p);
            vzgrid(p) = sum(vz)+vzgrid(p);          

            
     end

% Update new location

       centrx = timestep*vxgrid+centrx;
       centry = timestep*vygrid+centry;
       centrz = timestep*vzgrid+centrz;     


%check if it completes one rotational cycle           

    if(u<=ceil(stoptime/timestep))
        ;
    elseif t~=rotend
        t=t+rotstep;
    else
        break;
    end
end

end

centrx_plot(:,run+1)= centrx;
centry_plot(:,run+1)= centry;
centrz_plot(:,run+1)= centrz;

dlmwrite('centrx_gel.txt',centrx_plot);
dlmwrite('centry_gel.txt',centry_plot);
dlmwrite('centrz_gel.txt',centrz_plot);

toc;
end

function [vx,vy,vz] = vsedfield1(x,y,z,centrx,centry,centrz,rlength,a,F,vsedtheta)
mu = 8.9e-4;
r = sqrt(x.^2+y.^2+z.^2);
h=centrz;
R = sqrt(x.^2+y.^2+(z+2.*h).^2);
z2 = z+2.*h;

r1 = sqrt(centrx.^2+centry.^2);
h1 = sqrt(((rlength./r1-1).*centrx).^2+((rlength./r1-1).*centry).^2);

F1 = F*cos(vsedtheta);
F2 = F*sin(vsedtheta);

  aa=find(centrx>0);
     theta(aa)=pi/2+atan(centry(aa)./centrx(aa));
  bb=find(centrx<0);
     theta(bb)=3*pi/2+atan(centry(bb)./centrx(bb));
  cc=find(centrx==0);
     theta(cc)=0;
     
  theta=theta'; %coordinate trnasformation
     
  F11=F1.*cos(theta)+F2.*sin(theta);  %transform to new coordinate
  F31=-F1.*sin(theta)+F2.*cos(theta);
  F21(1:size(F11),1)=0;
  
  F11=(1-9*a./(16.*h1)).*F11;
  F31=(1-9*a./(8.*h1)).*F31;
  
  F1=F11.*cos(theta)-F31.*sin(theta);
  F2=F11.*sin(theta)+F31.*cos(theta);

% long range- point force- solid wall image reflection
  vx = (F1.*((1./r-1./R)+(x.^2)./r.^3-(x.^2)./R.^3) + ...
      F2.*(x.*y./r.^3-x.*y./R.^3)+ ...
      2*F1.*h.*(h./R.^3-3.*h.*x.^2./R.^5-z2./R.^3+3.*z2.*x.^2./R.^5)+...
      2*F2.*h.*(-3.*h.*x.*y./R.^5+3.*x.*y.*z2./R.^5))/(8*pi*mu);
  
  vy = (F2.*((1./r-1./R)+(y.^2)./r.^3-(y.^2)./R.^3) + ...
      F1.*(x.*y./r.^3-x.*y./R.^3)+ ...
      2*F2.*h.*(h./R.^3-3.*h.*y.^2./R.^5-z2./R.^3+3.*z2.*y.^2./R.^5)+...
      2*F1.*h.*(-3.*h.*x.*y./R.^5+3.*x.*y.*z2./R.^5))/(8*pi*mu);
  
  vz = (F1.*(x.*z./r.^3-x.*z2./R.^3) + ...
      F2.*(y.*z./r.^3-y.*z2./R.^3)+ ...
      2*F1.*h.*(-3.*h.*z2.*x./R.^5+x./R.^3+3.*x.*z2.^2./R.^5)+...
      2*F2.*h.*(-3.*h.*z2.*y./R.^5+y./R.^3+3.*y.*z2.^2./R.^5))/(8*pi*mu);
    
  %short range- source dipole
  A = (F1.*x+F2.*y)./r;
  
  vx1 = (((a^2)./(3*(r.^2))).*(F1-3.*A.*x./r))./(8*pi*mu.*r);
  vy1 = (((a^2)./(3*(r.^2))).*(F2-3.*A.*y./r))./(8*pi*mu.*r);
  vz1 = -(((a^2)./(3*(r.^2))).*(3.*A.*z./r))./(8*pi*mu.*r);
  
  vxx = vx1+vx;
  vyy = vy1+vy;
  vzz = vz1+vz;
  
  %short range- particle rotation
  omega_x = (y.*vzz-z.*vyy)./(a.*r);
  omega_y = (z.*vxx-x.*vzz)./(a.*r);
  omega_z = (x.*vyy-y.*vxx)./(a.*r);
  
  vx2 = (a^3)*(omega_y.*z-y.*omega_z)./r.^3;
  vy2 = (a^3)*(omega_z.*x-z.*omega_x)./r.^3;
  vz2 = (a^3)*(omega_x.*y-x.*omega_y)./r.^3;
  
  ux=vx1+vx2;
  uy=vy1+vy2;
  uz=vz1+vz2;

  kk=find(centrz<2*a);
  ux(kk)=(vx1(kk)+vx2(kk)).*centrz(kk)/(2*a);
  uy(kk)=(vy1(kk)+vy2(kk)).*centrz(kk)/(2*a);
  uz(kk)=(vz1(kk)+vz2(kk)).*centrz(kk)/(2*a);

  vx = vx+ux;
  vy = vy+uy;
  vz = vz+uz;  

          m=find(r==0);
    vx(m) = (1-9*a/(16*h(m))).*F1(m)/(6*pi*mu*a);
    vy(m) = (1-9*a/(16*h(m))).*F2(m)/(6*pi*mu*a);
    vz(m) = 0;
    
    clear m
end


