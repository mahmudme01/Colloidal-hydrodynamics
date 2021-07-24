% This code is to generate non overlapping circular spheres
% near a solid surface
%%
clear all
 a=2e-6;
 rlength=100e-6; %radius
 rheight=4e-6;
 centrx=rand(1,12000)*2*rlength-rlength;
 centrx=centrx';
 centry=rand(1,12000)*2*rlength-rlength;
 centry=centry';
 centrz=rand(1,12000)*(rheight-2*a)+a;

 centrz=centrz';
 rad(1:size(centrx,1),1) = a;

r = sqrt(centrx.^2+centry.^2);
        nn = find(r >(rlength));   
        centrx(nn)=[];
        centry(nn)=[];
        rad(nn)=[];
        centrz(nn)=[];
clear nn

 tic;
 
writerObj = VideoWriter('particle_gen2.avi');
open(writerObj);

f = figure;
% Set a size if desired
    width = 1466;
    height = 1130;
set(f,'Position',[5 5 width height])

%Change the renderer to avoid bugs due to OpenGL
set(f,'Renderer','ZBuffer')

zscaled = centrz(:,end)/1e-6;
cn = ceil(max(zscaled));
cm = colormap(jet(cn));
S = repmat(350,numel(centrx(:,1)),1);

for pp = 1:1:100
    
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
    centrx_check(:,pp+1)=centrx;
    
    %Clear the axes.
    cla;


    %Set the axis aspect ratio to 1:1.
    axis equal;

    %Set a title.
    title([num2str(pp) ' iteration'])
    hold on
    %Display the circles.
    scatter3(centrx, centry,centrz, S,centrz, 'filled')
    colorbar
    caxis([2e-6 20e-6]) %cn*1e-6
    view(0,90)
    grid off
    set(gca,'FontSize',18)
     
      xlim([-rlength-3*a rlength+3*a]);
      ylim([-rlength-3*a rlength+3*a]);
    
    %Pause for 1 second.
    pause(0.1);
    %movievector(i)= getframe;
    frame = getframe(f);
    writeVideo(writerObj,frame);

    if max(centrx_check(:,pp+1)-centrx_check(:,pp))<1e-8
        break;
    end
end
toc;

r = sqrt(centrx.^2+centry.^2);
        nn = find(r >(rlength));   
        centrx(nn)=[];
        centry(nn)=[];
        rad(nn)=[];
        centrz(nn)=[];
clear nn
