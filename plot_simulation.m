clear all
a = 2.0e-6;   %particle radius
rlength = 180e-6;
nfacets=15;

centrx_plot=readmatrix('centrxx12.txt');
centry_plot=readmatrix('centryy12.txt');
centrz_plot=readmatrix('centrzz12.txt');

size1 = size(centrx_plot(1,:),2);

centrx_plot2=readmatrix('centrxx3.txt');
centry_plot2=readmatrix('centryy3.txt');
centrz_plot2=readmatrix('centrzz3.txt');

centrx_plot(:,size(centrx_plot(1,:),2)+1:size(centrx_plot(1,:),2)+size(centrx_plot2(1,:),2))=centrx_plot2;
centry_plot(:,size(centry_plot(1,:),2)+1:size(centry_plot(1,:),2)+size(centry_plot2(1,:),2))=centry_plot2;
centrz_plot(:,size(centrz_plot(1,:),2)+1:size(centrz_plot(1,:),2)+size(centrz_plot2(1,:),2))=centrz_plot2;


centrx_plot(:, size1)=[];
centry_plot(:, size1)=[];
centrz_plot(:, size1)=[];

rad(1:size(centrx_plot,1),1) = a;

%% Create xy VideoWriter object (single color)
% writerObj = VideoWriter('gel1to6_3.avi');
% open(writerObj);
% 
% f = figure;
%  %Set a size if desired
%     width = 1466;
%     height = 1130;
%     set(f,'Position',[0 0 width height])
% 
%  %Change the renderer to avoid bugs due to OpenGL
%     set(f,'Renderer','ZBuffer')
% 
% % zscaled = centrz_plot(:,end)/1e-6;
% % cn = ceil(max(zscaled));
% % cm = colormap(jet(cn));
% %S = repmat(100,numel(centrx_plot(:,1)),1);
% S= rad;
% C= ones(size(centrx_plot,1),1)*[0 0 0];
% 
% [sx,sy,sz]= sphere(nfacets);
% 
% for i=1:1:size(centrx_plot,2)
%     %Clear the axes.
%     cla;
% 
%     %Set the axis aspect ratio to 1:1.
%     axis equal;
% 
%     %Set a title.
%     title([num2str(round((i-1)*9.8/60,1)) ' min'])
% 
%     %Display the circles.
%     %viscircles([centrx_plot(:,i), centry_plot(:,i)], rad,'Color','w');
%     hold on
%     c = gray(256);
%     tic;
%     for j= 1:size(centrx_plot,1)
%     	surface(sx*S(j)+centrx_plot(j,i), sy*S(j)+centry_plot(j,i), sz*S(j)+centrz_plot(j,i),...
% 		'LineStyle','none',...
% 		'FaceColor',C(j,:),...
%         'FaceLighting','gouraud');
% %         colormap(c);
% %         shading interp;
%     end
%     %light('Position',[0 0 1],'Style','local','Color',[1 1 1]);
%     %lighting gouraud;
% 
%     %scatter3sph(centrx_plot(1:10,i), centry_plot(1:10,i),centrz_plot(1:10,i))
%     toc;
% 
%     view(0,90)
%     grid off
%     
%     xlim([-rlength-3*a rlength+3*a]);
%     ylim([-rlength-3*a rlength+3*a]);
%     
%     %Pause for 1 second.
%     pause(0.3);
%     %movievector(i)= getframe;
%     frame = getframe(f);
%     writeVideo(writerObj,frame);
% end  
% 
% close(writerObj);

%% Create xy VideoWriter object (multicolor)
writerObj = VideoWriter('vortexx_multicolr2.avi');
open(writerObj);

f = figure;
% Set a size if desired
    width = 1466;
    height = 1130;
set(f,'Position',[5 5 width height])

%Change the renderer to avoid bugs due to OpenGL
set(f,'Renderer','ZBuffer')

zscaled = centrz_plot(:,end)/1e-6;
cn = ceil(max(zscaled));
cm = colormap(jet(cn));
S = repmat(35,numel(centrx_plot(:,1)),1);

for i=1:1:size(centrx_plot,2)

    %Clear the axes.
    cla;


    %Set the axis aspect ratio to 1:1.
    axis equal;

    %Set a title.
    title([num2str(round((i-1)*3/60,1)) ' min'])
    hold on
    %Display the circles.
%     viscircles([0, 0], rlength,'Color','c');
%     hold on
    %viscircles([centrx_plot(:,i), centrz_plot(:,i)], rad,'Color','b');
    scatter3(centrx_plot(:,i), centry_plot(:,i),centrz_plot(:,i), S,centrz_plot(:,i), 'filled')
    colorbar
    caxis([0 8.5e-6]) %cn*1e-6
    view(0,90)
    grid off
    
     xlim([-rlength-3*a rlength+3*a]);
     ylim([-rlength-3*a rlength+3*a]);
    
    %Pause for 1 second.
    pause(0.1);
    %movievector(i)= getframe;
    frame = getframe(f);
    writeVideo(writerObj,frame);

end  

close(writerObj);



%% Create xz VideoWriter object
% writerObj = VideoWriter('vortezz12.avi');
% open(writerObj);
% 
% f = figure;
% % Set a size if desired
% width = 800;
% height = 600;
% set(f,'Position',[5 5 width height])
% 
% %Change the renderer to avoid bugs due to OpenGL
% set(f,'Renderer','ZBuffer')
% 
% zscaled = centrz_plot(:,end)/1e-6;
% cn = ceil(max(zscaled));
% cm = colormap(jet(cn));
% S = repmat(30,numel(centrx_plot(:,1)),1);
% 
% for i=1:1:size(centrx_plot,2)
% 
%     %Clear the axes.
%     cla;
% 
% 
%     %Set the axis aspect ratio to 1:1.
%     axis equal;
% 
%     %Set a title.
%     title([num2str(round((i-1)*3/60,1)) ' min'])
%     hold on
%     %Display the circles.
% %     viscircles([0, 0], rlength,'Color','c');
% %     hold on
%     %viscircles([centrx_plot(:,i), centrz_plot(:,i)], rad,'Color','b');
%     scatter3(centrx_plot(:,i), centrz_plot(:,i),centrz_plot(:,i), S,centrz_plot(:,i), 'filled')
%     colorbar
%     caxis([0 100e-6]) %cn*1e-6
%     view(0,90)
%     grid off
%     
%     xlim([-rlength-10*a rlength+10*a]);
%     ylim([-10*a 100*a]);
%     
%     %Pause for 1 second.
%     pause(0.1);
%     %movievector(i)= getframe;
%     frame = getframe(f);
%     writeVideo(writerObj,frame);
% 
% end  
% 
% close(writerObj);