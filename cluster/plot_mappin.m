function work_quiver_new(clustername, MPSec, clims, usesubplots, texttowrite, rvir)

dilute = 8;
thick = 2;
mapcolormap

cube = S(MPSec);
[Vxx Vyy Vzz VcmX VcmY VcmZ] = V_Vcm(MPSec);

img = zeros([256,256,3]);
for projection = 1:3
    switch projection
        case 1
            img(:,:,projection) = log10(squeeze(sum(cube(129-thick:128+thick,:,:),projection)/(thick*2)));
        case 2
            img(:,:,projection) = transpose(log10(squeeze(sum(cube(:,129-thick:128+thick,:),projection)/(thick*2))));
        case 3
            img(:,:,projection) = transpose(log10(squeeze(sum(cube(:,:,129-thick:128+thick),projection)/(thick*2))));
    end
end
if (length(clims)==0)
     clims(1) = min(img(:));
     clims(2) = max(img(:));
end

figure;

if (usesubplots && exist('texttowrite'))
	subplot(2,2,2);
	title(sprintf('%s data', strrep(clustername, '_','\_')));
	set(gca, 'XTick', [], 'YTick', []);
	text(0.1,0.9, texttowrite, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    caxis(clims);
    colormap(maxcolormap);
    colorbar;
end
for projection = 1:3

    switch projection
        case 1
            v_x = squeeze(mean_nan(Vzz(129-thick:128+thick,:,:),projection));
            v_y = squeeze(mean_nan(Vyy(129-thick:128+thick,:,:),projection));
            if (usesubplots) 
                subplot(2,2,4);    
            else
                figure;
            end
            plot(0,0);
            xlabel('Z (h^{-1}Mpc)');ylabel('Y (h^{-1}Mpc)');title(sprintf('ZY entropy (thickness=%d)',thick*2));
        case 2
            v_y = transpose(squeeze(mean_nan(Vzz(:,129-thick:128+thick,:),projection)));
            v_x = transpose(squeeze(mean_nan(Vxx(:,129-thick:128+thick,:),projection)));
            if (usesubplots)
                subplot(2,2,1);
            else
                figure;
            end

            plot(0,0);
            xlabel('X (h^{-1}Mpc)');ylabel('Z (h^{-1}Mpc)');title(sprintf('XZ entropy (thickness=%d)',thick*2));
        case 3
            v_y = transpose(squeeze(mean_nan(Vyy(:,:,129-thick:128+thick),projection)));
            v_x = transpose(squeeze(mean_nan(Vxx(:,:,129-thick:128+thick),projection)));
            if (usesubplots)
                subplot(2,2,3);
            else
                figure;
            end
            plot(0,0);
            xlabel('X (h^{-1}Mpc)');ylabel('Y (h^{-1}Mpc)');title(sprintf('XY entropy (thickness=%d)',thick*2));
    end


    hold on
    imagesc(-MPSec/2:MPSec/2,-MPSec/2:MPSec/2,squeeze(img(:,:,projection)),clims);
    colormap(maxcolormap);
    
    if ~(usesubplots && exist('texttowrite'))
        colorbar;
    end

    %add quiver plot
    diluted_len = length(1:dilute:256);
    %diluted_jump = round(256/(diluted_len-1));
    %[xx yy] = meshgrid(1:diluted_jump:256,1:diluted_jump:256);
    diluted_jump = MPSec/(diluted_len-1);
    notdiluted_jump = MPSec/(256-1);
    [xx yy] = meshgrid(-MPSec/2:diluted_jump:MPSec/2, -MPSec/2:diluted_jump:MPSec/2);
    resampled_v_x = imresize(v_x, [diluted_len diluted_len], 'bilinear');
    resampled_v_y = imresize(v_y, [diluted_len diluted_len], 'bilinear');
    quiver(xx,yy,resampled_v_x,resampled_v_y, 1.5, 'w');

    %draw circle
    if (exist('rvir'))
        circle = rsmak('circle', rvir, [0 0]);
        lines = fnplt(circle,  1);
        plot(lines(1,:), lines(2,:), 'w', 'LineWidth', 1);
    end

    %add streamlines
    [xx yy] = meshgrid(-MPSec/2:notdiluted_jump:MPSec/2,-MPSec/2:notdiluted_jump:MPSec/2);
    h = streamslice(xx,yy,v_x,v_y,1.5);
    set(h,'color','k')
    %%
    xlim([-MPSec/2 MPSec/2]);ylim([-MPSec/2 MPSec/2]);

    %set(gca,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1], 'XTick', 0:64:256 + 1,'YTick', 0:64:256 + 1, 'TickLength',[-0.015 -0.015]);
    set(gca,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1]);
    tickjump = MPSec/8;
    set(gca,'XTick', -MPSec/2:tickjump:MPSec/2,'YTick', -MPSec/2:tickjump:MPSec/2, 'TickLength',[-0.015 -0.015]);
    
    
    %save plot without subplots
	if (~usesubplots)
		print('-dpng','-r400', getresultsdir(sprintf('%s - Quiver Entropy Plot (%dMpc,dir=%d).png', clustername, MPSec,projection)))
	end
end

%save the big plot with subplots
if (usesubplots)
	print('-dpng','-r400', getresultsdir(sprintf('%s - Quiver Entropy Plot (%dMpc).png',clustername, MPSec)))
end


%%% temp: 
% XY
          subplot(2,2,3);
          imagesc(side,side,squeeze(log10(slice(:,:,3))),clims);%axis equal tight;
          xlabel('X [Mpc/h]','Fontsize',12);
          ylabel('Y [Mpc/h]','Fontsize',12);
set(gca,'Ydir','normal');
set(gca,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1],'XTick', -boxx/2:tickjump:boxx/2,'YTick', -boxx/2:tickjump:boxx/2, 'TickLength',[-0.015 -0.015],'Fontsize',12)
zone_circles(zone);

resampled_v_x = imresize(squeeze(v_x(:,:,3)), [diluted_len diluted_len], 'bilinear');
resampled_v_y = imresize(squeeze(v_y(:,:,3)), [diluted_len diluted_len], 'bilinear');
h = streamslice(xxs,yys,squeeze(v_x(:,:,3)),squeeze(v_y(:,:,3)),1.5);
quiver(xxv,yyv,resampled_v_x,resampled_v_y, 1.5, 'w');
set(h,'color','k')


% ZX
subplot(2,2,1);
imagesc(side,side,squeeze(log10(slice(:,:,2))),clims);%axis equal tight;
set(gca,'Ydir','normal');
%%xlabel('X [Mpc/h]','Fontsize',12);
ylabel('Z [Mpc/h]','Fontsize',12);
set(gca,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1],'XTick', -boxx/2:tickjump:boxx/2,'YTick', -boxx/2:tickjump:boxx/2, 'TickLength',[-0.015 -0.015],'Fontsize',12)
zone_circles(zone);

resampled_v_x = imresize(squeeze(v_x(:,:,2)), [diluted_len diluted_len], 'bilinear');
resampled_v_y = imresize(squeeze(v_y(:,:,2)), [diluted_len diluted_len], 'bilinear');
h = streamslice(xxs,yys,squeeze(v_x(:,:,2)),squeeze(v_y(:,:,2)),1.5);
quiver(xxv,yyv,resampled_v_x,resampled_v_y, 1.5, 'w');
set(h,'color','k')

% ZY
subplot(2,2,4);
imagesc(side,side,squeeze(log10(slice(:,:,1))),clims);%axis equal tight;
set(gca,'Ydir','normal');
xlabel('Z [Mpc/h]','Fontsize',12);
%ylabel('Z [Mpc/h]','Fontsize',12);grid;
set(gca,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1],'XTick', -boxx/2:tickjump:boxx/2,'YTick', -boxx/2:tickjump:boxx/2, 'TickLength',[-0.015 -0.015],'Fontsize',12)
zone_circles(zone);

resampled_v_x = imresize(squeeze(v_x(:,:,1)), [diluted_len diluted_len], 'bilinear');
resampled_v_y = imresize(squeeze(v_y(:,:,1)), [diluted_len diluted_len], 'bilinear');
h = streamslice(xxs,yys,squeeze(v_x(:,:,1)),squeeze(v_y(:,:,1)),1.5);
quiver(xxv,yyv,resampled_v_x,resampled_v_y, 1.5, 'w');
set(h,'color','k')

% title box
subplot(2,2,2);
title(sprintf('%s Entropy Slicw Map Thickness=%s Mpc/h',clustername,num2str(thick,3)),'Fontsize',12);
	set(gca, 'XTick', [], 'YTick', []);
	%text(0.1,0.9, texttowrite, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    caxis(clims);
    colormap(maxcolormap);
    colorbar;


set(gcf,'Colormap',maxcolormap);%'PaperOrientation','landscape')
