function res = generateUnitSphereSampling(nSample,plotarg)
%% function for generating uniform sampling points on a unit sphere, such that each
% point represents an equal area on the sphere surface*.
%
%
% * - Areas are almost equal - a difference of ~ 0.5% between maximal and
% are minimal for n=100. 2% of the area is represented by points of area
% <0.998*max(area)
%
plotFlag=false;
if ~exist('nSample','var')
    nSample=100;
end

if exist('plotarg','var')
    plotFlag=strcmp(plotarg,'show');
end


% create uniform sphere of rectangular polygons (with triangles at the
% poles)
[xx,yy,zz] = uni_sphere(nSample);

% find centers of polygons
cnt=1;
for i=1:nSample
    for j=1:nSample
        % extract polygon
        xp=[xx(i,j) xx(i+1,j) xx(i+1,j+1) xx(i,j+1)];
        yp=[yy(i,j) yy(i+1,j) yy(i+1,j+1) yy(i,j+1)];
        zp=[zz(i,j) zz(i+1,j) zz(i+1,j+1) zz(i,j+1)];
        % calculate center point
        
        if( xp(1)==xp(4) && yp(1)==yp(4) && zp(1)==zp(4) )
            xp=xp(1:3);
            yp=yp(1:3);
            zp=zp(1:3);
        end
        
        xc(cnt)=mean(xp);
        yc(cnt)=mean(yp);
        zc(cnt)=mean(zp);
        
        %calculate area
        pp=cat(1,xp,yp,zp)';
        area(cnt) = polygonArea3d(pp);
        cnt=cnt+1;
    end
end

if sum(area)/(4*pi)<0.9
    fprintf('Total area is off by over 10% - consider encreasing sampling')
end

% convert to spherical
%rr=hypot3(xc,yc,zc)
theta=acos(zc./hypot3(xc,yc,zc));
phi=atan2(yc,xc)+pi;

% extend to actual unit sphere.
xs=sin(theta).*cos(phi);
ys=sin(theta).*sin(phi);
zs=cos(theta);

res.theta=theta;
res.phi=phi;

res.xs=xs;
res.ys=ys;
res.zs=zs;

res.area=area;

%% plot diagnostics
if plotFlag
    figure
    tt=tiledlayout(3,2);
    
   
    nexttile([2 2])
    plot3(xs,ys,zs,'x')
    hold on
    surf(xx,yy,zz,'faceAlpha',0.75)
    drawSphere(0,0,0,1,'FaceAlpha',0.5)
    xlabelmine('$X$');
    ylabelmine('$Y$');
    zlabelmine('$Z$');
    axis equal
    
    
    nexttile
    plot(theta/pi,area./(4*pi/length(area)))
    xlabelmine('$\theta/\pi$');
    ylabelmine('$\mathrm{area}/(4\pi/N))$');
    
    %figure
    nexttile
    histogram(area./(4*pi/length(area)),50)
    xlabelmine('$\mathrm{area}/(4\pi/N))$');
    ylabelmine('$N$');
    
    
    
   
      
end