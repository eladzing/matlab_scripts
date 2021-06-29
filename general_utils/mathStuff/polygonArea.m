function res=polygonArea(xx,yy)
%% calculate the area of a polygon 

len=length(xx); % number of points in polygon
i=1:len;
ip=i+1-len.*(i+1>len); % periodic i+1;

res=0.5.*abs(sum(yy(i).*xx(ip)-yy(ip).*xx(i)));