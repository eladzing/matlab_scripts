function [sig rms vbar]=find_dispersion(var,weight)

%% find the dispersion for a given value (velocity conmponent) 
% based on  sigma=sqrt(<v^2>-<v>^2)
% var - is an array of values. 
% weight  - a weighting field;

sig=ones(1,size(var,1));
rms=ones(1,size(var,1));
vbar=ones(1,size(var,1));
% finding <v^2>  - assume v is on  spherical shell? 
for i=1:size(var,1)
    w=squeeze(weight(i,:,:));
    v=squeeze(var(i,:,:));
    vsq=sum(sum(w.*(v.^2)))./(sum(sum(w)));
    vb=sum(sum((w.*v)))./(sum(sum(w)));
    
    vbar(i)=vb;
    sig(i)=sqrt(vsq-vb.^2);
    rms(i)=sqrt(vsq);
end
