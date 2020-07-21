function [meen stnd]=stax(stack,r_p,xind,yind,excl)

slen=size(stack,1);

arr=zeros(slen,length(r_p));

for i=1:slen
    if i~= excl
        arr(i,:)=interp1(stack{i,xind},stack{i,yind},r_p,'spline');
    end
end

meen=mean(arr,1);
stnd=std(arr,0,1);
end


