epTh=0.5:0.01:1.2;
s=[];

for i=1:length(epTh)
    s(i)=sum(cata.eps<epTh(i));
end

