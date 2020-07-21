function result = vec2dhist(t,ro,xx,weight)
ltm=log10(t);
lro=log10(ro);

wt=ones(size(ltm));
wt=wt.*weight;

len=200;
mnro=floor(min(lro));%  min(min(lro))));
mxro=ceil(max(lro));%  max(max(lro))));
mnt=floor(min(ltm)); %min(min(ltm))));
mxt=ceil(max(ltm));  %max(max(ltm))));
%mnt=mnt-((mxt-mnt)./len);
%mxt=mxt+((mxt-mnt)./len);
tro=zeros(len,len,2);
for i=1:size(ltm,1)   %.*size(ltm,2).*size(ltm,3))
    indx=ceil(((ltm(i)-mnt)./(mxt-mnt)).*len);
    indy=ceil(((lro(i)-mnro)./(mxro-mnro)).*len);  
    tro(indx,indy,1)=tro(indx,indy,1)+xx(i).*wt(i);
    tro(indx,indy,2)=tro(indx,indy,2)+wt(i);
end
result=tro;
%tro(:,1)=ltm(256^3:-1:1);
%tro(:,2)=lro(256^3:-1:1);
%n=hist3(tro,[100 100]);

