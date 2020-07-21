
function result = zinger2dhist(t,ro,xx,weight)
ltm=log10(t);
lro=log10(ro);

wt=ones(size(ltm));
wt=wt.*weight;
ltm=reshape(ltm,numel(ltm),1);
lro=reshape(lro,numel(lro),1);
wt=reshape(wt,numel(wt),1);
len=200;
mnro=floor(min(lro));%  min(min(lro))));
mxro=ceil(max(lro))-mnro;%  max(max(lro))));
mnt=floor(min(ltm)); %min(min(ltm))));
mxt=ceil(max(ltm))-mnt;  %max(max(ltm))));
%mnt=mnt-((mxt-mnt)./len);
%mxt=mxt+((mxt-mnt)./len);
tro=zeros(len,len,2);
ltm=ltm-mnt;
lro=lro-mnro;
for i=1:size(ltm,1)   %.*size(ltm,2).*size(ltm,3))
    indx=ceil((ltm(i)./mxt).*len);
    indy=len+1-ceil((lro(i)./mxro).*len);  
    tro(indy,indx,1)=tro(indy,indx,1)+xx(i).*wt(i);
    tro(indy,indx,2)=tro(indy,indx,2)+wt(i);
end
result=tro;
