
function [bird binx biny]=basic_bird(xx,yy,vv,weight,xxlim,yylim,len)

% calculate 2d-histogram on xx-yy plane for values vv
% values have a weight of weight. xxlim yylim are 2-vectors
% with user supplied limits. if both values are zero limits are 
% calculated automaticaly. len is 2-vector of the no. of bins in 
% each axis. if empty then default size is 200X200. binx & biny
% are size of bins needed for histogram normalization



ltm=xx(:);clear xx;
lro=yy(:);clear yy;
wt=weight(:); clear weight
val=vv(:); clear vv;

% prescribed limits
if any(xxlim)
    mnt=xxlim(1);mxt=xxlim(2);
else
    mnt=min(ltm(:));mxt=max(ltm(:));
end

if any(yylim)    
    mnro=yylim(1);mxro=yylim(2);
else
    mnro=min(lro(:));mxro=max(lro(:));
end

indd=find(((lro>mnro)&(lro<mxro))&((ltm>mnt)&(ltm<mxt)));
ltm=ltm(indd);
lro=lro(indd);

if all(len)
    lenx=len(1);leny=len(2);
else
    lenx=200;leny=200;
end

binx=(mxt-mnt)./lenx;
biny=(mxro-mnro)./leny;

bird=zeros(leny,lenx,2);

%indx=ceil(((ltm-mnt)./(mxt-mnt)).*len);
%indy=ceil(((lro-mnro)./(mxro-mnro)).*len); 
for i=1:size(ltm,1)   
    indx=ceil(((ltm(i)-mnt)./(mxt-mnt)).*lenx);
    indy=ceil(((lro(i)-mnro)./(mxro-mnro)).*leny);
    if(indx>=1 && indx<=lenx) && (indy>=1 && indy<=leny)
        bird(indy,indx,1)=bird(indy,indx,1)+val(i).*wt(i);
        bird(indy,indx,2)=bird(indy,indx,2)+wt(i);
    end
end


end