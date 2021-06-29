eps=zeros(size(mms));
qmin=eps;
qsub=eps;

for i=1:length(mms)
    ind=find(res.rr(i,:)./res.rd(i)<15);
    % eps a al Efstathiou et al 1982
    vmax=max(res.vc(i,ind)); % in km/sec
    eps(i)=vmax/sqrt(GG*mms(i)/res.rd(i));
    
    %Q Toomre
    Q=QToomre(res.rr(i,ind),res.vc(i,ind),mms(i),res.rd(i),0.7,zeta(i),GG);
    qmin(i)=min(Q);
    subInd=find(Q<1);
    if ~isempty(subInd)
        qsub(i)=diff(res.rr(i,[subInd(1)+1 subInd(end)+1]))./res.rd(i);
    end
    % loglog(res.rr(i,ind(1)+1:ind(end)-1)./res.rd(i),Q)
    %pause
end
