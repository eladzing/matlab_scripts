
% eta=0.01:0.01:50;
% epsi=[];
% for i=1:length(cat.sigma)
%     
%     eps(i)=max(sqrt( eta.^2.*(bfunc(eta,1)+...
%     cat.fg(i).*cat.beta(i).^3.*bfunc(eta,cat.beta(i))+...
%     2.*cat.fb(i)./eta.^3)./(1+cat.fg(i)+cat.fb(i))));
% end


eta=0.01:0.01:50;
q=[];
for i=1:length(newcat.sigma)
    
    vc=vcSquare(eta,newcat.fg(i),newcat.beta(i),newcat.fb(i),newcat.xi(i));
    epi=epicyclic_frequency(eta,newcat.fg(i),newcat.beta(i),newcat.fb(i),newcat.xi(i));
    q=QToomre(eta,1,newcat.fg(i),newcat.beta(i),newcat.fb(i),newcat.xi(i));
    
    
    
    qmin(i)=min(q);
    qmax(i)=max(q);
    qsub1(i)=sum(q<1);
    
%     loglog(eta,q)
%     pause



end

