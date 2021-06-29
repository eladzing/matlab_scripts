function ri=solve_MMW_NFW_arr(r,Mv,Rv,c,Ms,rd,fg,beta,Mb,xi)

sz=size(r);
len=sz(2);
% MvA=repmat(Mv,1,len);
% RvA=repmat(Rv,1,len);
% cA=repmat(c,1,len);
% MsA=repmat(Ms,1,len);
% rdA=repmat(rd,1,len);
% fgA=repmat(fg,1,len);
% betaA=repmat(beta,1,len);
% MbA=repmat(Mb,1,len);
% xiA=repmat(xi,1,len);

%ri=arrayfun(@(i) fzero(@(x)myfunc(x,r,MvA,RvA,cA,MsA,rdA,fgA,betaA,MbA,xiA),r(i)),...
%    1:numel(r));


%arrayfun(@(i) fzero(@(x) minme(y(i),x),1),1:numel(y))

%ri=zeros(size(r));
for j=1:sz(1)
    for i=1:sz(2);
        x=fzero(@(x)myfunc(x,r(j,i),Mv(j),Rv(j),c(j),Ms(j),rd(j),fg(j),beta(j),Mb(j),xi(j)),r(j,i));
        ri(j,i)=x;
    end
end

end

function f=myfunc(x,r,Mv,Rv,c,Ms,rd,fg,beta,Mb,xi)

mb=Mb./Mv;
md=Ms.*(1+fg)./Mv;

A=Mv./(log(1+c)-c./(1+c));
B=c./Rv;
C=r.*(Ms.*(1-exp(-1.*r./rd).*(1+r./rd)+fg.*(1-exp(-1.*r./rd.*beta).*(1+r./rd.*beta)))+...
    Mb.*r.^2./(rd./xi+r).^2);

f=A.*(log(1+B.*x)-B.*x/(1+B.*x)).*(x-(1-md-mb).*r)-C;

end