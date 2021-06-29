function ri=solve_MMW_NFW(r,Mv,Rv,c,Ms,rd,fg,beta,Mb,xi)

ri=zeros(size(r));
for i=1:length(r);
   x=fzero(@(x)myfunc(x,r(i),Mv,Rv,c,Ms,rd,fg,beta,Mb,xi),r(i));
   
   ri(i)=x;
end
    

end

function f=myfunc(x,r,Mv,Rv,c,Ms,rd,fg,beta,Mb,xi)

mb=Mb/Mv;
md=Ms*(1+fg)/Mv;

A=Mv/(log(1+c)-c/(1+c));
B=c/Rv;
C=r*(Ms*(1-exp(-1*r/rd)*(1+r/rd)+fg*(1-exp(-1*r/rd*beta)*(1+r/rd*beta)))+...
    Mb*r^2/(rd/xi+r)^2);

f=A*(log(1+B*x)-B*x/(1+B*x))*(x-(1-md-mb)*r)-C;

end