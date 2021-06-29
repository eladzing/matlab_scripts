function ri=solve_MMW_NFW_one(r,Mv,Rv,c,Ms,rd,fg,beta,Mb,xi)

%ri=zeros(size(r));
%for i=1:length(r);

try
   ri=fzero(@(x)myfunc(x,r,Mv,Rv,c,Ms,rd,fg,beta,Mb,xi),r);
catch
    fprintf('r= %f \n',r);
    fprintf('Mv= %f \n',Mv);
    fprintf('Rv= %f \n',Rv);
    fprintf('c= %f \n',c);
    fprintf('Ms= %f \n',Ms);
    fprintf('rd= %f \n',rd);
    fprintf('fg= %f \n',fg);
    fprintf('beta= %f \n',beta);
    fprintf('Mb= %f \n',Mb);
    fprintf('xi= %f \n',xi);
    
    
end
   
   %ri(i)=x;
%end
    

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