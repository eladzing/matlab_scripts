 % scratch sheet to try Mo Mao & white model on disks
function [Q rdf]=work_disk_model(p)
%r=1e-3:1e-3:1;

Mv=1e12;
[Rv,~,~,~]=calculate_virials('mv',Mv);

global md
global fg
global fb
global beta
global rh
global rd

md=p(1);%0.1;
fg=p(2);%0.0;
beta=p(3);%1;
lambda=p(4);%0.03;
xi=1;
jd=p(5);%1;
fb=p(6);%0;%$0; %.25;
rh=p(7);%0.2;

rd2=[];
rd0=sqrt(2)*lambda*md/jd*(1+fg+fb)*xi;
dr=1;
i=0;
while dr>1e-5
    i=i+1;
    rd=rd0;
    Q=quad(@get_integ,0,1);
    rd2(i)=sqrt(2)*lambda*md/jd*(1+fg+fb)*Q;
    dr=abs(rd2(i)./rd0-1);
    rd0=rd2(i);
    
end
%i
rdf=rd2(end);

end


function res=get_integ(r)

%global xi
%global jd
global md
global fg
global fb
global beta
%global lambda
global rh
global rd



mdsk=exp_disk_mass(r./rd,1)+fg.*exp_disk_mass(r./rd,beta);
mb=bulge_mass(r./rd,'fb',fb,'hern','rs',rh./rd);
mdisk=md./(1+fg+fb).*(mdsk+mb);


a=1;
b=-r.*(1-md);
c=-r.*mdisk;

ri=(-b+sqrt(b.^2-4.*a.*c))./(2.*a);
%ri2=(-b-sqrt(b.^2-4.*a.*c))./(2.*a);

mfin=ri.^2./r;
%mf2=ri2.^2./r;

vcdm2=(mfin-mdisk)./r;
g=expdisk_accel(r./rd,'fg',fg,'beta',beta);
b=bulge_accel(r./rd,'fb',fb,'hern','rs',rh./rd);
vcdk2=r.*(md/(1+fg+fb)/2/rd^2).*(g+b);

vc2=vcdm2+vcdk2;
integrand=vc2.*expdisk_density(r./rd,'both','fg',fg,'beta',beta).*r.^2./rd.^3;

res=integrand;
end