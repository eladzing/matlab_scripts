list1=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
typ='csf';% 'adiabatic');
aexp='a1';

%cvir0=10;

ngrid=256;
max=4./0.7;
min=1./0.7./ngrid;
np=max/min-1;

rmin=[min log10(min)];
rmax=[max log10(max)];
dr=2.*([min ((rmax(2)-rmin(2))./np)]);

rp=rmin(1):dr(1):rmax(1);
lrp=rmin(2):dr(2):rmax(2); 
lrp=10.^lrp;
figure;
for id=1:length(list1)
    clustername=sprintf('CL%d',list1(id));
    new_env(clustername,typ,aexp);
    
     rv=get_rvir();
    mv=get_mvir();  
     mt=read_MTOT_Profile(rp);
    lmt=read_MTOT_Profile(lrp);
    
    
    loglog(rp,mt);
    hold on
    loglog([rmin(1) rmax(1)],[mv mv])
    loglog([rv rv],[mv./1000 mv])
    hold off
    
   pause
end