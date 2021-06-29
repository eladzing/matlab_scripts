%function [qb1 qb2 qb4 qb8] = cooling_profs(cm) 
function [qb1 qb2 qb4 qb8] = cooling_profs(clustername,typ,aexp,cm,cen_fac) 
%function [qb nc lmb voll] = cooling_profs(clustername,typ,aexp,cm,cen_fac) 

[tlam zlam lambda]=read_lambda; %%read interpolation data for cooling function (in log10) 




if ~exist('cm')
    cm = [0,0,0];
end

%%some global definition
Ms=1.989e33; %gr
mu=0.5926; %mass per particle for primordial composition
mp=1.672649e-24; % gram
yr=3.155815e7; %  sec 
pc=3.0856e18; % cm
km=1e5; %cm
xi=0.5185185; 

%%numerical factors
f_cl=(xi./(mu.*mp)).^2*yr*Ms/(1e6*pc)^3/km^2*1.e-22;

new_env(clustername,typ,aexp);    
%    qstack{id,1}=clustername;
    
global hub;
global NCELL;

qb1=2;
qb2=2;
qb4=4;
qb8=4;


for i=[1 4] %1:4  %loop over boxes
    boxx=2^(i-1)
    vol=(boxx./hub./NCELL).^3;
    %load(sprintf('%s/virial%d_%s', HALO_PATH, boxx,aexpn));
    % load relevant cubes
    n=RHOG(boxx);
    tt=log10(T(boxx));
    zmet=ZIa(boxx)+ZII(boxx);
    
    % ignore inner part which in which overcooling takes place
    rc=mk_rcube(boxx,n,cm);
    mask=((rc>cen_fac*get_rvir) & (tt>4.0));  
    clear rc 
            
    lamb=zeros(size(tt));
    len=size(zmet,1);
    
    % create cube of cooling rate
    for ii=1:len    
        zm=squeeze(zmet(ii,:,:));
        tm=squeeze(tt(ii,:,:));  
        lamb(ii,:,:)=interp2(zlam,tlam,lambda,zm,tm,'spline');
        clear tm zm   
    end
        
    lamb2=10.^((lamb+22).*mask).*mask;
    clear mask lamb
    
    coolin=n.^2.*lamb2.*vol.*f_cl; %total energy cooling rate in units Ms*(km/sec)^2 / yr
        
    q=MAKE_PROFILE_FROM_CUBE(coolin);
        
    %lmb=lamb2;
    %voll(i)=vol;
    
    
    clear n coolin
    
    switch i
        case 1
            qb1=q;
        case 2
            
            qb2=q;
            
        case 3
            
            qb4=q;
            
        case 4
            
            qb8=q;
            
    end
    
end
%save(sprintf('mat_files/stk_cool_profs_%s.mat',aexp),'qstack');

