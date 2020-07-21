%% draw cluster profiles for density, temperature & entropy


printoutdir='printout/profiles';

titletag='%s %s, $z=%d$';



list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
mask=ones(size(list));
%mask =[ 0   0   0   0   0   0   0  0 0 1 0 0 1  0  0  0];

ll=length(list); 
%energyprofs=cell(ll,10);
   
%% global defineitions 
typ='csf';
aexp='a06';

%hubflag='hub';
%vcenterflag='defualt';
%center=[0,0,0];


% kb=1.38066e-16; %erg/K
% mu=0.5926; %mass per particle for primordial composition
% mp=1.672649e-24; % gram
% yr=3.155815e7; %  sec 
% pc=3.0856e18; % center
% km=1e5; %center
% Ms=1.989e33; %gr
% xi=0.5185185; 

%%numerical factors converting to units of M_sun*(km/sec)^2 per Mpc^3
% f_th=1.5.*kb./(mu.*mp)/km.^2 ;%E_th factor converts to units cited above
% f_pr=kb./(mu.*mp)/km.^2 ;%E_th factor converts to units cited above
% f_pv=kb./(mu.*mp.*km.*1.0e6.*pc).*yr; %units of M_sun*(km/sec)^2 per Mpc^3 per yr
% fu=km./(1e6.*pc).*yr;
% f_cl=(xi./(mu.*mp)).^2*yr*Ms/(1e6*pc)^3/km^2*1.e-22;

% in units of Rv 
%rsh=[0.1 1];


for id=1:length(list)
    if(~mask(id)) 
        continue
    end

    
    cluster=sprintf('CL%d',list(id))   

    %begin calculation 
    new_env(cluster,typ,aexp,'win');

    %global HALO_PATH
    %global aexpn;
    global hub;
    global NCELL;   
    global zred;
    global CLUSTER;
        
    cellsize=1/hub/NCELL;
    rv=get_rvir();tv=get_tvir();
    %%ri=1:0.5:  (4/(1/NCELL));
    %%rr=ri.*cellsize;
    rrr=cellsize:cellsize/2:10/hub;
    rrr=rrr./(1+zred);
    
    %[roprof rotot_prof] = read_RHO_Profiles(rrr);
    %tprof = read_T_Profile(rrr);
    sprof = read_S_Profile(rrr);
       
    % find derivatives
    ds=diff(sprof);
    ds1=ds(1:end-1);
    ds2=ds(2:end);
        
    rr=rrr;%./rv;
    dr=diff(rr);
    dr1=dr(1:end-1);
    dr2=dr(2:end);
    ddr=rr(2:end-1);
    
    ddr2=rr(1:end-1)+0.5.*dr;
        
    dsdr=(ds2+ds1)./(dr2+dr1); %1st derivative
    dsdr2=ds./dr;
    ds2dr2=(ds2-ds1)./(dr2.*dr1); %second derivative
    
    % define edge by maximal negative gradient
    s=ds2dr2;
    ss=dsdr(1:end-1);   
    
    m1=(s(1:end-1).*s(2:end)<0 & ss<0);
    
    m2=(ss==min(ss(m1))); 
    ss(m2)=1e10;
    m3=(ss==min(ss(m1)));
    
    f=logical(m2+m3);
    ind=find(f>0)+1;
    l=ddr(ind);
    if ~strcmp(CLUSTER,'CL10')
        l=l(l>1);
    else
        l=l(l>1 & l<10);
    end
        
        
    ll=cat(1,l,l);
    
    % define edge by maxima of entropy 
    sss=sprof(1:end-1);
    mm1=(dsdr2(1:end-1).*dsdr2(2:end)<0);
    mm2=(sss==max(sss(mm1)));
    sss(mm2)=-1;
    mm3=(sss==max(sss(mm1)));
    
    ff=logical(mm2+mm3);
    indd=find(ff>0)+1;
    l2=ddr(indd);
    l2=l2(l2>1); 
    
    ll2=cat(1,l2,l2);
    
       
    
    
    % plot 
    reslim=10.*find_inner_reslim(1,1,3);
    
    xl=[reslim rr(end) ./rv];
    
    name=sprintf('%s/%s_%s_%s.%s',printoutdir,CLUSTER,'%s',aexp,'%s');
       
    
    %figure
    subplot(3,1,1)
    loglog(rr./rv,sprof,'linewidth',2);
    xlim(xl);
    yl=cat(1,ylim,ylim);
    hold on
     semilogx(ll./rv,yl','--r')
     semilogx(ll2./rv,yl','--b')
    hold off
    grid
    ylabelmine('$S\, [\mathrm{Kev\,cm^2}]$')
    titlemine(sprintf(titletag,CLUSTER,'Entropy \& derivatives',zred));
    set(gca,'Fontsize',12)
    
    subplot(3,1,2)
    semilogx(ddr./rv,dsdr,'-','linewidth',2);
    hold on
    semilogx(ddr(f)./rv,dsdr(f),'.r','markersize',20);
    semilogx(ddr(ff)./rv,dsdr(ff),'.b','markersize',20);
    hold off
    xlim(xl);
    grid
    ylabelmine('$dS/dr$')
    set(gca,'Fontsize',12)
    
    subplot(3,1,3)
    semilogx(ddr./rv,ds2dr2,'-','linewidth',2);
    hold on
    semilogx(ddr(f)./rv,ds2dr2(f),'.r','markersize',20);
    semilogx(ddr(ff)./rv,ds2dr2(ff),'.b','markersize',20);
    hold off
    ylabelmine('$d^2S/dr^2$')
    xlabelmine('$r/R_{vir}\, [\mathrm{Mpc}]$')
    xlim(xl);
    grid
    set(gca,'Fontsize',12)
    %exportfig(gcf,sprintf(name,'shockedge','png'),'format','png');
    %exportfig(gcf,sprintf(name,'shockedge','eps'))
    
    %close  all
    %pause
    shockedge{id,1}=CLUSTER;
    shockedge{id,2}=rv;
    shockedge{id,3}=get_mvir;
    
    shockedge{id,4}=l;
    shockedge{id,5}=l2;
    
    
end
