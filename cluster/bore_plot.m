%% read 
%list1=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
cl=6;cltype='csf';aexp='a1';
pflag='print';
result_dir='/home/carina/eladzing_a/cold_flows/tit3/printout/';

clustername=sprintf('CL%d',cl)
new_env(clustername,cltype,aexp);

boxx=1;

global aexpn
global HALO_PATH
load(sprintf('%s/virial%d_%s', HALO_PATH, 8,aexpn));
    
%vir=interp1(R_Profile,S_Profile,RVIR,'spline');

SVIR=TVIR.*RVIR.^2./(3.*MVIR./(4.*pi)).^(2./3);
rho_mean=3.*MVIR./(4.*pi*RVIR.^3);


%fluxnorm=(0.056.*(MVIR./1e13).^0.15)./(4.*pi);


ro=RHOG(boxx)./rho_mean;
tm=T(boxx)./TVIR;
s=S(boxx)./SVIR;
cson=sqrt(1.52e8.*tm.*TVIR)./1e5; %speed of sound in km/sec
pre=(ro.*tm); %/(TVIR.*rho_mean);

%%find Vcm
h=0.7;
rvfac=0.2;
rvcm=RVIR*rvfac;
vbox=ceil(2^ceil(log2(2*rvcm*h)));

rc=mk_rcube(vbox,ones(size(ro)));
rv_ind=find(rc<=rvcm);
[VcmX VcmY VcmZ] = V_Vcm_r(vbox,rv_ind); 
clear rv_ind rc vbox;

if ~exist('cm')
    cm = [0,0,0];
end

vx = Vx(boxx)-VcmX;   vy = Vy(boxx)-VcmY;   vz = Vz(boxx)-VcmZ;   
clear VcmX VcmY VcmZ;

%%Vtot=sqrt(Vx.^2+Vy.^2+Vz.^2);
[meshY, meshX, meshZ] = meshgrid(1:size(vy,1), 1:size(vx,2), 1:size(vz,3));
%convert to center origin coordinates
meshX = meshX - (size(vx,1)+1)/2 -cm(1);
meshY = meshY - (size(vy,2)+1)/2 -cm(2);
meshZ = meshZ - (size(vz,3)+1)/2 -cm(3);
% Fix Units (to be in Mpc)
%h=0.7;
meshX = meshX * ((boxx/h)/256);
meshY = meshY * ((boxx/h)/256);
meshZ = meshZ * ((boxx/h)/256);

%rcube2=meshX.^2+meshY.^2+meshZ.^2 ; % r^2 cube in Mpc
%vr=(vx.*meshX+vy.*meshY+vz.*meshZ)./sqrt(rcube2) ; %radial velocity 

%flux=ro.*vr.*rcube2./Mg200.*(1e5./3.0856e24.*3.15e7.*1e9)./fluxnorm ; 

Mach=abs(sqrt(vx.^2+vy.^2+vz.^2)./cson); clear cson;

[gsx gsy gsz]=gradient(s);

grads=sqrt(gsx.^2+gsy.^2+gsz.^2);   
%gradsrad=(gsx.*meshX+gsy.*meshY+gsz.*meshZ)./sqrt(rcube2) ;

clear gsx gsy gsz meshY meshX meshZ; % rcube2;


%%
bp1=-0.1; %[-0.3;-0.25;-0.2;-0.15;-0.1;-0.05;0];
bp2=0;
dir='Y';
borelen=[-boxx/2 boxx/2];tickjump=0.1;
smook=1; % smoothing length in cells:L smoothing is over box of 2k+1
for id1=1:length(bp1)
    for id2=1:length(bp2)
        bind1=ceil((bp1(id1)-(-boxx./2))./(boxx./256));
        bind2=ceil((bp2(id2)-(-boxx./2))./(boxx./256));
        ind1=bind1-smook:1:bind1+smook;
        ind2=bind2-smook:1:bind2+smook;
            
        %% make bores
        bro=mk_bore(ro,ones(size(ro)),ind1,ind2,dir);
        btm=mk_bore(tm,ro,ind1,ind2,dir);
        bs=mk_bore(s,ro,ind1,ind2,dir);
        bmach=mk_bore(Mach,ones(size(Mach)),ind1,ind2,dir);
        bpre=mk_bore(pre,ones(size(pre)),ind1,ind2,dir);
        bgrads=mk_bore(grads,ro,ind1,ind2,dir);
        switch dir
            case {'x','X','1','yz','YZ'}
                vpar=vx; vper=sqrt(vy.^2+vz.^2);
            case {'y','Y','2','zx','ZX'}
                vpar=vy; vper=sqrt(vx.^2+vz.^2);
            case {'z','Z','3','xy','XY'}
                vpar=vz; vper=sqrt(vy.^2+vx.^2);  
        end
        bvpar=mk_bore(vpar./VVIR,ro,ind1,ind2,dir);
        bvper=mk_bore(vper./VVIR,ro,ind1,ind2,dir);
        
        
        %clear ro tm s Mach pre grads
        blen=length(bs);
        iax=1:blen;
        baxis=(boxx./blen).*(iax-0.5)-(boxx./2); 
        
        %% plot bores
        figure;
        %subplot(2,1,1)
        semilogy(baxis,[bro btm bs bpre],'linewidth',3);grid;
        xlim([-0.07 0.44]);ylim([1e-1 1e1]);%;xlim(borelen);ylim([1e-3 1e3]);
        xlabel(sprintf('%s [Mpc/h]',dir),'Fontsize',12);ylabel(sprintf('Profiles'),'Fontsize',12);
        set(gca,'XTick',borelen(1):tickjump:borelen(2),'Fontsize',12)
        %legend('\rho/\rho_{mean}','T/T_{vir}', 'S/S_{vir}','P/(T_{vir}\rho_{mean})','|v|/c_s','|grad(S)|','Location','Southwest');
        legend('\rho/\rho_{mean}','T/T_{vir}', 'S/S_{vir}','P/(T_{vir}\rho_{mean})','Location','Southwest');
        title(sprintf( '%s Bore Profile %s (%0.3g,%0.3g)',clustername,dir,bp1(id1),bp2(id2)),'Fontsize',12);
         if strcmp(pflag,'print')
            saveas(gcf,sprintf('%s/%s_asl_bore%s_%0.3g-%0.3g_sm%d.png',result_dir,clustername,dir,bp1(id1),bp2(id2),2*smook+1));
         end
        
        figure;
        %subplot(2,1,2)
        plot(baxis,[bvpar bvper],'linewidth',3);grid;
        xlim([-0.07 0.44]);%;ylim([5e-3 30]);xlim(borelen);%;ylim([5e-3 30]);
        xlabel(sprintf('%s [Mpc/h]',dir),'Fontsize',12);ylabel(sprintf('Profiles'),'Fontsize',12);
        set(gca,'XTick',borelen(1):tickjump:borelen(2),'Fontsize',12)
        %subplot(2,4,8)
        legend('v_{par}/V_{vir}','v_{per}/V_{vir}','Location','Southwest');
        %set(gcf,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0;1 0 1])
        if strcmp(pflag,'print')
            saveas(gcf,sprintf('%s/%s_asl_vbore%s_%0.3g-%0.3g_sm%d.png',result_dir,clustername,dir,bp1(id1),bp2(id2),2*smook+1));
        end
    end
end

    
