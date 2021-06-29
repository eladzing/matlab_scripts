function full_flux_profile(MPSec,halopath,clustername)

%global FILE_FORMAT;
%FILE_FORMAT = sprintf('%s/%s', halopath, haloformat);
global FILE_FORMAT_SPHERE;
FILE_FORMAT_SPHERE = sprintf('%s/%s', halopath, '%s_sphere_%d.mat');

load(sprintf('%s/virial%d', halopath, 8));
full_mgas = Mgas_sphere(8);
glen=ceil(256*RVIR/R_Profile(length(R_Profile)));

M_gas = zeros(glen,1);
for ridx = 1:glen
    if ridx < 2  
        mg = squeeze(full_mgas(ridx,:,:)).*(R_Profile(ridx));
    else
        mg = squeeze(full_mgas(ridx,:,:)).*(R_Profile(ridx)-R_Profile(ridx-1));
    end
    M_gas(ridx)=sum(mg(:));
end
MGAS=sum(M_gas)
MVIR
m_norm=560*(MVIR./1e13).^0.15.*(MGAS./1e13)  %virial accretion rate in Msun/yr

load(sprintf('%s/virial%d', halopath, MPSec));

full_ff = flux_sphere(MPSec);
%full_ee = energy_sphere(MPSec);
%full_mgas = Mgas_sphere(MPSec);

%IN_FLOW_TEMP = [];
%OUT_FLOW_TEMP = [];
Min_DOT  = [];
Mout_DOT = [];
M_DOT    = [];
%M_gas    = [];


%rp=R_Profile(10:size(R_Profile));

for ridx = 1:length(R_Profile)
    ff = squeeze(full_ff(ridx,:,:));
    %TT = squeeze(full_TT(ridx,:,:));

    %IN_FLOW_TEMP(end+1)  = sum(TT(ff<0).*ff(ff<0))/sum(ff(ff<0));
    %OUT_FLOW_TEMP(end+1) = sum(TT(ff>0).*ff(ff>0))/sum(ff(ff>0));
    Min_DOT(end+1)       = sum(ff(ff<0));
    Mout_DOT(end+1)      = sum(ff(ff>0));
    M_DOT(end+1)         = sum(ff(:));
   
end

l=length(R_Profile);lb=floor(0.01*RVIR*256/R_Profile(length(R_Profile))); dx=0.01*(R_Profile(l)- R_Profile(lb))./RVIR;
%%
used_mdot = 1.0; %M_DOT(RVIR_IDX);


result_dir='/home/titan3/eladzing/cold_flows/printout';

%figure; plot(R_Profile/RVIR, [IN_FLOW_TEMP; OUT_FLOW_TEMP]', '.-'); title('Temperature profile of inflowing and outflowing fluxes');xlabel('r / R_{v}'); ylabel('Mean Temperature (K)'); legend('In Flow', 'Out Flow'); saveas(gcf, sprintf('%s/%s_InOut_Flux_Temperature_%dMPc.png',result_dir,clustername, MPSec));
figure; 
subplot(2,1,1);
%plotyy(R_Profile(lb:l)/RVIR, [M_DOT(lb:l)/abs(used_mdot); Min_DOT(lb:l)/abs(used_mdot); Mout_DOT(lb:l)/abs(used_mdot)], R_Profile(lb:l)/RVIR, [M_DOT(lb:l)/abs(m_norm); Min_DOT(lb:l)/abs(m_norm); Mout_DOT(lb:l)/abs(m_norm)],'.-'); title(sprintf('Mass Flux for Cluster %s in %d Mpc/h box',clustername, MPSec));xlabel('r / R_{v}'); ylabel('Mdot [Msun/yr'); line([min(R_Profile/RVIR) max(R_Profile/RVIR)], [0 0]); legend('Total', 'In', 'Out'); 
%[AX,h1,h2]=plotyy(R_Profile(lb:l)/RVIR,M_DOT(lb:l)/abs(m_norm),R_Profile(lb:l)/RVIR, M_DOT(lb:l)/abs(used_mdot), 'semilogx','semilogx');xlim ([(R_Profile(lb)/RVIR-dx) (R_Profile(l)/RVIR+dx)]);

semilogx(R_Profile(lb:l)/RVIR, M_DOT(lb:l)/abs(m_norm),'.-');xlim ([(R_Profile(lb)/RVIR-dx) (R_Profile(l)/RVIR+dx)]);

%set(get(AX(2),'Ylabel'),'String','dotM [M_{sun}/yr]') 
%set(get(AX(1),'Ylabel'),'String','dotM normalized') 

xlabel('r / R_{v}');ylabel('Mdot normalised'); line([min(R_Profile/RVIR) max(R_Profile/RVIR)], [0 0]);legend(sprintf('M_{norm}=%d',m_norm));

%set(h1,'LineStyle','-');
%set(h2,'LineStyle','-');



title(sprintf('Mass Flux for Cluster %s in %d Mpc/h box',clustername, MPSec));

subplot(2,1,2);
semilogx(R_Profile(lb:l)/RVIR, [M_DOT(lb:l)/abs(used_mdot); Min_DOT(lb:l)/abs(used_mdot); Mout_DOT(lb:l)/abs(used_mdot)]','.-'); xlim ([(R_Profile(lb)/RVIR-dx) (R_Profile(l)/RVIR+dx)])
ylabel('Mdot [Msun/yr]');xlabel('r / R_{v}'); line([min(R_Profile/RVIR) max(R_Profile/RVIR)], [0 0]);legend('Total', 'In', 'Out');



%subplot(2,1,2);
%plot(R_Profile(lb:l)/RVIR, [M_DOT(lb:l)/abs(m_norm); Min_DOT(lb:l)/abs(m_norm); Mout_DOT(lb:l)/abs(m_norm)]', '.-'); xlabel('r / R_{v}'); ylabel('Mdot normalized'); line([min(R_Profile/RVIR) max(R_Profile/RVIR)], [0 0]); legend(sprintf('M_{norm}=%d',m_norm)); 
saveas(gcf,sprintf('%s/%s_fflux_%dMPc.png',result_dir,clustername, MPSec));


%save(sprintf('%s/fluxes%d',halopath, MPSec), 'IN_FLOW_TEMP', 'OUT_FLOW_TEMP', 'M_DOT','Min_DOT', 'Mout_DOT', 'used_mdot');
