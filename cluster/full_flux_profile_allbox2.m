function full_flux_profile_allbox2(halopath,clustername)

%global FILE_FORMAT;
%FILE_FORMAT = sprintf('%s/%s', halopath, haloformat);
global FILE_FORMAT_SPHERE;
FILE_FORMAT_SPHERE = sprintf('%s/%s', halopath, '%s_sphere_%d.mat');

load(sprintf('%s/virial%d', halopath, 2));
full_ff = flux_sphere(2);
M_DOT    = [];

for ridx = 1:(length(R_Profile)-1)
    ff = squeeze(full_ff(ridx,:,:));
    M_DOT(end+1) = sum(ff(:));
end

rp=R_Profile(1:(length(R_Profile)-1));

%%
for ii=2:3
    boxx=2^ii
    load(sprintf('%s/virial%d', halopath, boxx));
    full_ff = flux_sphere(boxx);
    
    for ridx =(length(R_Profile)/2):(length(R_Profile)-1)
        ff = squeeze(full_ff(ridx,:,:));
        M_DOT(end+1) = sum(ff(:));
    end
    rp=cat(2,rp,R_Profile(length(R_Profile)/2:length(R_Profile)-1));
end

ind = find(rp<(0.01*RVIR));
l=length(rp);lb=ind(length(ind))+1; %dx=0.01*(rp(l)- rp(lb))./RVIR;

%floor(0.01*RVIR*256/R_Profile(length(R_Profile))); dx=0.01*(R_Profile(l)- R_Profile(lb))./RVIR;
%%
used_mdot = 1.0; %M_DOT(RVIR_IDX);

MGAS=read_MGAS_Profile(halopath,RVIR)
MVIR
m_norm=-560*(MVIR./1e13).^(0.15).*(MGAS./1e13)  %virial accretion rate in Msun/yr


result_dir='/home/titan3/eladzing/cold_flows/printout';

%figure; plot(R_Profile/RVIR, [IN_FLOW_TEMP; OUT_FLOW_TEMP]', '.-'); title('Temperature profile of inflowing and outflowing fluxes');xlabel('r / R_{v}'); ylabel('Mean Temperature (K)'); legend('In Flow', 'Out Flow'); saveas(gcf, sprintf('%s/%s_InOut_Flux_Temperature_%dMPc.png',result_dir,clustername, MPSec));
figure; 
%subplot(1,2,1);
%plotyy(R_Profile(lb:l)/RVIR, [M_DOT(lb:l)/abs(used_mdot); Min_DOT(lb:l)/abs(used_mdot); Mout_DOT(lb:l)/abs(used_mdot)], R_Profile(lb:l)/RVIR, [M_DOT(lb:l)/abs(m_norm); Min_DOT(lb:l)/abs(m_norm); Mout_DOT(lb:l)/abs(m_norm)],'.-'); title(sprintf('Mass Flux for Cluster %s in %d Mpc/h box',clustername, MPSec));xlabel('r / R_{v}'); ylabel('Mdot [Msun/yr'); line([min(R_Profile/RVIR) max(R_Profile/RVIR)], [0 0]); legend('Total', 'In', 'Out'); 
%[AX,h1,h2]=plotyy(R_Profile(lb:l)/RVIR,M_DOT(lb:l)/abs(m_norm),R_Profile(lb:l)/RVIR, M_DOT(lb:l)/abs(used_mdot), 'semilogx','semilogx');xlim ([(R_Profile(lb)/RVIR-dx) (R_Profile(l)/RVIR+dx)]);

semilogx(rp(lb:l)/RVIR, M_DOT(lb:l)/abs(used_mdot),'.-'); xlim ([(0.8.*rp(lb)/RVIR) (1.2.*rp(l)/RVIR)])

%semilogx(R_Profile(lb:l)/RVIR, M_DOT(lb:l)/abs(m_norm),'.-');xlim ([(R_Profile(lb)/RVIR-dx) (R_Profile(l)/RVIR+dx)]);

%set(get(AX(2),'Ylabel'),'String','dotM [M_{sun}/yr]') 
%set(get(AX(1),'Ylabel'),'String','dotM normalized') 

xlabel('r / R_{v}');ylabel('Mdot [M_{sun}/yr]'); line([min(rp/RVIR) max(rp/RVIR)], [0 0]); line([min(rp/RVIR) max(rp/RVIR)], [m_norm m_norm]);         %legend(sprintf('M_{norm}=%d',m_norm));

%set(h1,'LineStyle','-');
%set(h2,'LineStyle','-');



title(sprintf('Mass Flux for Cluster %s',clustername));

pos=get(gca,'Position');

pos(1)=pos(1)+pos(3)/9.5;
pos(2)=pos(2)+pos(4)/10;
pos(3:4)=pos(3:4)./2.75;
axes('Position',pos);

%subplot(1,2,2);
ind=find(rp>(0.1*RVIR)); le=ind(1)+2;
semilogx(rp(lb:le)/RVIR, M_DOT(lb:le)/abs(used_mdot),'.-'); xlim ([(0.8.*rp(lb)/RVIR) (1.3.*rp(le)/RVIR)])
%semilogx(R_Profile(lb:l)/RVIR, [M_DOT(lb:l)/abs(used_mdot);
%Min_DOT(lb:l)/abs(used_mdot); Mout_DOT(lb:l)/abs(used_mdot)]','.-'); xlim ([(R_Profile(lb)/RVIR-dx) (R_Profile(l)/RVIR+dx)])
%ylabel('Mdot [M_{sun}/yr]');xlabel('r / R_{v}'); 
line([min(rp/RVIR) max(rp/RVIR)], [0 0]);%legend('Total', 'In', 'Out');



%subplot(2,1,2);
%plot(R_Profile(lb:l)/RVIR, [M_DOT(lb:l)/abs(m_norm); Min_DOT(lb:l)/abs(m_norm); Mout_DOT(lb:l)/abs(m_norm)]', '.-'); xlabel('r / R_{v}'); ylabel('Mdot normalized'); line([min(R_Profile/RVIR) max(R_Profile/RVIR)], [0 0]); legend(sprintf('M_{norm}=%d',m_norm)); 
saveas(gcf,sprintf('%s/%s_fullflux.png',result_dir,clustername));


%save(sprintf('%s/fluxes%d',halopath, MPSec), 'IN_FLOW_TEMP', 'OUT_FLOW_TEMP', 'M_DOT','Min_DOT', 'Mout_DOT', 'used_mdot');
