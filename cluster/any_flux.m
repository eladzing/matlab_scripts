function any_flux(MPSec,halopath,clustername)

%global FILE_FORMAT;
%FILE_FORMAT = sprintf('%s/%s', halopath, haloformat);
global FILE_FORMAT_SPHERE;
FILE_FORMAT_SPHERE = sprintf('%s/%s', halopath, '%s_sphere_%d.mat');

load(sprintf('%s/virial%d', halopath, 8));

full_ff = flux_sphere(MPSec);
full_TT = T_sphere(MPSec);

IN_FLOW_TEMP = [];
OUT_FLOW_TEMP = [];
Min_DOT  = [];
Mout_DOT = [];
M_DOT    = [];

for ridx = 1:length(R_Profile)
    ff = squeeze(full_ff(ridx,:,:));
    TT = squeeze(full_TT(ridx,:,:));

    IN_FLOW_TEMP(end+1)  = sum(TT(ff<0).*ff(ff<0))/sum(ff(ff<0));
    OUT_FLOW_TEMP(end+1) = sum(TT(ff>0).*ff(ff>0))/sum(ff(ff>0));
    Min_DOT(end+1)       = sum(ff(ff<0));
    Mout_DOT(end+1)      = sum(ff(ff>0));
    M_DOT(end+1)         = sum(ff(:));
end

%%
used_mdot = M_DOT(RVIR_IDX);

result_dir='/home/titan3/eladzing/cold_flows/printout';

figure; plot(R_Profile/RVIR, [IN_FLOW_TEMP; OUT_FLOW_TEMP]', '.-'); title('Temperature profile of inflowing and outflowing fluxes');xlabel('r / R_{v}'); ylabel('Mean Temperature (K)'); legend('In Flow', 'Out Flow'); saveas(gcf, sprintf('%s/%s_InOut_Flux_Temperature_%dMPc.png',result_dir,clustername, MPSec));
figure; plot(R_Profile/RVIR, [M_DOT/abs(used_mdot); Min_DOT/abs(used_mdot); Mout_DOT/abs(used_mdot)]', '.-'); title('Flux decomposition to inflowing and outflowing');xlabel('r / R_{v}'); ylabel('Mdot / Mdot_{v}'); line([min(R_Profile/RVIR) max(R_Profile/RVIR)], [0 0]); legend('Total', 'In', 'Out'); saveas(gcf,sprintf('%s/%s_Mdot_InOut_%dMPc.png',result_dir,clustername, MPSec));

%save(sprintf('%s/fluxes%d',halopath, MPSec), 'IN_FLOW_TEMP', 'OUT_FLOW_TEMP', 'M_DOT','Min_DOT', 'Mout_DOT', 'used_mdot');