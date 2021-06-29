function any_flux_coldhot(MPSec,halopath,clustername)

%global FILE_FORMAT;
%FILE_FORMAT = sprintf('%s/%s', halopath, haloformat);
global FILE_FORMAT_SPHERE;
FILE_FORMAT_SPHERE = sprintf('%s/%s', halopath, '%s_sphere_%d.mat');

load(sprintf('%s/virial%d', halopath, 8));

full_ff = flux_sphere(MPSec);
full_TT = T_sphere(MPSec);
full_RHOG=RHOG_sphere(MPSec);
full_SS = full_TT./(full_RHOG.^(2./3));


% Using simple threshoolds
Tlim=10^((log10(TVIR)-4.5).*(2./3)+4.5)


full_cold_mask = bitand(full_TT<Tlim,full_RHOG>1e11);

IN_FLOW_TEMP = [];
OUT_FLOW_TEMP = [];

Min_DOT  = [];
Mout_DOT = [];
Mcold_DOT = [];
Mhot_DOT = [];
M_DOT    = [];
RR = [];
Mcold_in = [];
Mcold_out = [];
Mhot_in = [];
Mhot_out = [];

[r_mesh phi_mesh theta_mesh] = sphere_grid(MPSec);
RR = squeeze(r_mesh(:,1,1));


for ridx = 1:length(R_Profile)
    ff = squeeze(full_ff(ridx,:,:));
    TT = squeeze(full_TT(ridx,:,:));

    %IN_FLOW_TEMP(end+1)  = sum(TT(ff<0).*ff(ff<0))/sum(ff(ff<0));
    %OUT_FLOW_TEMP(end+1) = sum(TT(ff>0).*ff(ff>0))/sum(ff(ff>0));
    %Min_DOT(end+1)       = sum(ff(ff<0));
    %Mout_DOT(end+1)      = sum(ff(ff>0));
    %M_DOT(end+1)         = sum(ff(:));

    IN_FLOW_TEMP  = [IN_FLOW_TEMP  sum(TT(ff<0).*ff(ff<0))/sum(ff(ff<0))];
    OUT_FLOW_TEMP = [OUT_FLOW_TEMP sum(TT(ff>0).*ff(ff>0))/sum(ff(ff>0))];
    Min_DOT       = [Min_DOT  sum(ff(ff<0))];
    Mout_DOT      = [Mout_DOT sum(ff(ff>0))];
    M_DOT         = [M_DOT    sum(ff(:))];
    
    cold_mask = squeeze(full_cold_mask(ridx,:,:));
    Mcold_DOT     = [Mcold_DOT   sum(ff(cold_mask))];
    Mhot_DOT      = [Mhot_DOT    sum(ff(~cold_mask))];
    
Mcold_in(end+1) = sum(ff(cold_mask & ff<0));
Mcold_out(end+1) = sum(ff(cold_mask & ff>0));
Mhot_in(end+1) = sum(ff(~cold_mask & ff<0));
Mhot_out(end+1) = sum(ff(~cold_mask & ff>0));



end

%%
used_mdot = M_DOT(RVIR_IDX);

result_dir='/home/titan3/eladzing/cold_flows/printout';

%figure; plot(R_Profile/RVIR, [IN_FLOW_TEMP; OUT_FLOW_TEMP]', '.-'); title('Temperature profile of inflowing and outflowing fluxes');xlabel('r / R_{v}'); ylabel('Mean Temperature (K)'); legend('In Flow', 'Out Flow'); saveas(gcf, sprintf('%s/%s_InOut_Flux_Temperature_%dMPc.png',result_dir,clustername, MPSec));
%figure; plot(R_Profile/RVIR, [M_DOT/abs(used_mdot); Min_DOT/abs(used_mdot); Mout_DOT/abs(used_mdot)]', '.-'); title('Flux decomposition to inflowing and outflowing');xlabel('r / R_{v}'); ylabel('Mdot / Mdot_{v}'); line([min(R_Profile/RVIR) max(R_Profile/RVIR)], [0 0]); legend('Total', 'In', 'Out'); saveas(gcf,sprintf('%s/%s_Mdot_InOut_%dMPc.png',result_dir,clustername, MPSec));


figure(1); plot(RR/RVIR, [IN_FLOW_TEMP; OUT_FLOW_TEMP]', '.-'); title('Temperature profile of inflowing and outflowing fluxes');xlabel('r / R_{v}'); ylabel('Mean Temperature (K)'); legend('In Flow', 'Out Flow');%saveas(gcf, getresultsdir(sprintf('InOutFlux (%d MPSec).png',MPSec)));
figure(2); plot(RR/RVIR, abs([Min_DOT/abs(used_mdot); Mout_DOT/abs(used_mdot)])', '.-'); title('Inflowing and outflowing components of the flux');xlabel('r / R_{v}'); ylabel('Mdot / Mdot_{v}'); legend('In', 'Out');% saveas(gcf, getresultsdir(sprintf('Mdot_InOut (%d MPSec).png',MPSec)));
figure(3); plot(RR/RVIR, [M_DOT/abs(used_mdot)]', '.-'); title('Total flux profile');line([min(RR/RVIR) max(RR/RVIR)], [0 0]); xlabel('r / R_{v}'); ylabel('Mdot / Mdot_{v}');% saveas(gcf, getresultsdir(sprintf('Mdot_Total (%d MPSec).png',MPSec)));
figure(4); plot(RR/RVIR, [Mcold_DOT/abs(used_mdot); Mhot_DOT/abs(used_mdot)]', '.-'); title('Cold flows decomposition of flux');xlabel('r / R_{v}'); ylabel('Mdot / Mdot_{v}'); legend('cold', 'hot');% saveas(gcf, getresultsdir(sprintf('Mdot_ColdHot (%d MPSec).png',MPSec)));



%save(sprintf('%s/fluxes%d',halopath, MPSec), 'IN_FLOW_TEMP', 'OUT_FLOW_TEMP', 'M_DOT','Min_DOT', 'Mout_DOT', 'used_mdot');