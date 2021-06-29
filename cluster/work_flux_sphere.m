%% Init
global FILE_FORMAT;
FILE_FORMAT = '../data/%s_a1.001L%dMpc.dat';

%%
MPSec = 8;load virial8
full_ff = flux_sphere(MPSec);
full_TT = T_sphere(MPSec);
full_SS = S_sphere(MPSec);

% Using simple threshoolds
%full_cold_mask = bitand(full_TT<1.5e7,full_SS<0.135);
%using birds:
load ('/home/alf/giorae/results/2008-01-10/coldflows/parametric 2008-01-10 16-20-49.mat', 'xv','yv')
full_cold_mask = inpolygon(log10(full_TT), log10(RHOG_sphere(MPSec)),xv,yv);


%%
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

for ridx = 1:length(RR)
    ff = squeeze(full_ff(ridx,:,:));
    TT = squeeze(full_TT(ridx,:,:));

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
used_mdot = M_DOT(RVIR_IDX{MPSec});
%  used_mdot = MDOTVIR;

figure(1); plot(RR/RVIR, [IN_FLOW_TEMP; OUT_FLOW_TEMP]', '.-'); title('Temperature profile of inflowing and outflowing fluxes');xlabel('r / R_{v}'); ylabel('Mean Temperature (K)'); legend('In Flow', 'Out Flow'); saveas(gcf, getresultsdir(sprintf('InOutFlux (%d MPSec).png',MPSec)));
figure(2); plot(RR/RVIR, abs([Min_DOT/abs(used_mdot); Mout_DOT/abs(used_mdot)])', '.-'); title('Inflowing and outflowing components of the flux');xlabel('r / R_{v}'); ylabel('Mdot / Mdot_{v}'); legend('In', 'Out'); saveas(gcf, getresultsdir(sprintf('Mdot_InOut (%d MPSec).png',MPSec)));
figure(3); plot(RR/RVIR, [M_DOT/abs(used_mdot)]', '.-'); title('Total flux profile');line([min(RR/RVIR) max(RR/RVIR)], [0 0]); xlabel('r / R_{v}'); ylabel('Mdot / Mdot_{v}'); saveas(gcf, getresultsdir(sprintf('Mdot_Total (%d MPSec).png',MPSec)));
figure(4); plot(RR/RVIR, [Mcold_DOT/abs(used_mdot); Mhot_DOT/abs(used_mdot)]', '.-'); title('Cold flows decomposition of flux');xlabel('r / R_{v}'); ylabel('Mdot / Mdot_{v}'); legend('cold', 'hot'); saveas(gcf, getresultsdir(sprintf('Mdot_ColdHot (%d MPSec).png',MPSec)));


figure(5); plot(RR/RVIR, abs([Mcold_in; Mcold_out; Mhot_in; Mhot_out]')/abs(used_mdot), '.-'); title('Cold flows decomposition of flux');xlabel('r / R_{v}'); ylabel('Mdot / Mdot_{v}'); legend('Cold in', 'Cold out', 'Hot in', 'Hot out'); saveas(gcf, getresultsdir(sprintf('Mdot_ColdHot inout (%d MPSec).png',MPSec)));
figure(6); loglog(RR/RVIR, abs([Mcold_in; Mcold_out; Mhot_in; Mhot_out]')/abs(used_mdot), '.-'); title('Cold flows decomposition of flux');xlabel('r / R_{v}'); ylabel('Mdot / Mdot_{v}'); legend('Cold in', 'Cold out', 'Hot in', 'Hot out'); saveas(gcf, getresultsdir(sprintf('Mdot_ColdHot inout (%d MPSec loglog).png',MPSec)));
figure(7); semilogy(RR/RVIR, abs([Mcold_in; Mcold_out; Mhot_in; Mhot_out]')/abs(used_mdot), '.-'); title('Cold flows decomposition of flux');xlabel('r / R_{v}'); ylabel('Mdot / Mdot_{v}'); legend('Cold in', 'Cold out', 'Hot in', 'Hot out'); saveas(gcf, getresultsdir(sprintf('Mdot_ColdHot inout (%d MPSec semilogy).png',MPSec)));

%% "AITOFF"
