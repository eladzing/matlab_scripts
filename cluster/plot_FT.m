function plot_FT(ff, TT, R, MPSec)

figure; 
h = 0.8;

semilogx(TT,ff,'.');
xlim([1e4 1e8]);
title(sprintf('F-T   (R=%.2f MPSec)', R/h*(MPSec)/256));
xlabel('T');
ylabel('Flux');

line([1e4 1e8]', [0 0]');
%line([1e5 1e5]', [min(ff(:)) max(ff(:))]');

%COLD_FLOW = sum(ff(TT<1e7))
%HOT_FLOW  = sum(ff(TT>1e7))
IN_FLOW_TEMP  = sum(TT(ff<0).*ff(ff<0))/sum(ff(ff<0))
OUT_FLOW_TEMP = sum(TT(ff>0).*ff(ff>0))/sum(ff(ff>0))
