function [multix multiy] = merge_multi_res(xs, ys)
%xs - the x values, sorted from small x to large x. each row os a different
%resolution, the first row is the finest resolution
CLIPFIRST = 3;
CLIPLAST  = 1;

multix = xs(1, CLIPFIRST+1:end-CLIPLAST);
multiy = ys(1, CLIPFIRST+1:end-CLIPLAST);
for res = 2:size(xs,1)
    idx = find(xs(res,:)>multix(end),1);
    multix = [multix xs(res, idx:end-CLIPLAST)];
    multiy = [multiy ys(res, idx:end-CLIPLAST)];
end


