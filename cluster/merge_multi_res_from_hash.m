function [multix multiy] = merge_multi_res_from_hash(hashx, datax, normx, hashy, datay, normy)

xs = [hashx{1}.(datax)(:)'; ...
 hashx{2}.(datax)(:)'; ...
 hashx{4}.(datax)(:)'; ...
 hashx{8}.(datax)(:)'];

if (length(normx) > 0)
     xs = xs / hashx{8}.(normx);
end

ys = [hashy{1}.(datay)(:)'; ...
 hashy{2}.(datay)(:)'; ...
 hashy{4}.(datay)(:)'; ...
 hashy{8}.(datay)(:)'];

if (length(normy) > 0)
    ys = ys / hashy{8}.(normy);
end
 
[multix multiy] = merge_multi_res(xs, ys);

end