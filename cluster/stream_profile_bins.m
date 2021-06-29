function [mean_result std_result minresult maxresult sum_result sum_count] = stream_profile_bins(verts, MIN_DIST, cube, MAX_LINES, bins_factor)

if (~exist('bins_factor'))
    bins_factor = 100
end

h = 0.7;
load virial8

num_verts = min(MAX_LINES,length(verts));

bins = zeros([num_verts bins_factor*3]);

for k = 1:num_verts
	vv = verts{k};
	if ~isempty(vv)
        dist = (vv-128.5);
        dist = ceil(min(sqrt(dist(:,1).^2 + dist(:,2).^2 + dist(:,3).^2)));
        if (dist > MIN_DIST) 
            continue;
        end
        
        vv = round(vv);
        %inds = sub2ind(size(mask), vv(:,1), vv(:,2), vv(:,3));
        inds = sub2ind(size(cube), vv(:,2), vv(:,1), vv(:,3));
        
        data = cube(inds);

        dist = (vv-128.5);
        dist = sqrt(dist(:,1).^2 + dist(:,2).^2 + dist(:,3).^2);
        dist = round(dist *8/256/h/RVIR*bins_factor);
        dist(dist==0)=1;
        
%         plot(log10((dist*(8/256)/h)/2.2),log10(data));
%          plot(((dist*(8/256)/h)/2.2),log10(data));
         for idx = 1:length(dist)
             bins(k, dist(idx)) = data(idx);
         end
        
        %plot(dist);
        %pause(1/2);
	end
end

%% average 
nonzeros = (bins ~= 0);
sum_result = sum(bins,1);
sum_count = sum(nonzeros,1);

mean_result = sum_result./sum_count;
bins(bins == 0) = NaN;
minresult = min(bins,[], 1);
maxresult = max(bins,[], 1);

std_result = zeros(size(mean_result));
for i = 1:size(bins,2)
    tmp = bins(:,i);
    tmp = tmp(~isnan(tmp));
    std_result(i) = std(tmp);
end

