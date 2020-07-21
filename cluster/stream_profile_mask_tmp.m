function mask = stream_profile_mask_tmp(verts, MIN_DIST, cube, MAX_LINES)

h = 0.7;

mask = zeros(size(cube));

for k = 1:min(MAX_LINES,length(verts))
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
        
        deltas = data(2:end)-data(1:end-1);
        
        if (max(deltas) > 0.1)
            continue
        end
        
        if (max(data)>0.5)
            continue
        end
        
        mask(inds) = 1;
    end
end
