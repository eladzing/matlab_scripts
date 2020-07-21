function new_verts = filter_verts(verts, MIN_DIST, cube, thresh_max_delta, thresh_max)

h = 0.7;

new_verts = {};
new_verts_idx = 1;

for k = 1:length(verts)
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

        deltas = data(2:end)-data(1:end-1);
        
        if (max(deltas) > thresh_max_delta)
            continue
        end
        
        if (max(data)>thresh_max)
            continue
        end
        
        new_verts{new_verts_idx} = verts{k};
        new_verts_idx = new_verts_idx + 1;
        
	end
end

