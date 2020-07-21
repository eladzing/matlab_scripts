function mask = f(verts, MIN_DIST)

mask = zeros([256 256 256]);

for k = 1:length(verts);
	vv = verts{k};
	if ~isempty(vv)
        dist = (vv-128.5);
        dist = ceil(min(sqrt(dist(:,1).^2 + dist(:,2).^2 + dist(:,3).^2)));
        if (dist > MIN_DIST) 
            continue;
        end
        
        vv = ceil(vv);
        %inds = sub2ind(size(mask), vv(:,1), vv(:,2), vv(:,3));
        inds = sub2ind(size(mask), vv(:,2), vv(:,1), vv(:,3));
        mask(inds) = 1;
	end
end