function stream_profile(verts, cube)

h = 0.7;

figure;
hold on;

for k = 1:length(verts)
	vv = verts{k};
	if ~isempty(vv)

        dist = (vv-128.5);
        dist = sqrt(dist(:,1).^2 + dist(:,2).^2 + dist(:,3).^2);

        vv = round(vv);
        inds = sub2ind(size(cube), vv(:,2), vv(:,1), vv(:,3));
        data = cube(inds);

        plot((dist*(8/256)/h)/2.2,log10(data));

%        plot((dist(2:end)+dist(1:end-1))/2*(8/256)/h/2.2, (data(2:end)-data(1:end-1))./data(1:end-1));
%        plot((dist(2:end)+dist(1:end-1))/2*(8/256)/h/2.2, data(2:end)-data(1:end-1));
        %plot(dist);
        pause(1/100);
	end
end

hold off;