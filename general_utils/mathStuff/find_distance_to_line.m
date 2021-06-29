function dist = find_distance_to_line(p0,p1,p2)
%FIND_DISTANCE_TO_LINE - find the distance from a point to a line.
%   Find the distance from point p0 to the line defined by points p1 and p2

if length(p0)<3
    p0=[p0(1) p0(2) 0];
end

if length(p1)<3
    p1=[p1(1) p1(2) 0];
end

if length(p2)<3
    p2=[p2(1) p2(2) 0];
end

dist=sqrt(sum(cross(p0-p1,p0-p2).^2)./sum((p2-p1).^2));

end

