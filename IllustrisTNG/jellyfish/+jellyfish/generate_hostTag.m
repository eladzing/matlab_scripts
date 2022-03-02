function galPropStr = generate_hostTag(galPropStr)
% GENERATE_HOSTTAG Add a host tag field to the galaxy property structure 
%   for a unique identifier of hosts across sims and snaps

for j=1:length(galPropStr.hostID)
    galPropStr.hostTag(j)=join([galPropStr.sim(j),"snp0",string(num2str(galPropStr.snap(j))), ...
        "hostID",string(num2str(galPropStr.hostID(j)))],'');
end

