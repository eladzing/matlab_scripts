function [range values] = hist_density(data, density, nbins, mindata, maxdata)

data = data(:);
density = density(:);
if (~exist('mindata'))
    maxdata = max(data);
    mindata = min(data);
end
jump = (maxdata-mindata)/nbins;
range = mindata:jump:maxdata-jump;
values = zeros([1 length(range)]);

count = 0;
for minbin = range
    count = count + 1
    mask = bitand(data>=minbin,data<minbin+jump);
    chosen_density = density(mask);
    values(count) = sum(chosen_density);
end

range = range+(jump/2);