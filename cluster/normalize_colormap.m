function [newimg ticks] = normalize_colormap(img)

newimg = zeros(size(img));
sortedimg = sort(img(:));
jumps = floor(length(sortedimg)/256);


ticks = {};

%last_value = 0;
count = 0;
for idx = 1:jumps:length(sortedimg)
    newimg(img > sortedimg(idx)) = count/255;
    count = count + 1;
    ticks{count} = sprintf('%.2e', sortedimg(idx));
%     last_value = img(idx);
end
