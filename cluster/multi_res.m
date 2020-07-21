function result = multi_res(SOURCE_EXP, SLICE_EXP)

tmp = eval(sprintf(SOURCE_EXP, 8)); %RHOG(8);
%tmp = eval(strcat('tmp', SLICE_EXP)); %tmp(:,:,128);
%result = reshape(tmp(:), [256 256]);
result = squeeze(eval(sprintf(SLICE_EXP,'tmp')));
imshow(log10(1+result), []); pause(1/2);

for MP = [4 2 1]
    result = imresize(result, 2);
    tmp = eval(sprintf(SOURCE_EXP, MP));
    %tmp = eval(strcat('tmp', SLICE_EXP));
    %tmp = reshape(tmp(:), [256 256]);
    tmp = squeeze(eval(sprintf(SLICE_EXP,'tmp')));
    
    dim = size(result, 1);    
    from = dim/2-128+1;
    to   = from+256-1;
    
    result(from:to, from:to) = tmp;    
    imshow(log10(1+abs(result)), []); pause(1/2);
end

end