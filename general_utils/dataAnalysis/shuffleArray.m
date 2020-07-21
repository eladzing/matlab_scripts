function res = shuffleArray(arr)
%SHUFFLEARRAY Shuffle an array in random order 
%   
[~,newInd]=sort(rand(size(arr)));
res=arr(newInd);
end

