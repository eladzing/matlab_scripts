function meaned = mean_nan(cube,dim)

no_nan = cube;
no_nan(isnan(cube)) = 0;
counts = sum(~isnan(cube),dim);
meaned = sum(no_nan,dim)./counts;
meaned(isinf(meaned))=NaN;
