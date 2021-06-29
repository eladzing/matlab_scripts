function result = read_cube_fix(filename)
%% a fix to re-order the second index. not sure what this is good for

    result = read_cube(filename);
    result.data = result.data(:,end:-1:1,:);
end
