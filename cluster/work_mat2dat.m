global FILE_FORMAT;
FILE_FORMAT = '%s_a1.001L%dMpc.dat';

files = dir('*.mat');
for idx = 1:length(files)
	idx
	load(files(idx).name);
	write_cube_fix(sprintf('%s.dat',files(idx).name(1:end-4)), res);
end
