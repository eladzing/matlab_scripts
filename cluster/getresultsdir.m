function str = getresultsdir(filename)

str = sprintf('../results/%s/%s', datestr(now, 'yyyy-mm-dd'), filename);