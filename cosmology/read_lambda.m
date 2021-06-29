function [temperature, metallicity lambda] = read_lambda()
%read the tabulated data file of the Lambda (cooling function)

sss = dlmread('sam_input_cooling_data.txt', ' ');
lambda = sss(:,2:9);
temperature = sss(:,1);
metallicity = [0 10.^[-3 -2 -1.5 -1 -0.5 0 0.5]];
