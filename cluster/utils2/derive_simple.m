function [derivative n_range] = derive_simple(data, dx, ORDER, WINDOW_LENGTH)
% Performs a simple discrete first derivative. Works on a +-1 window.
% Works only on a spherical 3D cube (derives along R).
% (ignores ORDER and WINDOW_LENGTH)

derivative = (-data(1:end-2,:,:) + data(3:end,:,:))/(2*dx);
n_range = 2:length(data)-1;