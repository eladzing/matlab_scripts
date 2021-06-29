function result = lamb22(T,Z)

% approximation of cooling function for high temperatures T>10^6 K
% units - 10^-22 * erg *cm^3 *sec^-1

T6=T./1e6;
Z3=Z./0.03;

result = 0.12.*(Z3.^0.7).*(T6.^-1)+0.02.*T6.^0.5;