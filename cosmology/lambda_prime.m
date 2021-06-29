function res=lambda_prime(arr,varargin)
% function for randomly selecting a spin parameter (a la Bullock et al
% 2001) value from the log-normal distribution 
% values are taken from Munoz-Cuartas et al 2011
% argument is an array to be filled with randomly selected valus 

lambda0=0.031;lambdaSig=0.57;

i=1;
while i<=length(varargin)
    switch varargin{i}
        case 'lambda0'
            i=i+1;
            lambda0=varargin{i};
            case {'lambdaSig','sigma','sig'}
            i=i+1;
            lambdaSig=varargin{i};
        otherwise
        error('lambda_prime: Illegal argument')
    end
end


res=lognrnd(log(lambda0),lambdaSig,size(arr)); % spin parameter 