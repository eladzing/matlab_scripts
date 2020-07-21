function result=mk_exponent_string(value)

% restructures a value to  val times 10^exp in a pretty latex string

% mvType='vir';
% 
% if exist('type','var')
%     mvType=type;
% end


expo=floor(log10(value));

sexpo=int2str(expo);
sm=num2str(value./10^expo,'%4.1f');

result=[sm '\times 10^{' sexpo '}'];

%ll=sprintf('$M_{\\mathrm{%s}}=',mvType);
%result=sprintf('%s %s %s 10^{%s} %s %s',ll,sm,'\times',sexpo,'\,','\mathrm{M_\odot}$');

