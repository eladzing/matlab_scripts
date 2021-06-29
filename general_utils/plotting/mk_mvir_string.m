function result=mk_mvir_string(mv,type)

% restructures the virial mass into a pretty latex string

mvType='vir';

if exist('type','var')
    mvType=type;
end


expo=floor(log10(mv));

sexpo=int2str(expo);
sm=num2str(mv./10^expo,'%4.1f');

ll=sprintf('$M_{\\mathrm{%s}}=',mvType);
result=sprintf('%s %s %s 10^{%s} %s %s',ll,sm,'\times',sexpo,'\,','\mathrm{M_\odot}$');

