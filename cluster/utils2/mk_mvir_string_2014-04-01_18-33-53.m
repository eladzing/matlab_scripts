function result=mk_mvir_string(mv)

% restructures the virial mass into a pretty latex string


expo=floor(log10(mv));

sexpo=int2str(expo);
sm=num2str(mv./10^expo,'%4.1f');

ll='$M_{\mathrm{vir}}=';
result=sprintf('%s %s %s 10^{%s} %s %s',ll,sm,'\cdot',sexpo,'\,','\mathrm{M_\odot}$');

