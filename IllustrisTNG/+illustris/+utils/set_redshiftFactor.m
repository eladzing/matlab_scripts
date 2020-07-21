function  zFactor=set_redshiftFactor(snap)
%SET_REDSHIFTFACTOR set the factor to fix comoving coordinates 



%global zFactor

aexp=illustris.utils.snap2aexpn(snap);

zFactor.zred=illustris.utils.snap2redshift(snap);
zFactor.length=aexp;
zFactor.density=aexp.^(-3);
%zFactor.velocityPart=sqrt(aexp);


end

