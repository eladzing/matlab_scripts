function subsInfo = build_sub_fof_connection(subs,fofs)%,snap)
%% prepare list of host halo and central for all the subs
% keep the convention that indices start from 0


% prepare list of host Fof and the central of that fof each sub
hostFof=subs.SubhaloGrNr;
hostCentral=fofs.GroupFirstSub(hostFof+1);

% prepare logical mask for fofs that have all their virial parameters
% defined. 
hostHasVirial=fofs.Group_M_Crit500(hostFof+1)~=0 & ...
    fofs.Group_M_Crit200(hostFof+1)~=0 & ...
    fofs.Group_M_Mean200(hostFof+1)~=0 & ...
    fofs.Group_M_TopHat200(hostFof+1)~=0;


subsInfo.hostFof=hostFof;
subsInfo.hostCentral=hostCentral;
subsInfo.hostHasVirial=hostHasVirial;

% prepare logical list of which subs are centrals
isCentral=false([1 subs.count]);

subIndx=1:subs.count;
isCentral(subIndx==hostCentral+1)=true;
subsInfo.isCentral=isCentral;

% determine if sub has gas/star/dm within the entire sub
subsInfo.hasGas= subs.SubhaloMassType(1,:)>0;
subsInfo.hasDM= subs.SubhaloMassType(2,:)>0;
subsInfo.hasStars= subs.SubhaloMassType(5,:)>0;
subsInfo.hasBH= subs.SubhaloMassType(6,:)>0;
%subsInfo.zred=illustris.utils.get_zred(snap);

subsInfo.DM10perc= subs.SubhaloMassInRadType(2,:)./subs.SubhaloMassInRad>0.1;
subsInfo.DM20perc= subs.SubhaloMassInRadType(2,:)./subs.SubhaloMassInRad>0.2;
subsInfo.DM30perc= subs.SubhaloMassInRadType(2,:)./subs.SubhaloMassInRad>0.3;
subsInfo.DM50perc= subs.SubhaloMassInRadType(2,:)./subs.SubhaloMassInRad>0.5;

end
