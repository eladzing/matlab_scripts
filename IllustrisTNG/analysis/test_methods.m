
snap35=4;

bp=illustris.set_env(35,'0000');
fofs35=illustris.groupcat.loadHalos(bp,snap35);
subs35=illustris.groupcat.loadSubhalos(bp,snap35);

subsInfo100 = illustris.infrastructure.build_sub_fof_connection(subs100,fofs100);
centralMask100= subsInfo100.isCentral(gasP100.galMask);

subsInfo35 = illustris.infrastructure.build_sub_fof_connection(subs35,fofs35);
centralMask35= subsInfo35.isCentral(gasP35.galMask);

ms100=gasP100.galMass(gasP100.galMask);
sk100=gasP100.inGal.meanEntMW(gasP100.galMask);

ms35=gasP35.galMass(gasP35.galMask);
sk35=gasP35.inGal.meanEntMW(gasP35.galMask);

ms35_bh=gasP35_bh.galMass(gasP35_bh.galMask);
sk35_bh=gasP35_bh.inGal.meanEntMW(gasP35_bh.galMask);

ms35_km=gasP35_km.galMass(gasP35_km.galMask);
sk35_km=gasP35_km.inGal.meanEntMW(gasP35_km.galMask);




