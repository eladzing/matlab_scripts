


tt=objectTableComp.tag(inds(22)).extractBefore('typ')

iiRand=find(objectTable.tag.contains(tt));
sidRand=objectTable.subject_ids(iiRand)
iiPref=find(objectTablePref.tag.contains(tt));
sidPref=objectTablePref.subject_ids(iiPref)