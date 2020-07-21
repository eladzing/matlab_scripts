load(sprintf('mat_files/stk_eng_profs_%s.mat','a1'),'eng_stack');

id=5;
rp=eng_stack{id,2};
ekr=eng_stack{id,3};
ekt=eng_stack{id,4};
ept=eng_stack{id,5};
eth=eng_stack{id,6};
qcl=eng_stack{id,7};
%[MVIR,RVIR,VVIR,TVIR]=eng_stack{id,8};

ff=find(rp>0.05);

rp1=rp(ff);
ekr1=ekr(ff);
ekt1=ekt(ff);
ept1=ept(ff);
eth1=eth(ff);
qcl1=qcl(ff);


figure
semilogx(rp1,ekr1,rp1,ekt1,rp1,ept1,rp1,eth1)