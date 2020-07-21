%check something:

typ='csf';aexp='a1';

list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 24];

for id=1:length(list)
   clname=sprintf('CL%d',list(id));
   new_env(clname,typ,aexp)
   
   global HALO_PATH;
   load(sprintf('%s/virial%d_%s', HALO_PATH, 8,aexp));
   
   global hub
   
   rv=RVIR
   vbox=ceil(2^ceil(log2(2*rv*hub)))
   
   pause
   
end