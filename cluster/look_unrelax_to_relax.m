%% look for unrelaxed clusters at z=0.6 which are relaxed at z=0 

cl=14;

new_env(cl,'a06')

mkmap(2,'type','flux','proj',[1 1 1],'clims',[-3 3],'marks',0.2)
mkmap(2,'type','rho','proj',[1 1 1],'log','marks',0.2)

new_env(cl,'a1')

mkmap(2,'type','flux','proj',[1 1 1],'clims',[-3 3],'marks',0.2)
mkmap(2,'type','rho','proj',[1 1 1],'log','marks',0.2)

