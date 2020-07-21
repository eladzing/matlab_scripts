

bp=illustris.set_env(100);

list=[194947,197109,343639,371014,397867,453754,517859,524875];


for i=1:length(list)
        id=list(i)-1
        [galHistory(i),galBirdHistory(i)]=get_galaxy_history_full(id,'zmax',4,'startSnap',99);
        
        [galProf(i),~]=get_halo_history_profiles(id,'zmax',4,'startSnap',99);
        
  
end

global DRACOFLAG
if DRACOFLAG
    global DEFAULT_MATFILE_DIR
    global simDisplayName
    fname=sprintf('list_history_%s',simDisplayName);
    %fnameBird=sprintf('list_history_birds_%s',simDisplayName);
    
    save([DEFAULT_MATFILE_DIR '/' fname],'galHistory','galBirdHistory','galProf','-v7.3')
    fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
    %fprintf(' *** Bird Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fnameBird]);
    
end
fprintf(' *** DONE!  *** \n');

