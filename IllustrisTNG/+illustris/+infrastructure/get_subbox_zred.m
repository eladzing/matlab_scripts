bp=illustris.set_env(100,'sub0');
global subBox
global BASEPATH


    
len=subBox.Nsnaps;

prc5=floor(0.05.*len);
for i=1:len
    sbmax
    head=illustris.snapshot.loadHeader(BASEPATH, i-1);
    
    zred(i)=head.Redshift;
    time(i)=head.Time;
     
    
    if mod(i,prc5)==0
        fprintf('%s %% completed \n',num2str(ceil(i/len*100)));
    end

end
    

global DEFAULT_MATFILE_DIR
global simDisplayName

fname=sprintf('redshift_%s',simDisplayName);
save([DEFAULT_MATFILE_DIR '/' fname],'zred','time','-v7.3')
fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
    

