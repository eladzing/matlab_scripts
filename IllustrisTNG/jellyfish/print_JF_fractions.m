%% load data
if setupFlag
    global DEFAULT_MATFILE_DIR
    % load object table
    load([ DEFAULT_MATFILE_DIR '\cosmic_jellyfish_objectTable.mat']);
    
    maskJF=objectTable.scoreWeighted>=0.8;
    maskNJF=objectTable.scoreWeighted<=0.2;
    
end


% list snaps and redshifts
snaps=unique(objectTable.snap);
zreds=illustris.utils.snap2redshift(snaps);

%% total numbers
nObject=height(objectTable);
fprintf('Total number of objects = %i \n',nObject);
fprintf('Total number of JF = %i (%4.2f %%) \n',sum(maskJF),sum(maskJF)/nObject*100);

%% break up by simulation33
mask100=objectTable.sim=="TNG100";
mask50=~mask100;
fprintf('TNG100 number of objects = %i (%4.2f %%) \n',sum(mask100),sum(mask100)/nObject*100);
fprintf('TNG50 number of objects = %i (%4.2f %%) \n',sum(mask50),sum(mask50)/nObject*100);

fprintf('TNG100 number of JFs = %i (%4.2f %%) \n',sum(maskJF & mask100),sum(maskJF & mask100)/sum(mask100)*100);
fprintf('TNG50 number of JFs = %i (%4.2f %%) \n',sum(maskJF & mask50),sum(maskJF & mask50)/sum(mask50)*100);

%% break up by redshift
fprintf('========================== \n');
for k=1:3
    switch k
        case 1
            fprintf('redshift range 0<z<0.5 \n')
            snps=snaps(zreds<=0.5);         
        case 2
            fprintf('redshift range 0.5<z<1 \n')
            snps=snaps(zreds>0.5 & zreds<=1);
        case 3
            fprintf('redshift range 1<z<2 \n')
            snps=snaps(zreds>1);
    end
    
    zmask=objectTable.snap>=min(snps) & objectTable.snap<=max(snps); 
    
    
    
    
    
    m100=mask100 & zmask;
    m50=~mask100 & zmask;
    
    fprintf('TNG100 number of objects = %i  \n',sum(m100));
    fprintf('TNG50 number of objects = %i  \n',sum(m50));
    
    
    fprintf('TNG100 number of JFs = %i (%4.2f %%) \n',...
        sum(maskJF & m100),sum(maskJF & m100)/sum(m100)*100);
    fprintf('TNG50 number of JFs = %i (%4.2f %%) \n',...
        sum(maskJF & m50),sum(maskJF & m50)/sum(m50)*100);
    fprintf('========================== \n');
end


