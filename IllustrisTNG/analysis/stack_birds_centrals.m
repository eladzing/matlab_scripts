%% load general stuff
if ~skipFlag
    bp=illustris.set_env(100);
    snap=99;
    
    fprintf('*** loading info *** \n');
    
    %subs=illustris.groupcat.loadSubhalos(bp,99);
    fofs=illustris.groupcat.loadHalos(bp,99);
    
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs,snap);
    
    mm=matfile('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/galaxy_birds_z0_TNG100.mat');
end

len=15923;

step=1000;


stop=0;

cnt=0;
while stop<len
    
    % set indexing
    start=stop+1;
    stop=min(stop+step,len);
    
    fprintf('*** start -stop: %s - %s ***\n',num2str(start),num2str(stop))
    
    gbird=mm.galBirds(1,start:stop);
    
    %% initialize
    if start==1
        fprintf('*** Initialize ***\n')
        bside=size(gbird(start).BirdInGal);
        bside=bside(1:2);
        
        
        galBird=zeros(bside);
        cgmBird=zeros(bside);
        subBird=zeros(bside);
        totBird=zeros(bside);
        countA=0;
        
        galBirdC=zeros(bside);
        cgmBirdC=zeros(bside);
        subBirdC=zeros(bside);
        countC=0;
        
        galBirdQC=zeros(bside);
        cgmBirdQC=zeros(bside);
        subBirdQC=zeros(bside);
        countQC=0;
        
        galBirdTC=zeros(bside);
        cgmBirdTC=zeros(bside);
        subBirdTC=zeros(bside);
        countTC=0;
        
        galBirdSFC=zeros(bside);
        cgmBirdSFC=zeros(bside);
        subBirdSFC=zeros(bside);
        countSFC=0;
        
        galBirdS=zeros(bside);
        cgmBirdS=zeros(bside);
        subBirdS=zeros(bside);
        countS=0;
        
        galBirdQS=zeros(bside);
        cgmBirdQS=zeros(bside);
        subBirdQS=zeros(bside);
        countQS=0;
        
        galBirdTS=zeros(bside);
        cgmBirdTS=zeros(bside);
        subBirdTS=zeros(bside);
        countTS=0;
        
        galBirdSFS=zeros(bside);
        cgmBirdSFS=zeros(bside);
        subBirdSFS=zeros(bside);
        countSFS=0;
        
        
    end
    
    for i=1:length(gbird)
        gb=gbird(i);
        gid=gb.id;
        
        if subsInfo.isCentral(gid+1)
            
            cnt=cnt+1;
            %         if gb.GasInCGM==0 || gb.GasInGal==0 || gb.GasInSub==0
            %
            %         end
            if gb.GasInGal.total>0
                bird1=gb.BirdInGal(:,:,1)./max(sum(sum(gb.BirdInGal(:,:,1))),1); %  gb.GasInGal;
            else
                bird1=zeros(200,200);
            end
            
            if gb.GasInCGM.total>0
                bird2=gb.birdInCGM(:,:,1)./max(sum(sum(gb.birdInCGM(:,:,1))),1); %  gb.GasInCGM;
            else
                bird2=zeros(200,200);
            end
            
            if gb.GasInSub.total>0
                bird3=gb.birdInSub(:,:,1)./max(sum(sum(gb.birdInSub(:,:,1))),1); %  gb.GasInSub;
            else
                bird3=zeros(200,200);
            end
            % birdTot=(gb.BirdInGal(:,:,1)+gb.birdInSub(:,:,1))./max(sum(sum(gb.BirdInGal(:,:,1)+gb.birdInSub(:,:,1))),1);
            
            %stack all
            
            galBird=galBird+bird1;
            cgmBird=cgmBird+bird2;
            subBird=subBird+bird3;
            % totBird=totBird+birdTot;
            countA=countA+1;
            
            % find ssfr
            ssfr=subs.SubhaloSFRinRad(gid+1)/gb.stellarMass;
            
            centraldata.stellarMass(cnt)=gb.stellarMass;
            centraldata.ssfr(cnt)=ssfr;
            centraldata.hotGal(cnt)=gb.GasInGal.hotTemp./gb.GasInGal.total;
            centraldata.warnmHotGal(cnt)=gb.GasInGal.warmHot./gb.GasInGal.total;
            centraldata.sfGal(cnt)=gb.GasInGal.starForming./gb.GasInGal.total;
            centraldata.coldGal(cnt)=gb.GasInGal.coldDense./gb.GasInGal.total;
            
            centraldata.hotCGM(cnt)=gb.GasInCGM.hotTemp./gb.GasInCGM.total;
            centraldata.warnmHotCGM(cnt)=gb.GasInCGM.warmHot./gb.GasInCGM.total;
            centraldata.sfCGM(cnt)=gb.GasInCGM.starForming./gb.GasInCGM.total;
            centraldata.coldCGM(cnt)=gb.GasInCGM.coldDense./gb.GasInCGM.total;
            
            centraldata.hotSub(cnt)=gb.GasInSub.hotTemp./gb.GasInSub.total;
            centraldata.warnmHotSub(cnt)=gb.GasInSub.warmHot./gb.GasInSub.total;
            centraldata.sfSub(cnt)=gb.GasInSub.starForming./gb.GasInSub.total;
            centraldata.coldSub(cnt)=gb.GasInSub.coldDense./gb.GasInSub.total;
            
            
            
            
            isSF=ssfr>1e-11;
            isQuenched=ssfr<1e-16;
            isQuenching=~isSF & ~isQuenched;
            
            %         %stack centrals
            %         if subsInfo.isCentral(gid+1)
            
            
            galBirdC=galBirdC+bird1;
            cgmBirdC=cgmBirdC+bird2;
            subBirdC=subBirdC+bird3;
            countC=countC+1;
            
            if isSF
                
                galBirdSFC=galBirdSFC+bird1;
                cgmBirdSFC=cgmBirdSFC+bird2;
                subBirdSFC=subBirdSFC+bird3;
                countSFC=countSFC+1;
                
            elseif isQuenched
                
                galBirdQC=galBirdQC+bird1;
                cgmBirdQC=cgmBirdQC+bird2;
                subBirdQC=subBirdQC+bird3;
                countQC=countQC+1;
                
            else
                galBirdTC=galBirdTC+bird1;
                cgmBirdTC=cgmBirdTC+bird2;
                subBirdTC=subBirdTC+bird3;
                countTC=countTC+1;
            end
            %         else
            %             % stack sats
            %
            %             galBirdS=galBirdS+bird1;
            %             cgmBirdS=cgmBirdS+bird2;
            %             subBirdS=subBirdS+bird3;
            %             countS=countS+1;
            %
            %             if isSF
            %
            %                 galBirdSFS=galBirdSFS+bird1;
            %                 cgmBirdSFS=cgmBirdSFS+bird2;
            %                 subBirdSFS=subBirdSFS+bird3;
            %                 countSFS=countSFS+1;
            %
            %             elseif isQuenched
            %
            %                 galBirdQS=galBirdQS+bird1;
            %                 cgmBirdQS=cgmBirdQS+bird2;
            %                 subBirdQS=subBirdQS+bird3;
            %                 countQS=countQS+1;
            %
            %             else
            %                 galBirdTS=galBirdTS+bird1;
            %                 cgmBirdTS=cgmBirdTS+bird2;
            %                 subBirdTS=subBirdTS+bird3;
            %                 countTS=countTS+1;
            %
            %             end
            %
            %
            %         end
            
            
        end
    end
    
end



fprintf('*** ***\n')



