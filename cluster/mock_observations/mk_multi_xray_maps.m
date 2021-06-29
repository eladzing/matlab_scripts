%% generate SZ maps

%% perliminaries
list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
%units;
zemit=0.03;

%%  get Chandra stuff
fname='/home/zinger/workProjects/cluster/xray_code/chandra.area';
%fname='C:\\Users\\eladzing\\Documents\\cluster\\xray_code\\chandra.area';

fid=fopen(fname);
aa=fscanf(fid,'%g');
area=(reshape(aa,[2 length(aa)/2]))';


% unit factor
sterad2arcsec=(180/pi*3600)^2; % convert steradian to arcsec^2
DaOverDl=(1+zemit).^-2;

prjs={'xy' 'yz' 'xz'};

for lm=[4] %8 10 11] %    1:length(list)
    cl=list(lm);
    %cl=6;
    new_env(cl,'a1');
    
    %% generate projections for the boxes
    global NCELL;
    %global hub;
    global CLUSTER;
    global zred;
    
    %mask=true(NCELL,NCELL,NCELL);
    %i1=floor(NCELL/4)+1;i2=i1+floor(NCELL/2)-1;
    
    coTag='full';
    
    
    for k=1:4
        boxx=2^(k-1);
        cellsize=get_cellsize(boxx,'cm');
        unitFac=cellsize*DaOverDl^2/(4*pi*sterad2arcsec);

        % read projection
        for l=1:3
            prj=prjs{l};
            
            proj=xrayProj(prj,boxx,zemit,coTag);
            
            %setup chandra stuff
            area2=interp1(area(:,1),area(:,2),proj.ebins);
            chanArea=zeros(size(proj.data));
            for i=1:size(chanArea,1)
                for j=1:size(chanArea,2)
                    chanArea(i,j,:)=area2;
                end
            end
            
            %fac=cellsize*pi/(360*3600)^2*(1+zemit)^-4;
            
            xRayProj(k).(prj)=trapz(proj.ebins,proj.data.*chanArea,3).*unitFac;
        end
        xRayProj(k).cl=cellsize;
        xRayProj(k).boxx=boxx;
        if k==1
            %mask(i1:i2,i1:i2,i1:i2)=false;
            coTag='co';
        end
    end
    
    clear cub mask proj
    
    %% generate integrated map
    len=xRayProj(end).boxx./xRayProj(1).boxx.*NCELL;
    
    projXY=zeros(len,len);
    projYZ=zeros(len,len);
    projXZ=zeros(len,len);
    
    indShift=(len-[1 2 4 8].*NCELL)/2;
    
    
    for k=1:length(xRayProj)
        
        boxx=xRayProj(k).boxx;
        for i=1:NCELL
            
            indx=indShift(k)+(((i-1)*boxx+1):boxx*i);
            for j=1:NCELL
                
                indy=indShift(k)+(((j-1)*boxx+1):boxx*j);
                
                projXY(indx,indy)=projXY(indx,indy) +  xRayProj(k).xy(i,j);
                projYZ(indx,indy)=projYZ(indx,indy) +  xRayProj(k).yz(i,j);
                projXZ(indx,indy)=projXZ(indx,indy) +  xRayProj(k).xz(i,j);
                
            end
        end
    end
    
    xRayMap(lm).cluster=CLUSTER;
    xRayMap(lm).zred=zred;
    xRayMap(lm).projXY=projXY;
    xRayMap(lm).projYZ=projYZ;
    xRayMap(lm).projXZ=projXZ;
    
    
end
