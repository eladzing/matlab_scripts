function result = read_cube_header(filename)
% Utility function to load a "Kravtsov Data Cube"
%
% @param filename  The file to load.
% @returns  A 256^3 matrix of singles.
%
    narginchk(1,1);


    inp = fopen(filename,'rb'); 
    
    JUNK = fread(inp, 1, 'int32');
    
    result.om0    = fread(inp, 1, 'float32');
    result.ol0    = fread(inp, 1, 'float32');
    result.ob0    = fread(inp, 1, 'float32');
    result.hubble = fread(inp, 1, 'float32');
    
    SPACE = fread(inp, 2, 'int32');
    
    result.aexpn = fread(inp, 1, 'float32');
    result.nx    = fread(inp, 3, 'int32');
    %result.nx
    
    SPACE = fread(inp, 2, 'int32');
    
    result.dx = fread(inp, 3, 'single');
    
    %SPACE = fread(inp, 2, 'int32');
    
    %result.data = single(reshape(fread(inp, result.nx(1)*result.nx(2)*result.nx(3), 'float32'), result.nx'));
    
    %JUNK = fread(inp, 1, 'int32');

    fclose(inp);
end

%%%% fortran sourcecode
% fread(&ob0,sizeof(float),1,inp); 
% fread(&hubble,sizeof(float),1,inp); 
% fread(space,sizeof(int),2,inp); 
% fread(&aexpn,sizeof(float),1,inp); 
% fread(nx,sizeof(int),3,inp); 
% fread(space,sizeof(int),2,inp); 
% fread(dx,sizeof(float),3,inp); 
% fread(space,sizeof(int),2,inp); 
% nsize=nx[0]; 
% sv=f3tensor(1,nsize,1,nsize,1,nsize); 
% %/* Read in Fortran output */ 
% for (k=1;k<=nsize;k++) for (j=1;j<=nsize;j++) for (i=1;i<=nsize;i++) 
%    fread(&sv[i][j][k],sizeof(float),1,inp); 
% fread(&junk,sizeof(int),1,inp); 
% fclose(inp); 
% 
% open ( 21, file='rhog_a1.001L2Mpc.dat', form='UNFORMATTED', status='UNKNOWN' ) 
% read(21) om0, oml0, omb0, hubble 
% read(21) aexpn, nx, ny, nz 
% read(21) ddx, ddy, ddz 
% read(21) sv 
% close ( 21 ) 
