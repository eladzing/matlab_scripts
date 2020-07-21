function write_cube(filename, cube_struct)
% Writes a cartesian cube in the same format as read by read_cube().
% cube_struct must contain all fields as returned by read_cube().

	filename    
    inp = fopen(filename,'wb'); 
    
    fwrite(inp, 16, 'int32');
    
    fwrite(inp, cube_struct.om0, 'float32');
    fwrite(inp, cube_struct.ol0, 'float32');
    fwrite(inp, cube_struct.ob0, 'float32');
    fwrite(inp, cube_struct.hubble, 'float32');
    
    fwrite(inp, 16, 'int32'); fwrite(inp, 16, 'int32');
    
    fwrite(inp, cube_struct.aexpn, 'float32');
    fwrite(inp, cube_struct.nx, 'int32');
    
    fwrite(inp, 16, 'int32'); fwrite(inp, 12, 'int32');
    
    fwrite(inp, cube_struct.dx, 'int32');
    
    fwrite(inp, 12, 'int32'); fwrite(inp, 67108864, 'int32');
    
    fwrite(inp, cube_struct.data(:), 'float32');
    
    fwrite(inp, 67108864, 'int32');

    fclose(inp);
end

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
