function result = read_xray_cube_temp(filename)
% Utility function to load an xray projection created with modified apec
% code. 
%
% @param filename  The file to load.
% @returns  A 256^3 matrix of singles.
%
    global NCELL
    ngrid=NCELL; 
    
    narginchk(1,1);


    inp = fopen(filename,'rb'); 
    
    
    
    %JUNK = fread(inp, 1, 'int32');
    
    % to be used after update 
    result.emin   = fread(inp, 1, 'float64');
    result.emax   = fread(inp, 1, 'float64');
    result.eres   = fread(inp, 1, 'float64');
    result.ebins0=result.emin:result.eres:result.emax;
    result.ebins = 0.5.*(result.ebins0(1:end-1)+result.ebins0(2:end));
    ne=length(result.ebins);
    
    
    %result.ebins = fread(inp, 60, 'float64');
    %ne=length(result.ebins);
    %SPACE = fread(inp, 2, 'int32');
    
    %result.aexpn = fread(inp, 1, 'float32');
    %result.nx    = fread(inp, 3, 'int32');
    %result.nx
    
    %SPACE = fread(inp, 2, 'int32');
    
    %result.dx = fread(inp, 3, 'single');
    
    %SPACE = fread(inp, 2, 'int32');
    
    %result.data = single(reshape(fread(inp, result.nx(1)*result.nx(2)*result.nx(3), 'float32'), result.nx'));
    data=fread(inp, ngrid*ngrid*ngrid*ne, 'float64');
    
%     %% testing 
%     data1=single(reshape(data,ne,ngrid,ngrid));
% 
%     data2=zeros(256,256,60);
%    
%    
%     for i=1:256
%         for j=1:256
%             for k=1:60
%                 ind=(i-1)*256*60+(j-1)*60+k;
%                 data2(i,j,k)=data(ind);
%             end
%         end
%     end
%     
%     dd1=squeeze(sum(data1,1));
%     figure;imagesc(log10(dd1))
%     set(gca,'Ydir','normal','Fontsize',14)
%    
%     dd2=sum(data2,3);
%     figure;imagesc(log10(dd2))
%     set(gca,'Ydir','normal','Fontsize',14)
%     
    result.data = single(reshape(data,ne,ngrid,ngrid,ngrid));
    %JUNK = fread(inp, 1, 'int32');

    fclose(inp);
end

