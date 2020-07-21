 function r = get_url(url,params)
 
   global apiKey 
  
   header1 = struct('name','api-key','value',apiKey);
   header2 = struct('name','accept','value','application/json,octet-stream');
   header = struct([header1,header2]);
 
   if exist('params','var')
     keys = fieldnames(params);
     query = '';
     for i=1:numel(keys), query = strcat( query, '&', keys(i), '=',  num2str(params.(keys{i})) ); end
     url = strcat(url, '?', query);
     url = url{1};
   end
 
   [r,extras] = urlread2(url, 'GET', '', header);
 
   if ~extras.isGood, error(['error: ' num2str(extras.status.value) ' ' extras.status.msg]); end
   
   if strcmp(extras.firstHeaders.('Content_Type'),'application/json'), r = loadjson(r); end
  if strcmp(extras.firstHeaders.('Content_Type'),'application/octet-stream')
     global TEMPFILEDIR
      filename = strsplit(extras.firstHeaders.('Content_Disposition'),'filename=');
      fullfilename=[TEMPFILEDIR '/' filename{2}];
     f = fopen(fullfilename,'w');%f = fopen(filename{2},'w');
       fwrite(f,r);
     fclose(f);
     r = fullfilename; %filename{2}; % return filename string
   end
 end