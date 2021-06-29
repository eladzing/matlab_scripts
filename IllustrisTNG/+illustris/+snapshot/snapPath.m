% Illustris Simulation: Public Data Release.

function [filePath] = snapPath(basePath,snapNum,chunkNum)
  % SNAPPATH  Return absolute path to a snapshot HDF5 file (modify as needed).
  if ~exist('chunkNum','var')
    chunkNum = 0;
  end
  
  k=strfind(basePath,'subbox');
  
  if isempty(k)  % check if we're dealing with a sub-box
      
      snapPath = [basePath '/snapdir_' num2str(snapNum,'%03d') '/'];
      filePath = [snapPath 'snap_' num2str(snapNum,'%03d') '.' num2str(chunkNum) '.hdf5'];
  else
      tag=basePath(k:k+6);
      snapPath = [basePath '/snapdir_' tag '_' num2str(snapNum,'%03d') '/'];
      filePath = [snapPath 'snap_' tag '_' num2str(snapNum,'%03d') '.' num2str(chunkNum) '.hdf5'];
      
  end
  
end
