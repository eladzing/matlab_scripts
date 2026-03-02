function newCoord = centerObject(coord,center,boxSize)
%centerObject - Center elements of an object (halo,gal etc)
%   This function centers the cooridinates of the elements in an object
%   around the center of the object Given as an argument In addition, it
%   fixes the periodic boundary conditions for the coordinates.
%   The function operates in simulation units: comoving kpc/h.


if ~exist('boxSize','var')
    global LBox
    boxSize=LBox;
end

% partial test to see that no coordinate is beyond the box - can be a
% problem of using the right units. 
if any(any(coord>LBox))
    error('%s - Coordinates found beyond box length. Coordinates should be in simulation units!',current_function().upper);
end



newCoord=coord;


for i=1:length(center)
    
    mask=abs(coord(i,:)-center(i))>0.5.*boxSize;
    
    if sum(mask)>0
        if center(i)>0.5*boxSize
            newCoord(i,mask)=newCoord(i,mask)+boxSize;
        else
            newCoord(i,mask)=newCoord(i,mask)-boxSize;
        end
    end
    
    newCoord(i,:)=newCoord(i,:)-center(i);
end

end

