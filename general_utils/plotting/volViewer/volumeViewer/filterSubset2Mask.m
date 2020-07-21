% Converts a filter subset string into a mask. If the filter string is
% invalid it is manipulated to form a valid string which is returned in 
% place of a mask.
%
% Usage:
% ------
% mask = filterSubset2Mask(fString) - Assumes a volume large enough to
%     include all the slices and volumes.
% mask = filterSubset2Mask(fString, dim) - Create a volume with given total
%     dimensions, or large enough to include all slices and volumes
%     whichever is larger. Specifying dim helps resolve : values in
%     fString.
%
% Examples of fstring are:
% (1-2,3-4,:) - volume with rows 1-2, columns 3-4, and all slices (!!! dim
% must be specified to resolve :)
% X10, S10, H10 - all data in row 10 (saggital, horizontal)
% Y5, C5, V5 - all data in column 5 (coronal, vertical)
% Z30, T30, SH30 - all data in slice 30 (transaxial, short-axis)
% A combination of these sets may be used such as:
% (1-40,10:120,:),X10 will give rows 1 to 40 cropped to columns 10 to 120 
% and all slices and row 10 will be displayed.
%
% See also: filterString2Mask

% By Ran Klein, University of Ottawa Heart Institute, 20-Oct-2005


function mask = filterSubset2Mask(fString, dim)


if nargin<2
	dim = nan;
	nodim = true;
elseif length(dim)~=3
	msgId = sprintf('FADS:%s:wrongParameterDimensions', mfilename);
	error(msgId, 'Dim must be an array of 3 integers.');
else
	nodim = false;
end

nonvalid = false;  % Assume valid filter string
fString = upper(deblank(fString)); % Ignore all blanks and make capital letters
fString = strrep(fString,'-',':'); % ##-## and ##:## formats are identical

if isempty(fString)
	if nodim
		msgId = sprintf('FADS:%s:wrongNumberOfInputParameters', mfilename);
		error(msgId, 'Cannot resolve empty subset string if no dimensions are provided.');
	end
	mask = {1:dim(1) 1:dim(2) 1:dim(3)};
else	
	valid = ((fString>='0' & fString<='9') | fString==':' | fString==',' |...
		fString=='(' | fString==')' | ...
		fString=='T' | fString=='S' | fString=='C' | fString=='V' | fString=='H' |...
		fString=='X' | fString=='Y' | fString=='Z'); % The valid character set
	if any(~valid) % Remove non valid characters from the string
		nonvalid = true;
		fString = fString(valid);
	end

	mask = {};
	corrString = [];
	if fString(1)=='(' % String specified a cropped volume
		if fString(end) ==')' % remove the brackets
			fString = fString(2:end-1);
		else
			fString = fString(2:end);
		end
		i = find(fString==',');
		if length(i)~=2
			nonvalid = true;
		else
			X = str2double(fString(1:i(1)-1));
			Y = str2double(fString(i(1)+1:i(2)-1));
			Z = str2double(fString(i(2)+1:end));
			corrString = '(';
			if isempty(X)
				if nodim, nonvalid = true; else	mask{1,1} = 1:dim(1); end
			elseif isnumeric(X)
				if X(1)<1, X = 1:max(1,X(end)); nonvalid=true; end
				if ~nodim & X(end)>dim(1), X = min(X(1),dim(1)):dim(1); nonvalid=true; end
				if length(X)>1
					corrString = [corrString num2str(X(1)) '-' num2str(X(end))];
				else
					corrString = [corrString num2str(X)];
				end
				mask{1,1} = X;
			else
				nonvalid = true;
			end
			if isempty(Y)
				if nodim, nonvalid = true; else	mask{1,2} = 1:dim(2); end
			elseif isnumeric(Y)
				if Y(1)<1, Y = 1:max(1,Y(end)); nonvalid=true; end
				if ~nodim & Y(end)>dim(2), Y = min(Y(1),dim(2)):dim(2); nonvalid=true; end
				if length(Y)>1
					corrString = [corrString ',' num2str(Y(1)) '-' num2str(Y(end))];
				else
					corrString = [corrString ',' num2str(Y)];
				end
				mask{1,2} = Y;
			else
				nonvalid = true;
			end
			if isempty(Z)
				if nodim, nonvalid = true; else	mask{1,3}= 1:dim(3); end
			elseif isnumeric(Z)
				if Z(1)<1, Z = 1:max(1,Z(end)); nonvalid=true; end
				if ~nodim & Z(end)>dim(3), Z = min(Z(1),dim(3)):dim(3); nonvalid=true; end
				if length(Z)>1
					corrString = [corrString ',' num2str(Z(1)) '-' num2str(Z(end))];
				else
					corrString = [corrString ',' num2str(Z)];
				end
				mask{1,3} = Z;
			else
				nonvalid = true;
			end
			corrString = [corrString ')'];
		end % 2 commas found
	else % String specifies 2D slices
		
		[token, fString] = strtok(fString,',');
		while ~isempty(token) | ~isempty(fString)
			if ~isempty(findstr(token,'('))
				[token1 fString] = strtok(fString,',');
				token = [token ',' token1];
			end
			dir=[]; plane=[]; crop1=[]; crop2=[]; % Clear the variables
			thisok = true; % No problems in this token yet
			thisCorrString = token(1); % start the corrected string
			switch token(1)
				case {'Y','V','C'}, dir = 2; % Coronal (2nd dimension)
				case {'X','H','S'}, if token(1:2)=='SH'
						dir = 3; % Transaxial (3rd dimension)
						thisCorrString = [thisCorrString 'H'];
						token = token(2:end);
					else
						dir = 1; % Sagittal (1st dimension)
					end
				case {'Z','T'},	dir = 3; % Transaxial (3rd dimension)
				otherwise
					thisCorrString = '';
					thisok = false;
			end
			if thisok & length(token)>1
				token = token(2:end);
			else
				thisok = false;
				thisCorrString = '';
			end

			% Next - choose the plane number in the given direction
			if thisok
				i = 1;
				while i<=length(token) & ((token(i)>='0' & token(i)<='9') | token(i)==':'), i=i+1; end
				plane = str2num(token(1:min(length(token),i-1))); % Have the plane number and direction
				if isempty(plane)
					thisok = false;
					thisCorrString = '';
				else
					if plane(1)<1, plane=1:max(plane(end),1); thisok = false; end
					if ~nodim & plane(end)>dim(dir), plane=min(plane(1),dim(dir)):dim(dir); thisok = false; end
					if length(plane)==1
						thisCorrString=[thisCorrString num2str(plane)];
					else
						thisCorrString=[thisCorrString num2str(plane(1)) '-' num2str(plane(end))];
					end
					if i>length(token) % Does token contain more data?
						token = '';
					else
						token = token(i:end);
					end
				end
				% Look for cropping
				if ~isempty(token)
					if token(1) ~= '('
						thisok = false;
					else
						if token(end)==')', token=token(2:end-1); end % ) is not critical
						[token, trem] = strtok(token,',');
						crop1 = str2num(token);
						[token, trem] = strtok(trem,',');
						crop2 = str2num(token);
						if isempty(crop1) | isempty(crop2) | ~isempty(trem)
							thisok = false;
						else
							switch dir
								case 1
									if crop1(1)<1, crop1=1:crop1(end); thisok = false; end
									if crop2(1)<1, crop2=1:crop2(end); thisok = false; end
									if ~nodim & crop1(end)>dim(2), crop1=crop1(1):dim(2); thisok = false; end
									if ~nodim & crop2(end)>dim(3), crop2=crop2(1):dim(3); thisok = false; end
								case 2
									if crop1(1)<1, crop1=1:crop1(end); thisok = false; end
									if crop2(1)<1, crop2=1:crop2(end); thisok = false; end
									if ~nodim & crop1(end)>dim(1), crop1=crop1(1):dim(1); thisok = false; end
									if ~nodim & crop2(end)>dim(3), crop2=crop2(1):dim(3); thisok = false; end
								case 3
									if crop1(1)<1, crop1=1:crop1(end); thisok = false; end
									if crop2(1)<1, crop2=1:crop2(end); thisok = false; end
									if ~nodim & crop1(end)>dim(1), crop1=crop1(1):dim(1); thisok = false; end
									if ~nodim & crop2(end)>dim(2), crop2=crop2(1):dim(2); thisok = false; end
							end % switch
							thisCorrString = [thisCorrString '(' num2str(crop1(1)) '-' num2str(crop1(end)) ...
								',' num2str(crop2(1)) '-' num2str(crop2(end)) ,')'];
						end % crop1 and crop2 not empty
					end % '(' at start of token
				else % token is empty
					if nodim
						msgId = sprintf('FADS:%s:wrongNumberOfInputParameters', mfilename);
						error(msgId, ['Cannot resolve uncropped plane if no dimensions are provided.']);
					end
					switch dir
						case 1, crop1=1:dim(2);  crop2=1:dim(3);
						case 2, crop1=1:dim(1);  crop2=1:dim(3);
						case 3, crop1=1:dim(1);  crop2=1:dim(2);
					end
				end % token empty
			end

			if thisok  % Token was valid
				switch dir
					case 1, smask = {plane crop1 crop2};
					case 2, smask = {crop1 plane crop2};
					case 3, smask = {crop1 crop2 plane};
				end
				mask=[mask; smask];
			else % Not a valid token - the entry is not valid and the corrected string will be returned
				nonvalid = true;
			end
			if ~isempty(thisCorrString)
				corrString = [corrString ',' thisCorrString];
			end
			[token, fString] = strtok(fString,',');
		end % of while loop
	end % if slices specified
end % is fString empty

if nonvalid % Was not a valid string then return the corrected string
	if ~isempty(corrString)  % Remove the first character as it is a comma
		corrString = corrString(2:end);
	end
	mask = corrString;
end  % Otherwise the mask is returned