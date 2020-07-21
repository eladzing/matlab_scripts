function result = MAKE_CUM_PROFILE_FROM_CUBE(data)
% Return a cumulative sum profile for cartesian cubes of extensive parameters. 
% Same as cumsum(MAKE_PROFILE_FROM_CUBE(data))
%
% @param data    The cartesian data cube.
%
% @returns     The cumulative sum profile.
%
error(nargchk(1,1,nargin,'struct'));

result = cumsum(MAKE_PROFILE_FROM_CUBE(data));