%% Fast sub2ind
% avoids for-loop (up to 3times faster than MATLAB's ind2sub.m)
%
% same as sub2ind
% BUT two different function calls (function automatically decides) :
%
% (a1) ind = sub2indFast(size,x,y); % for 2D
% (a2) ind = sub2indFast(size,x,y,z); % for 3D
% (b)  ind = sub2indFast(size,sub); % where sub = NX3 (3D) or Nx2 (2D)
%
% added on 2010-10-19 by evadittrich
function ind = sub2indFast(siz,varargin)

m = siz(1);
n = siz(2);
mn = m*n;

if length(varargin)==1
	if size(varargin{1},2)==3
		x = varargin{1}(:,1);
		y = varargin{1}(:,2);
		z = varargin{1}(:,3);
	elseif size(varargin{1},2)==2
		x = varargin{1}(:,1);
		y = varargin{1}(:,2);
	end
elseif length(varargin)==2
	x = varargin{1};
	y = varargin{2};
elseif length(varargin)==3	
	x = varargin{1};
	y = varargin{2};
	z = varargin{3};
else
	error('sub2indFast:illegalarguments','Illegal number of arguments!');
end

switch numel(siz)
	%2D
	case 2 
		ind = ((y-1).*m)+x;
	%3D
	case 3 
		ind = ((y-1).*m)+x+(z-1)*mn;
end
