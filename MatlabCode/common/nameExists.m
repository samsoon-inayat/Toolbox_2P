function existence = nameExists (varargin)

if nargin == 1
    nameWithPath = varargin{1};
end
if nargin == 2
    nameWithPath = makeName(varargin{1},varargin{2});
end
E = dir(nameWithPath);
existence = length(E);