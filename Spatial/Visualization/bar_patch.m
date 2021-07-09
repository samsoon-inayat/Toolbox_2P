function hp = bar_patch(x,m,varargin)
if nargin == 1
    data = x;
    x = 1:size(data,2);
    [m,s] = findMeanAndStandardError(data);
    barwidth = 0.7;
    for ii = 1:length(x)
        hp(ii) = bar_patch(x(ii),m(ii),s(ii),barwidth);
    end
    return;
end
if nargin == 2
    barwidth = 0.7;
    s = [];
end
if nargin == 3
    s = varargin{1};
    barwidth = 0.7;
end
if nargin == 4
    s = varargin{1};
    barwidth = varargin{2};
end
% barwidth = 0.7;
X = [(x-barwidth/2) (x+barwidth/2) (x+barwidth/2) (x-barwidth/2)];
Y = [(0) (0) (m) (m)];
hp = patch(X,Y,'c');