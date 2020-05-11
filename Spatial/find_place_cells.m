function cellList = find_place_cells(S,varargin)

p = inputParser;
addRequired(p,'S',@isstruct);
default_width_threshold = 30; addOptional(p,'width_threshold',default_width_threshold,@isnumeric);
default_RSquare_threshold = 0.7; addOptional(p,'RSquare_threshold',default_width_threshold,@isnumeric);
default_adjRSquare_threshold = 0.7; addOptional(p,'adjRSquare_threshold',default_width_threshold,@isnumeric);
parse(p,S,varargin{:});
width_threshold = p.Results.width_threshold;
adjRSquare_threshold = p.Results.adjRSquare_threshold;
RSquare_threshold = p.Results.RSquare_threshold;
f = S.f;
gof = S.gof;
output = S.output;
cellList = [];
for ii = 1:length(f)
    if gof{ii}.adjrsquare > adjRSquare_threshold
        cellList = [cellList ii];
    end
end