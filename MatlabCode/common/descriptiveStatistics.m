function out = descriptiveStatistics (data,varargin)

p = inputParser;
addRequired(p,'data',@isnumeric);
addOptional(p,'dimension',1,@isnumeric);
addOptional(p,'decimal_places',3,@isnumeric);
% addOptional(p,'max',inf,@isnumeric);
% addOptional(p,'colors',default_colors,@iscell);
% addOptional(p,'maxY',100,@isnumeric);
% addOptional(p,'cumPos',0,@isnumeric);
% addOptional(p,'barGraph',{},@iscell);
% addOptional(p,'legend',{},@iscell);
% addOptional(p,'BaseValue',0.2,@isnumeric);
addOptional(p,'do_mean','Yes');
parse(p,data,varargin{:});

dimension = p.Results.dimension;
decimal_places = p.Results.decimal_places;

[avg,se,sd,med] = findMeanAndStandardError(data,dimension);

out.avg = avg;
out.med = med;
out.sem = se;
out.sd = sd;
minD = min(data,[],dimension);
maxD = max(data,[],dimension);
out.min = minD;
out.max = maxD;

if isvector(data)
    txt = sprintf('(average: %%.%df; median: %%.%df; standard deviation: %%.%df; range: %%.%df,%%.%df)',decimal_places,decimal_places,decimal_places,decimal_places,decimal_places);
% %     if decimal_places == 0
% %         avg = round(avg);sd = round(sd);med = round(med);minD = round(minD); maxD = round(maxD);
% %     end
    cmdTxt = sprintf('out.txt = sprintf(''%s'',avg,sd,med,minD,maxD);',txt);
    eval(cmdTxt)
end