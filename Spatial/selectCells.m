function selCells = selectCells(data,mData,selText,option)

% selCells = selectCells_1(data,mData,selText);
% return;
for ii = 1:length(data)
    c{ii} = selectCells_1(data,mData,ii);
end

for jj = 1:length(data)
    sci = jj;
    selCells = c{sci};
    ci = 1:length(data); ci(sci) = [];
    for ii = 1:length(ci)
        selCells = setdiff(selCells,c{ci(ii)});
    end
    only_c{jj} = selCells;
end

if strcmp(selText,'Context')
    if ischar(option)
        number = str2num(option(1));
    else
        number = option;
    end
    if length(option) == 2
        cmdTxt = sprintf('selCells = only_c{number};');
    else
        cmdTxt = sprintf('selCells = c{number};');
    end
    eval(cmdTxt);
    return;
end


if strcmp(selText,'Common')
    list = option;
    selCells = c{list(1)};
    for ii = 2:length(list)
        selCells = intersect(selCells,c{list(ii)});
    end
    return;
end
n = 0;


% 
% 
% for ii = 1:4
%     c{ii} = selectCells(data,mData,ii);
% end
% for jj = 1:4
%     sci = jj;
%     selCells = c{sci};
%     ci = 1:4; ci(sci) = [];
%     for ii = 1:length(ci)
%         selCells = setdiff(selCells,c{ci(ii)});
%     end
%     only_c{jj} = selCells;
% end
% purePIs = selectCells(data,mData,[1]);
% selCells = purePIs;
% originalSC = purePIs;
% for ii = 2:4
%     selCells = intersect(selCells,data{ii}.sel);
% end
% persistedC = selCells;
% cCells = setxor(originalSC,persistedC);
% selCells = setxor(allCells,purePIs);
% selCells = intersect(selCells,data{3}.sel);
% 
% 
% ii = 2;
% selCells = c{ii};
% % temp = [1 2 3 4];
% % temp(ii) = [];
% % for jj = 1:length(temp)
% %     selCells = setdiff(selCells,c{temp(jj)});
% % end
% % selCells = persistedC;
% selCells = intersect(c{2},c{3});


function selCells = selectCells_1(data,mData,list)
totalCells = sum(mData.cellsN);
selCells = 1:totalCells;

if ~exist('list','var')
    return;
end

if ~exist('sith','var')
    sith = 3;
end
% if ~exist('sith','var')
%     sith = mData.sith;
% else
%     temp = nan(1:length(list));
%     temp(list) = sith;
%     sith = temp;
% end


for ii = 1:length(list)
%     if gl(ii)
%         theseCells = find(data{list(ii)}.SI > sith(list(ii)));
%     else
%         theseCells = find(data{list(ii)}.SI < sith(list(ii)));
%     end
    theseCells = find(data{list(ii)}.SI > sith);
    theseCells = intersect(theseCells,data{list(ii)}.sel);
    if strcmp(data{list(ii)}.name,'motion onsets')
        thisData = data{list(ii)};
        ctd = [];
        for jj = 1:length(theseCells)
            thisRaster = thisData.rasters(:,:,theseCells(jj));
            cols = size(thisRaster,2);
            firstR = mean2(thisRaster(:,1:(cols/2)));
            secondR = mean2(thisRaster(:,(cols/2):end));
            if secondR < firstR
                ctd = [ctd jj];
            end
        end
        theseCells(ctd) = [];
    end
    if strcmp(data{list(ii)}.name,'motion offsets')
        thisData = data{list(ii)};
        ctd = [];
        for jj = 1:length(theseCells)
            thisRaster = thisData.rasters(:,:,theseCells(jj));
            cols = size(thisRaster,2);
            firstR = mean2(thisRaster(:,1:(cols/2)));
            secondR = mean2(thisRaster(:,(cols/2):end));
            if secondR < firstR
                ctd = [ctd jj];
            end
        end
        theseCells(ctd) = [];
    end
    selCells = intersect(selCells,theseCells);
end


% p = inputParser;
% default_colors = distinguishable_colors(20);
% default_ySpacingFactor = 10;
% addRequired(p,'distD',@iscell);
% addOptional(p,'incr',0,@isnumeric);
% addOptional(p,'colors',default_colors,@iscell);
% addOptional(p,'maxY',100,@isnumeric);
% addOptional(p,'cumPos',0,@isnumeric);
% addOptional(p,'legend',0,@iscell);
% parse(p,distD,varargin{:});
% 
% cols = p.Results.colors;
% maxY = p.Results.maxY;
% cumPos = p.Results.cumPos;
% temp = p.Results.legend;
% legs = temp(1:(end-1));
% specs = temp{end};