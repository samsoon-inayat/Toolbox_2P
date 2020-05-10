function pt = getPositionTuning(ei,type,N)

if iscell(ei)
    for ii = 1:length(ei)
        if ~exist('N','var')
            temp = getPositionTuning(ei{ii},type);
        else
            temp = getPositionTuning(ei{ii},type,N);
        end
        if ii == 1
            pt = temp;
        else
            pt = [pt;temp];
        end
    end
    return;
end

if ~exist('N','var')
    NN = 0;
else
    NN = 1;
end

if ~exist('type','var')
    type = 'sp';
end

if strcmp(lower(type),'sp')
    rasters = ei.rasters.rasters;
end
if strcmp(lower(type),'ca')
    rasters = ei.rasters.caRasters;
end
w = gausswin(3);
for ii = 1:size(rasters,3)
    thisRaster = rasters(:,:,ii);
    if NN
        thisRaster = thisRaster/max(thisRaster(:));
    end
    temp = nanmean(thisRaster);
%     pt(ii,:) = filter(w,1,temp);
    pt(ii,:) = temp;
%     figure(1000);clf;
%     plot(pt(ii,:));
%     n = 0;
end
