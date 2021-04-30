function [rasters,status,diffbs] = correct_for_size(rasters)

status = 0;
binsizes = [];
for rr = 1:size(rasters,1)
    for cc = 1:size(rasters,2)
        binsizes(rr,cc) = size(rasters{rr,cc},2);
    end
end
diffbs = diff(binsizes,[],2);
if sum(diffbs) ~= 0
    status = 1;
    [rows,cols] = find(diffbs ~=0);
    urows = unique(rows);
    for ii = 1:length(urows)
        msz = min(binsizes(urows(ii),:));
        for cc = 1:size(rasters,2)
            rasters{urows(ii),cc} = rasters{urows(ii),cc}(:,1:msz,:);
        end
    end
end
n = 0;
