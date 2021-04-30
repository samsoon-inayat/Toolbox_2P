function [popVs,cross_pos_corr,cross_cell_corr,cell_nums] = calc_pop_vector_corr(meanRs,ccs,cellNumsI)

if ~isempty(ccs)
    if ~iscell(ccs)
        for rr = 1:size(meanRs,1)
            for cc = 1:size(meanRs,2)
                thisMR = meanRs{rr,cc};
                meanRs{rr,cc} = thisMR(ccs,:);
            end
        end
    else
        for rr = 1:size(meanRs,1)
            for cc = 1:size(meanRs,2)
                thisMR = meanRs{rr,cc};
                meanRs{rr,cc} = thisMR(ccs{rr,cc},:);
            end
        end
    end
end

if ~isempty(cellNumsI)
    if length(cellNumsI) > 2
        cell_nums = cellNumsI;
    else
        [~,~,cell_nums] = findPopulationVectorPlot(meanRs{cellNumsI(1),cellNumsI(2)},[]);
    end
else
    cell_nums = [];
end

for rr = 1:size(meanRs,1)
    for cc = 1:size(meanRs,2)
        thisMR = meanRs{rr,cc};
        [tempP,corrV,CNs] = findPopulationVectorPlot(thisMR,[],cell_nums);
        corrV = fillmissing(corrV,'linear',2,'EndValues','nearest');
        corrV = fillmissing(corrV,'linear',1,'EndValues','nearest');
        thispv.popV = tempP;
        thispv.corrV = corrV;
        thispv.cell_nums = CNs;
        popVs{rr,cc} = thispv;
    end
end
n = 0;
if size(meanRs,2) == 1 || nargout == 1
    cross_pos_corr = [];
    cross_cell_corr = [];
    return;
end
for rr = 1:size(meanRs,1)
    this_cross_pos = [];
    this_cross_cell = [];
    for cc1 = 1:size(meanRs,2)
        set1 = meanRs{rr,cc1};
        a_corrV = NaN(size(set1,2),size(set1,2),size(meanRs,2));;
        a_corrCV = NaN(size(set1,1),size(set1,1),size(meanRs,2));;
        for cc2 = 1:size(meanRs,2)
            set2 = meanRs{rr,cc2};
            [pv1,~,cellnums] = findPopulationVectorPlot(set1,[]);
            [pv2,~,~] = findPopulationVectorPlot(set2,[],cellnums);
            [corrV,pV] = corr(pv1,pv2);
            corrV = fillmissing(corrV,'linear',2,'EndValues','nearest');
            corrV = fillmissing(corrV,'linear',1,'EndValues','nearest');
            a_corrV(:,:,cc2) = corrV;
%             this_cross_pos{cc1,cc2} = corrV;
            [corrCV,pCV] = corr(pv1',pv2');
            corrCV = fillmissing(corrCV,'linear',2,'EndValues','nearest');
            corrCV = fillmissing(corrCV,'linear',1,'EndValues','nearest');
            a_corrCV(:,:,cc2) = corrCV;
%             this_cross_cell{cc1,cc2} = corrCV;
        end
        for cc2 = 1:size(meanRs,2)
            this_cross_pos{cc1,cc2} = a_corrV(:,:,cc2);
            this_cross_cell{cc1,cc2} = a_corrCV(:,:,cc2);
        end
    end
    cross_pos_corr{rr,1} = this_cross_pos;
    cross_cell_corr{rr,1} = this_cross_cell;
end

