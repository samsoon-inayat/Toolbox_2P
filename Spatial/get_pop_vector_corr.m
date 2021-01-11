function [allP_an,allC_an,avg_C_conds,mRs_trimmed,all_corr_an,all_corr_cell_an] = get_pop_vector_corr(sel_out,conditionsAndRasterTypes,ncols,cellSel_C)

for an = 1:size(sel_out.sz,1)
    for ii = 1:length(conditionsAndRasterTypes)
        tcond = abs(conditionsAndRasterTypes(ii));
        Ndigits = dec2base(tcond,10) - '0';
        mRsi = sel_out.mean_rasters{an,ii};
        if size(mRsi,2) < ncols
            cncols = size(mRsi,2);
            mRsi(:,(cncols+1:ncols)) = nan(size(mRsi,1),length(cncols+1:ncols));
        end
        mRs_trimmed{an,ii} = mRsi(:,1:ncols);
        if isempty(mRs_trimmed{an,ii})
            n = 0;
        end
        if size(mRs_trimmed{an,ii},1) == 0
            allC_an{an,ii} = NaN;
            avg_C_conds{ii}(:,:,an) = allC_an{an,ii};
            continue;
        end
        [allP_an{an,ii},corrV,cell_nums{an,ii}] = findPopulationVectorPlot(mRs_trimmed{an,ii},[]);
        corrV = fillmissing(corrV,'linear',2,'EndValues','nearest');
        corrV = fillmissing(corrV,'linear',1,'EndValues','nearest');
        allC_an{an,ii} = corrV;
        avg_C_conds{ii}(:,:,an) = allC_an{an,ii};
    end
end
if ~ischar(cellSel_C)
% combs = nchoosek(1:4,2);
    for an = 1:size(sel_out.sz,1)
        all_corr = [];
        for ii = 1:size(sel_out.sz,2)
            for jj = 1:size(sel_out.sz,2)
                set1 = mRs_trimmed{an,ii}; set2 = mRs_trimmed{an,jj};
                [pv1,~,cellnums] = findPopulationVectorPlot(set1,[]);
                [pv2,~,~] = findPopulationVectorPlot(set2,[],cellnums);
                [corrV,pV] = corr(pv1,pv2);
                corrV = fillmissing(corrV,'linear',2,'EndValues','nearest');
                corrV = fillmissing(corrV,'linear',1,'EndValues','nearest');
                all_corr{ii,jj} = corrV;
                [corrCV,pCV] = corr(pv1',pv2');
                corrCV = fillmissing(corrCV,'linear',2,'EndValues','nearest');
                corrCV = fillmissing(corrCV,'linear',1,'EndValues','nearest');
                all_corr_cell{ii,jj} = corrCV;
    %         figure(100);clf;subplot 121;imagesc(pv1); subplot 122;imagesc(pvc1);%set(gca,'Ydir','Normal');
    %         figure(101);clf;subplot 121;imagesc(pv2); subplot 122;imagesc(pvc2);%set(gca,'Ydir','Normal');
    %         figure(103);clf;imagesc(corrV);colorbar;%set(gca,'Ydir','Normal');colorbar
    %         n = 0;
            end
        end
        all_corr_an{an} = all_corr;
        all_corr_cell_an{an} = all_corr_cell;
    end
else
    for an = 1:size(sel_out.sz,1)
        all_corr_an{an} = NaN;
        all_corr_cell_an{an} = NaN;
    end
end

    
n = 0;