function [allP_an,allC_an,avg_C_conds,mRs_trimmed,all_corr_an] = get_pop_vector_corr(sel_out,conditionsAndRasterTypes,ncols)

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
        [allP_an{an,ii},corr,cell_nums{an,ii}] = findPopulationVectorPlot(mRs_trimmed{an,ii},[]);
        corr = fillmissing(corr,'linear',2,'EndValues','nearest');
        corr = fillmissing(corr,'linear',1,'EndValues','nearest');
        allC_an{an,ii} = corr;
        avg_C_conds{ii}(:,:,an) = allC_an{an,ii};
        
    end
end
an = 1; 
if length(cell_nums{an,1}) == length(cell_nums{an,2}) & length(cell_nums{an,3}) == length(cell_nums{an,4})
% combs = nchoosek(1:4,2);
    for an = 1:size(sel_out.sz,1)
        all_corr = [];
        for ii = 1:4
            for jj = 1:4
                set1 = mRs_trimmed{an,ii}; set2 = mRs_trimmed{an,jj};
                [pv1,~,cellnums] = findPopulationVectorPlot(set1,[]);
                [pv2,~,~] = findPopulationVectorPlot(set2,[],cellnums);
                corr = find_pv_corr(pv1,pv2);
                corr = fillmissing(corr,'linear',2,'EndValues','nearest');
                corr = fillmissing(corr,'linear',1,'EndValues','nearest');
                all_corr{ii,jj} = corr;
    %         figure(100);clf;subplot 121;imagesc(pv1); subplot 122;imagesc(pvc1);%set(gca,'Ydir','Normal');
    %         figure(101);clf;subplot 121;imagesc(pv2); subplot 122;imagesc(pvc2);%set(gca,'Ydir','Normal');
    %         figure(103);clf;imagesc(corr);colorbar;%set(gca,'Ydir','Normal');colorbar
    %         n = 0;
            end
        end
        all_corr_an{an} = all_corr;
    end
else
    for an = 1:size(sel_out.sz,1)
        all_corr_an{an} = NaN;
    end
end

    
n = 0;