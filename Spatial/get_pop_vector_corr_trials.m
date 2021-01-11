function out = get_pop_vector_corr_trials(sel_out)

n = 0;
n_an = size(sel_out{1}.mean_rasters,1);
n_cn = size(sel_out{1}.mean_rasters,2);
for rr = 1:length(sel_out)
    for cc = 1:length(sel_out)
        for an = 1:n_an
            for cn = 1:n_cn
                set1 = sel_out{rr}.mean_rasters_T{an,cn}; set2 = sel_out{cc}.mean_rasters_T{an,cn};
                [pv1,~,cellnums] = findPopulationVectorPlot(set1,[]);
                [pv2,~,~] = findPopulationVectorPlot(set2,[],cellnums);
                [corrV,pV] = corr(pv1,pv2);
                corrV = fillmissing(corrV,'linear',2,'EndValues','nearest');
                corrV = fillmissing(corrV,'linear',1,'EndValues','nearest');
                all_corr{an,cn} = corrV;
                [corrCV,pCV] = corr(pv1',pv2');
                corrCV = fillmissing(corrCV,'linear',2,'EndValues','nearest');
                corrCV = fillmissing(corrCV,'linear',1,'EndValues','nearest');
                all_corr_cell{an,cn} = corrCV;
            end
        end
        all_corr_trials{rr,cc} = all_corr;
        all_corr_cell_trials{rr,cc} = all_corr_cell;
    end
end

out.pos_corr = all_corr_trials;
out.cell_corr = all_corr_cell_trials;
disp('done pos cell corr');