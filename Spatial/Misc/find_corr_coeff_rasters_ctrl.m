function out = find_corr_coeff_rasters(in)

out = in;

combs = [1 2;
        2 3;
        3 4];
for an = 1:length(selAnimals_C)
    ccs = [];
    for ii = 1:size(combs,1)
        cond1 = combs(ii,1); cond2 = combs(ii,2);
        set1 = mean_rasters_C{an,cond1}; set2 = mean_rasters_C{an,cond2};
        cc = [];
        for cn = 1:size(set1,1)
            cell_sig1 = set1(cn,:); cell_sig1 = fillmissing(cell_sig1,'linear',2,'EndValues','nearest');
            cell_sig2 = set2(cn,:); cell_sig2 = fillmissing(cell_sig2,'linear',2,'EndValues','nearest');
            temp = corrcoef(cell_sig1,cell_sig2);
            cc(cn) = temp(1,2);
        end
        ccs(ii,:) = cc;
        gAllVals = [gAllVals cc];
        
        [pv1,pvc1,cellnums] = findPopulationVectorPlot(set1,[]);
        [pv2,pvc2,~] = findPopulationVectorPlot(set2,[],cellnums);
        corr = find_pv_corr(pv1,pv2);
        figure(100);clf;imagesc(pv1);%set(gca,'Ydir','Normal');
        figure(101);clf;imagesc(pv2);%set(gca,'Ydir','Normal');
        figure(103);clf;imagesc(corr);colorbar;%set(gca,'Ydir','Normal');colorbar
        n = 0;
    end
    all_ccs_C(an) = {ccs};
    mean_over_cells_C(an,:) = nanmean(ccs,2)';
end
out.all_css = all_ccs_C;
out.gAllVals = gAllVals;
out.mean_over_cells = mean_over_cells_C;
out.mean_rasters = mean_rasters_C;
out.xs = xs;
