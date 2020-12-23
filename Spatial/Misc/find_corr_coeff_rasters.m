function out = find_corr_coeff_rasters(pMs_C,paramMs_C,selAnimals_C,ei_C,conditionsAndRasterTypes,selC,cpMs,trials)

if ~exist('trials','var')
    trials = 1:10;
end

all_conds = []; all_rts = [];
for rr = 1:size(pMs_C,1)
    for cc = 1:size(pMs_C,2)
        tcond = conditionsAndRasterTypes(rr,cc);
        nds = dec2base(tcond,10) - '0';
        varNames{rr,cc} = sprintf('C%dR%d',nds(1),nds(2));
        all_conds = [all_conds nds(1)]; all_rts = [all_rts nds(2)];
        xticklabels{cc,rr} = sprintf('%s-%s',paramMs_C.stimMarkers{nds(2)},paramMs_C.rasterTypes{nds(2)}(1));
        for an = 1:length(selAnimals_C)
            tei = ei_C(an); conditionNumber = nds(1); rasterType = paramMs_C.rasterTypes{nds(2)}; 
            stimMarker = paramMs_C.stimMarkers{nds(2)}; maxDistTime = [140 Inf];%paramMs_C.maxDistTime;
%             selCells = pMs_C{conditionNumber}.cellSel{an};
            selCells = cpMs.cellSel{an};
            cns = paramMs_C.all_cns{an};
            [temp_rasters ~] = getParamValues('rasters',tei,selC.plane_number,conditionNumber,stimMarker,rasterType,cns(selCells,2:3),maxDistTime);
            this_mean_rasters = squeeze(nanmean(temp_rasters(trials,:,:),1))';
            mean_rasters_C{an,cc} = this_mean_rasters;
        end
    end
end

combs = [1 2;
        2 3;
        3 4];
all_ccs = cell(1);
gAllVals = [];
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
out.all_css_C = all_ccs_C;
out.gAllVals = gAllVals;
out.mean_over_cells = mean_over_cells_C;
out.mean_rasters_C = mean_rasters_C;