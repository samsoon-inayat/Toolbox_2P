function out = get_mean_rasters(pMs_C,paramMs_C,selAnimals_C,ei_C,conditionsAndRasterTypes,selC,cpMs,trials)

if ~exist('trials','var')
    trials = 1:10;
end
[szr,szc] = size(pMs_C);
if szc == 1
    pMs_C = pMs_C';
    conditionsAndRasterTypes = conditionsAndRasterTypes';
end
all_conds = []; all_rts = [];
for rr = 1:size(pMs_C,1)
    for cc = 1:size(pMs_C,2)
        tcond = conditionsAndRasterTypes(rr,cc);
        nds = dec2base(tcond,10) - '0';
        varNames{rr,cc} = sprintf('C%dR%d',nds(1),nds(2));
        all_conds = [all_conds nds(1)]; all_rts = [all_rts nds(2)];
        xticklabels{cc,rr} = sprintf('%s-%s',paramMs_C.stimMarkers{nds(2)},paramMs_C.rasterTypes{nds(2)}(1));
        for ani = 1:length(selAnimals_C)
            an = selAnimals_C(ani);
            tei = ei_C(an); conditionNumber = nds(1); rasterType = paramMs_C.rasterTypes{nds(2)}; 
            stimMarker = paramMs_C.stimMarkers{nds(2)}; maxDistTime = [150 Inf];%paramMs_C.maxDistTime;
            if ischar(cpMs)
                selCells = pMs_C{conditionNumber}.cellSel{an};
            end
            if isstruct(cpMs)
                selCells = cpMs.cellSel{an};
            end
            cns = paramMs_C.all_cns{an};
            [temp_rasters ~] = getParamValues('rasters',tei,selC.plane_number,conditionNumber,stimMarker,rasterType,cns(selCells,2:3),maxDistTime);
            this_mean_rasters = squeeze(nanmean(temp_rasters(trials,:,:),1))';
            this_mean_rasters = fillmissing(this_mean_rasters,'linear',2,'EndValues','none');
            [nrr,ncc] = find(isnan(this_mean_rasters));
            mcc = min(ncc);
            if ~isempty(mcc)
            if mcc == 1
                incc = ncc == 1;
                nrrs = nrr(incc);
                check_mat = this_mean_rasters(nrrs,:);
                if sum(isnan(check_mat(:))) == numel(check_mat)
                    this_mean_rasters(nrrs,:) = [];
                    [nrr,ncc] = find(isnan(this_mean_rasters));
                    mcc = min(ncc);
                else
                    error;
                end
            end
            this_mean_rasters = this_mean_rasters(:,1:mcc-1);
            end
            if isempty(this_mean_rasters)
                n = 0;
            end
%             this_mean_rasters = fillmissing(this_mean_rasters,'linear',2,'EndValues','nearest');%inpaint_nans(thisR(rnan(jj),:),4);
            mean_rasters_C{ani,cc} = this_mean_rasters;
            [temp ~] = getParamValues('',tei,selC.plane_number,conditionNumber,stimMarker,rasterType,cns(selCells,2:3),maxDistTime);
            if strcmp(rasterType,'dist')
                tempxs =  1.5:3:1000;
                xs{ani,cc} = tempxs(1:size(this_mean_rasters,2));
                sz(ani,cc) = length(xs{ani,cc});
                xs{ani,cc} = tempxs;
            end
            if strcmp(rasterType,'time')
                tempxs =  0.1:0.2:100;
                xs{ani,cc} = tempxs(1:size(this_mean_rasters,2));
                sz(ani,cc) = length(xs{ani,cc});
                xs{ani,cc} = tempxs;
            end
            if size(this_mean_rasters,2) == 0
                n = 0;
            end
        end
    end
end
out.mean_rasters = mean_rasters_C;
out.xs = xs;
out.sz = sz;

