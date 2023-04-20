function FR_C = get_spike_rate_ext(ei_C,thr,cellList)

for ii = 1:length(ei_C)
    ei = ei_C{ii};
    for jj = 1:4
        aironsets = ei.plane{1}.contexts(jj).markers.air_onsets;
        airoffsets = ei.plane{1}.contexts(jj).markers.air_offsets;
        FR_T = [];
        for tt = 1:10
            ton = aironsets(tt); toff = airoffsets(tt);
            spSigAll = [];
            all_spSigAll = [];
            for pp = 1:length(ei.plane)
                tspSigAll = [];
                this_cell_list = logical(ei.plane{pp}.tP.iscell(:,1));
%                 tspSigAll = ei.plane{pp}.tP.deconv.spSigAll(this_cell_list,:);
                tspSigAll = get_calcium_data_raw(ei,pp); tspSigAll = tspSigAll(this_cell_list,:);
                mask = tspSigAll > thr;
                tspSigAll(~mask) = NaN;
                spSigAll{pp} = tspSigAll;
                frames_f = ei.plane{pp}.b.frames_f;
                inds = find(frames_f > ton & frames_f < toff);
                all_spSigAll = [all_spSigAll;nanmean(tspSigAll(:,inds),2)];
            end
            FR_T(:,tt) = all_spSigAll;
        end
        FR_C{ii,jj} = nanmean(FR_T(cellList{ii,jj},:),2);
    end
    for jj = 5:8
        aironsets = ei.plane{1}.contexts(jj-4).markers.airI_onsets;
        airoffsets = ei.plane{1}.contexts(jj-4).markers.airI_offsets;
        FR_T = [];
        for tt = 1:10
            ton = aironsets(tt); toff = airoffsets(tt);
            spSigAll = [];
            all_spSigAll = [];
            for pp = 1:length(ei.plane)
                tspSigAll = [];
                this_cell_list = logical(ei.plane{pp}.tP.iscell(:,1));
%                 tspSigAll = ei.plane{pp}.tP.deconv.spSigAll(this_cell_list,:);
                tspSigAll = get_calcium_data_raw(ei,pp); tspSigAll = tspSigAll(this_cell_list,:);
                mask = tspSigAll > thr;
                tspSigAll(~mask) = NaN;
                spSigAll{pp} = tspSigAll;
                frames_f = ei.plane{pp}.b.frames_f;
                inds = find(frames_f > ton & frames_f < toff);
                all_spSigAll = [all_spSigAll;nanmean(tspSigAll(:,inds),2)];
            end
            FR_T(:,tt) = all_spSigAll;
        end
        FR_C{ii,jj} = nanmean(FR_T(cellList{ii,jj},:),2);
    end
end
