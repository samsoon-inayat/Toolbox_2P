function FR_C = get_spike_rate_ext_dist_dur(ei_C,thr,cellList,selCon)

for ii = 1:length(ei_C)
    ei = ei_C{ii};
    for jj = 1:2:10
        jjj = selCon(jj);
        if ii == 1
        disp(ei.plane{1}.contexts(jjj).name);
        end
        if strcmp(ei.plane{1}.contexts(jjj).name,'Air - Brake')
            aironsets = ei.plane{1}.contexts(jjj).markers.airR_onsets;
            airoffsets = ei.plane{1}.contexts(jjj).markers.airR_offsets;
        else
            aironsets = ei.plane{1}.contexts(jjj).markers.air_onsets;
            airoffsets = ei.plane{1}.contexts(jjj).markers.air_offsets;
        end
        FR_T = [];
        for tt = 1:10
            ton = aironsets(tt); toff = airoffsets(tt);
            spSigAll = [];
            all_spSigAll = [];
            for pp = 1:length(ei.plane)
                tspSigAll = [];
                this_cell_list = logical(ei.plane{pp}.tP.iscell(:,1));
                tspSigAll = ei.plane{pp}.tP.deconv.spSigAll(this_cell_list,:);
                mask = tspSigAll > thr;
                tspSigAll(~mask) = NaN;
                spSigAll{pp} = tspSigAll;
                frames_f = ei.plane{pp}.b.frames_f;
                inds = find(frames_f > ton & frames_f < toff);
                all_spSigAll = [all_spSigAll;max(tspSigAll(:,inds),[],2)];
            end
            FR_T(:,tt) = all_spSigAll;
        end
        FR_C{ii,jj} = max(FR_T(cellList{ii,jj},:),[],2);
    end
    for jj = 2:2:10
        jjj = selCon(jj);
        if ii == 1
        disp(ei.plane{1}.contexts(jjj).name);
        end
        if strcmp(ei.plane{1}.contexts(jjj).name,'Air - Brake')
            aironsets = ei.plane{1}.contexts(jjj).markers.airIR_onsets;
            airoffsets = ei.plane{1}.contexts(jjj).markers.airIR_offsets;
        else
            aironsets = ei.plane{1}.contexts(jjj).markers.airI_onsets;
            airoffsets = ei.plane{1}.contexts(jjj).markers.airI_offsets;
        end
        FR_T = [];
        for tt = 1:10
            ton = aironsets(tt); toff = airoffsets(tt);
            spSigAll = [];
            all_spSigAll = [];
            for pp = 1:length(ei.plane)
                tspSigAll = [];
                this_cell_list = logical(ei.plane{pp}.tP.iscell(:,1));
                tspSigAll = ei.plane{pp}.tP.deconv.spSigAll(this_cell_list,:);
                mask = tspSigAll > thr;
                tspSigAll(~mask) = NaN;
                spSigAll{pp} = tspSigAll;
                frames_f = ei.plane{pp}.b.frames_f;
                inds = find(frames_f > ton & frames_f < toff);
                all_spSigAll = [all_spSigAll;max(tspSigAll(:,inds),[],2)];
            end
            FR_T(:,tt) = all_spSigAll;
        end
        FR_C{ii,jj} = max(FR_T(cellList{ii,jj},:),[],2);
    end
end
