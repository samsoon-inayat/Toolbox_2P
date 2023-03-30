function FR_C = get_spike_rate_ext_transients(ei_C,thr,cellList)

for ii = 1:length(ei_C)
    ei = ei_C{ii};
    Fso = ei.thorExp.frameRate;
    for jj = 1:4
        aironsets = ei.plane{1}.contexts(jj).markers.air_onsets;
        airoffsets = ei.plane{1}.contexts(jj).markers.air_offsets;
        FR_T = [];
        for tt = 1:10
            ton = aironsets(tt); toff = airoffsets(tt);
            spSigAll = [];
            a_ppm = [];
            for pp = 1:length(ei.plane)
                tspSigAll = [];
                this_cell_list = logical(ei.plane{pp}.tP.iscell(:,1));
                tspSigAll = ei.plane{pp}.tP.deconv.spSigAll(this_cell_list,:);
                mask = tspSigAll > thr;
                tspSigAll(~mask) = NaN;
                spSigAll{pp} = tspSigAll;
                frames_f = ei.plane{pp}.b.frames_f;
                inds = find(frames_f > ton & frames_f < toff);
                all_spSigAll = tspSigAll(:,inds);
                ppm = NaN(size(all_spSigAll,1),1);maxTime = size(all_spSigAll,2)/Fso;
                parfor cni = 1:size(all_spSigAll,1)
                    X = all_spSigAll(cni,:);
                    [pks,locs,w,p] = findpeaks(X,Fso);
                    ppm(cni) = length(pks)/(maxTime/60);
                end
                a_ppm = [a_ppm;ppm];
            end
            FR_T(:,tt) = a_ppm;
        end
        FR_C{ii,jj} = mean(FR_T(cellList{ii,jj},:),2);
%         err
    end
    for jj = 5:8
        aironsets = ei.plane{1}.contexts(jj-4).markers.airI_onsets;
        airoffsets = ei.plane{1}.contexts(jj-4).markers.airI_offsets;
        FR_T = [];
        for tt = 1:10
            ton = aironsets(tt); toff = airoffsets(tt);
            spSigAll = [];
            a_ppm = [];
            for pp = 1:length(ei.plane)
                tspSigAll = [];
                this_cell_list = logical(ei.plane{pp}.tP.iscell(:,1));
                tspSigAll = ei.plane{pp}.tP.deconv.spSigAll(this_cell_list,:);
                mask = tspSigAll > thr;
                tspSigAll(~mask) = NaN;
                spSigAll{pp} = tspSigAll;
                frames_f = ei.plane{pp}.b.frames_f;
                inds = find(frames_f > ton & frames_f < toff);
                all_spSigAll = tspSigAll(:,inds);
                ppm = NaN(size(all_spSigAll,1),1);maxTime = size(all_spSigAll,2)/Fso;
                parfor cni = 1:size(all_spSigAll,1)
                    X = all_spSigAll(cni,:);
                    [pks,locs,w,p] = findpeaks(X,Fso);
                    ppm(cni) = length(pks)/(maxTime/60);
                end
                a_ppm = [a_ppm;ppm];
            end
            FR_T(:,tt) = a_ppm;
        end
        FR_C{ii,jj} = mean(FR_T(cellList{ii,jj},:),2);
    end
end
