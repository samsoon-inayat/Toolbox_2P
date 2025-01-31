function [FR_C,FR_C_motion,FR_C_rest] = get_spike_rate_ext_dist_dur(ei_C,thr,cellList,selCon)

for ii = 1:length(ei_C)
    ei = ei_C{ii};
    speed = ei.b.fSpeed;
    inds_motionM = find(speed ~= 0);
    inds_restM = find(speed == 0);
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
        FR_T = []; FR_T_motion = []; FR_T_rest = [];
        for tt = 1:10
            ton = aironsets(tt); toff = airoffsets(tt);
            spSigAll = [];
            all_spSigAll = [];            all_spSigAll_motion = [];            all_spSigAll_rest = [];
            for pp = 1:length(ei.plane)
                tspSigAll = [];
                this_cell_list = logical(ei.plane{pp}.tP.iscell(:,1));
                tspSigAll = ei.plane{pp}.tP.deconv.spSigAll(this_cell_list,:);
                tspSigAll_mean = nanmean(tspSigAll,2); tspSigAll_mean1 = repmat(tspSigAll_mean,1,size(tspSigAll,2));
                tspSigAll_std = nanstd(tspSigAll,[],2); tspSigAll_std1 = repmat(tspSigAll_std,1,size(tspSigAll,2));
                tspSigAll = (tspSigAll - tspSigAll_mean1)./tspSigAll_std1;
                % mask = tspSigAll > thr;
                % tspSigAll(~mask) = NaN;
                spSigAll{pp} = tspSigAll;
                frames_f = ei.plane{pp}.b.frames_f;
                frames_f_motion = intersect(frames_f,inds_motionM); frames_f_rest = intersect(frames_f,inds_restM);
                
                inds = find(frames_f > ton & frames_f < toff);  
                all_spSigAll = [all_spSigAll;nanmean(tspSigAll(:,inds),2)];

                if isempty(frames_f_motion)
                    all_spSigAll_motion = [all_spSigAll_motion;nan(size(tspSigAll,1),1)];
                else
                    inds_motion = find(frames_f_motion > ton & frames_f_motion < toff);
                    all_spSigAll_motion = [all_spSigAll_motion;nanmean(tspSigAll(:,inds_motion),2)];
                end

                if isempty(frames_f_rest)
                    all_spSigAll_rest = [all_spSigAll_rest;nan(size(tspSigAll,1),1)];
                else
                    inds_rest = find(frames_f_rest > ton & frames_f_rest < toff);
                    all_spSigAll_rest = [all_spSigAll_rest;nanmean(tspSigAll(:,inds_rest),2)];
                end

            end
            FR_T(:,tt) = all_spSigAll;
            FR_T_motion(:,tt) = all_spSigAll_motion;
            FR_T_rest(:,tt) = all_spSigAll_rest;
        end
        FR_C{ii,jj} = nanmean(FR_T(cellList{ii,jj},:),2);
        FR_C_motion{ii,jj} = nanmean(FR_T_motion(cellList{ii,jj},:),2);
        FR_C_rest{ii,jj} = nanmean(FR_T_rest(cellList{ii,jj},:),2);
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
        FR_T = []; FR_T_motion = []; FR_T_rest = [];
        for tt = 1:10
            ton = aironsets(tt); toff = airoffsets(tt);
            spSigAll = [];
            all_spSigAll = []; all_spSigAll_motion = []; all_spSigAll_rest = [];
            for pp = 1:length(ei.plane)
                tspSigAll = [];
                this_cell_list = logical(ei.plane{pp}.tP.iscell(:,1));
                tspSigAll = ei.plane{pp}.tP.deconv.spSigAll(this_cell_list,:);
                tspSigAll_mean = nanmean(tspSigAll,2); tspSigAll_mean1 = repmat(tspSigAll_mean,1,size(tspSigAll,2));
                tspSigAll_std = nanstd(tspSigAll,[],2); tspSigAll_std1 = repmat(tspSigAll_std,1,size(tspSigAll,2));
                tspSigAll = (tspSigAll - tspSigAll_mean1)./tspSigAll_std1;
                % mask = tspSigAll > thr;
                % tspSigAll(~mask) = NaN;
                spSigAll{pp} = tspSigAll;
                frames_f = ei.plane{pp}.b.frames_f;
                frames_f_motion = intersect(frames_f,inds_motionM); frames_f_rest = intersect(frames_f,inds_restM);

                inds = find(frames_f > ton & frames_f < toff);
                all_spSigAll = [all_spSigAll;nanmean(tspSigAll(:,inds),2)];
                
                if isempty(frames_f_motion)
                    all_spSigAll_motion = [all_spSigAll_motion;nan(size(tspSigAll,1),1)];
                else
                    inds_motion = find(frames_f_motion > ton & frames_f_motion < toff);
                    all_spSigAll_motion = [all_spSigAll_motion;nanmean(tspSigAll(:,inds_motion),2)];
                end

                if isempty(frames_f_rest)
                    all_spSigAll_rest = [all_spSigAll_rest;nan(size(tspSigAll,1),1)];
                else
                    inds_rest = find(frames_f_rest > ton & frames_f_rest < toff);
                    all_spSigAll_rest = [all_spSigAll_rest;nanmean(tspSigAll(:,inds_rest),2)];
                end

            end
            FR_T(:,tt) = all_spSigAll;
            FR_T_motion(:,tt) = all_spSigAll_motion;
            FR_T_rest(:,tt) = all_spSigAll_rest;
        end
        FR_C{ii,jj} = nanmean(FR_T(cellList{ii,jj},:),2);
        FR_C_motion{ii,jj} = nanmean(FR_T_motion(cellList{ii,jj},:),2);
        FR_C_rest{ii,jj} = nanmean(FR_T_rest(cellList{ii,jj},:),2);
    end
end


% out_C.FR = FR_C;
% out_C.sp = sp_C;
% out_C_motion.FR = FR_C_motion;
% out_C_motion.sp = sp_C_motion;
% out_C_rest.FR = FR_C_rest;
% out_C_rest.sp = sp_C_rest;
