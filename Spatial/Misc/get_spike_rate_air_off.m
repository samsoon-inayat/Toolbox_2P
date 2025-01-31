function out = get_spike_rate(ei_C,thr,cellList)
allVals = []; allVals_motion = []; allVals_rest = []; allVals_th = [];
allVals_an = [];
for ii = 1:length(ei_C)
    ei = ei_C{ii};
    Fso = ei.thorExp.frameRate;
%     speed = ei.b.air_puff_raw;
%     inds_motion = find(speed > 0.5);
%     inds_rest = find(speed < 0.5);
    
    speed = ei.b.fSpeed;
    air_offsets = ei.b.air_puff_f(11:39);
    air_onsets = ei.b.air_puff_r(12:40);
    
   
%     inds_motion = find(speed > 0);
    inds_motion = find(speed ~= 0);
    inds_rest = find(speed == 0);

    inds_motion2 = [];     inds_rest2 = [];
    for aoii = 1:length(air_onsets)
        inds_motion2 = [inds_motion2 intersect(inds_motion,air_offsets(aoii):air_onsets(aoii))];
        inds_rest2 = [inds_rest2 intersect(inds_rest,air_offsets(aoii):air_onsets(aoii))];
    end
    inds_motion = inds_motion2;     inds_rest = inds_rest2;

    onset_first_trial = ei.b.air_puff_r(1)-1000;
    offset_last_trial = ei.b.air_puff_f(end)+1000;
    inds_motion1 = inds_motion(find(inds_motion > onset_first_trial & inds_motion < offset_last_trial));
    inds_rest1 = inds_rest(find(inds_rest > onset_first_trial & inds_rest < offset_last_trial));
    percent_motion(ii) = 100*length(inds_motion1)/(offset_last_trial - onset_first_trial);
    percent_rest(ii) = 100*length(inds_rest1)/(offset_last_trial - onset_first_trial);
    
    spSigAll = [];
    all_spSigAll = [];
    for pp = 1:length(ei.plane)
        tspSigAll = [];
        this_cell_list = logical(ei.plane{pp}.tP.iscell(:,1));
        tspSigAll = ei.plane{pp}.tP.deconv.spSigAll(this_cell_list,:);
        % tspSigAll = fillmissing(tspSigAll, 'linear', 2);
        % Find rows where all elements are NaN and replace them with a constant (e.g., 0)
        % nan_rows = all(isnan(tspSigAll), 2);
        % tspSigAll(nan_rows, :) = 0;  % Replace entire NaN rows with 0
        tspSigAll_mean = nanmean(tspSigAll,2); tspSigAll_mean1 = repmat(tspSigAll_mean,1,size(tspSigAll,2));
        tspSigAll_std = nanstd(tspSigAll,[],2); tspSigAll_std1 = repmat(tspSigAll_std,1,size(tspSigAll,2));
        tspSigAll = (tspSigAll - tspSigAll_mean1)./tspSigAll_std1;
%         for cn = 1:length(this_cell_list)
%             tspSigAll(cn,:) = ei.plane{pp}.tP.deconv.spSigAll{this_cell_list(cn)}';
%         end
        % mask = tspSigAll > thr;
        % tspSigAll(~mask) = NaN;
        spSigAll{pp} = tspSigAll;
        all_spSigAll = [all_spSigAll;tspSigAll];
        frames_f = ei.plane{pp}.b.frames_f;
        [~,ind_frames_motion{pp},~] = intersect(frames_f,inds_motion);
        [~,ind_frames_rest{pp},~] = intersect(frames_f,inds_rest);
    end
    tm_sp_animal_motion = []; tm_sp_animal_rest =[]; a_ppmR = []; a_ppmM = [];
    for pp = 1:length(ei.plane)
        tspSigAll = spSigAll{pp};
        temp_frames = ind_frames_motion{pp}; temp_frames = temp_frames(temp_frames <= size(tspSigAll,2));
        tm_sp_animal_motion = [tm_sp_animal_motion;nanmean(tspSigAll(:,temp_frames),2)]; sig_motion = tspSigAll(:,temp_frames);
        temp_frames = ind_frames_rest{pp}; temp_frames = temp_frames(temp_frames <= size(tspSigAll,2));
        tm_sp_animal_rest = [tm_sp_animal_rest;nanmean(tspSigAll(:,temp_frames),2)]; sig_rest = tspSigAll(:,temp_frames);
        
        ppm = NaN(size(sig_rest,1),1);maxTime = size(sig_rest,2)/Fso;
        parfor cni = 1:size(sig_rest,1)
            X = sig_rest(cni,:);
            [pks,locs,w,p] = findpeaks(X,Fso);
            ppm(cni) = length(pks)/(maxTime/60);
        end
        a_ppmR = [a_ppmR;ppm];
        
        ppm = NaN(size(sig_motion,1),1);maxTime = size(sig_motion,2)/Fso;
        parfor cni = 1:size(sig_motion,1)
            X = sig_motion(cni,:);
            [pks,locs,w,p] = findpeaks(X,Fso);
            ppm(cni) = length(pks)/(maxTime/60);
        end
        a_ppmM = [a_ppmM;ppm];
    end
    m_sp_animal{ii} = nanmean(all_spSigAll,2);
    m_sp_animal_resp_cells{ii} = nanmean(all_spSigAll(cellList{ii},:),2);
    m_sp_animal_motion{ii} = tm_sp_animal_motion;
    m_sp_animal_rest{ii} = tm_sp_animal_rest ;
    m_tr_animal_rest{ii} = a_ppmR;
    m_tr_animal_motion{ii} = a_ppmM;

    allVals = [allVals;m_sp_animal{ii}];
    allVals_an{ii} = all_spSigAll;
    allVals_motion = [allVals_motion;m_sp_animal_motion{ii}];
    allVals_rest = [allVals_rest;m_sp_animal_rest{ii}];
    m_sp_animal_level(ii) = nanmean(m_sp_animal{ii});
    m_sp_animal_level_resp_cells(ii) = nanmean(m_sp_animal_resp_cells{ii});
    m_sp_animal_level_motion(ii) = nanmean(m_sp_animal_motion{ii});
    m_sp_animal_level_rest(ii) = nanmean(m_sp_animal_rest{ii});
%     thr = nanmean(spSigAll,2) + 3*nanstd(spSigAll,[],2);
%     for cn = 1:size(spSigAll,1)
%         m_sp_animal_th{ii}(cn,1) = nanmean(spSigAll(cn,spSigAll(cn,:) > thr(cn)));
%     end
%     m_sp_animal_level_th(ii) = nanmean(m_sp_animal_th{ii});
%     allVals_th = [allVals;m_sp_animal_th{ii}];
end
out.m_sp_animal = m_sp_animal;
out.allVals = allVals;
out.m_sp_animal_level = m_sp_animal_level;
out.m_sp_animal_level_resp_cells = m_sp_animal_level_resp_cells;
out.m_sp_animal_motion = m_sp_animal_motion;
out.allVals_motion = allVals_motion;
out.m_sp_animal_level_motion = m_sp_animal_level_motion;
out.m_sp_animal_rest = m_sp_animal_rest;
out.allVals_rest = allVals_rest;
out.m_sp_animal_level_rest = m_sp_animal_level_rest;
out.percent_motion = percent_motion;
out.allVals_an = allVals_an;
out.m_tr_animal_rest = m_tr_animal_rest;
out.m_tr_animal_motion = m_tr_animal_motion;

% out.m_sp_animal_th = m_sp_animal_th;
% out.m_sp_animal_level_th = m_sp_animal_level_th;
% out.allVals_th = allVals_th;