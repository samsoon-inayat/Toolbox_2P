function out = find_shifts(allresp_trials,allpeakL_trials,respDT,si_names)
if ~exist('si_names','var')
    si_names = repmat({'-'},1,size(allresp_trials,2));
end
trial_th = 2;
% for all active cells in a trial, I want to see the peak locations and average shift of peak locations for animals
for an = 1:5
    for cn = 1:size(allresp_trials,2)
        pLs = allpeakL_trials{an,cn};
        resp = allresp_trials{an,cn};
        if ~isempty(respDT)
            if size(respDT,2) == 1
                resp = sum(resp,2)>trial_th & respDT{an};
            else
                resp = sum(resp,2)>trial_th & respDT{an,cn};
            end
        else
            resp = sum(resp,2)>trial_th;
        end
        pLs = pLs(resp,:);
        
        tr_seq = [];
        for trn = 1:10
            [~,inds] = sort(pLs(:,trn));
            [~,inds1] = sort(inds);
            
            nan_temp = double(isnan(pLs(:,trn)));
            nan_temp(nan_temp == 1) = NaN; nan_temp(nan_temp == 0) = 1;
            tr_seq(:,trn) = (inds1/nansum(nan_temp)) .* nan_temp;
        end
        all_pos_seq{an,cn} = tr_seq;
        
%         [~,inds] = sort(pLs(:,1));
    %     figure(1000);clf;imagesc(pLs(inds,:));colorbar;
%             figure(1000);clf;imagesc(pLs);colorbar;
        100*sum(~isnan(pLs))/size(pLs,1);
        prob_participation = [];
        % find distribution of pLs
        bins = 1:1:50;
        [bar1 xs] = hist(pLs(:),bins); bar1 = 100*bar1/sum(bar1);
        dists{an,cn} = bar1;
%         figure(1000);clf;plot(xs,bar1);
%         title(sprintf('Animal %d,%s',an,si_names{cn}));
%         pause(0.3);
        mshift = []; 
        for ii = 1:size(pLs,1)
            tcell = pLs(ii,:);
            prob_participation(ii,1) = sum(~isnan(tcell))/length(tcell);
            ntcell = tcell(~isnan(tcell));
            mshift(ii,1) = mean(diff(ntcell));
        end
        mshifts(an,cn) = mean(mshift);
        mshiftsB(an,cn) = mean(mshift(mshift<0));
        mshiftsF(an,cn) = mean(mshift(mshift>0));
        fshiftB(an,cn) = sum(mshift<0)/length(mshift);
        fshiftF(an,cn) = sum(mshift>0)/length(mshift);
        fshiftZ(an,cn) = sum(mshift==0)/length(mshift);
    end
end
   
out.mshifts = mshifts;
out.mshiftsB = mshiftsB;
out.mshiftsF = mshiftsF;
out.fshiftsB = fshiftB;
out.fshiftsF = fshiftF;
out.fshiftsZ = fshiftZ;
out.xs = xs;
out.dists = dists;
out.ens_pos_shifts = all_pos_seq;
out.ens = get_quants(all_pos_seq);


function out = get_quants(all_pos_seq)

for an = 1:5
    for cn = 1:size(all_pos_seq,2)
            tr_seq = all_pos_seq{an,cn};
%             figure(1000);clf;imagesc(tr_seq);colorbar
            seq_shift = abs(diff(tr_seq,[],2));
            m_seq_shift = nanmean(seq_shift,2); % this is mean over trials;
            m_f_seq_shift_trials_LZ(an,cn) = nanmean(sum(seq_shift<0,2)./sum(~isnan(seq_shift),2));
            m_f_seq_shift_trials_GZ(an,cn) = nanmean(sum(seq_shift>0,2)./sum(~isnan(seq_shift),2));
            mm_seq_shift(an,cn) = nanmean(m_seq_shift); % this is mean over cells
            f_seq_shift_LZ(an,cn) = sum(m_seq_shift<0);
            f_seq_shift_GZ(an,cn) = sum(m_seq_shift>0);
            
            seq_shift = abs(diff(tr_seq(:,1:5),[],2));
            m_seq_shift1 = nanmean(seq_shift,2); % this is mean over trials;
            seq_shift = abs(diff(tr_seq(:,6:10),[],2));
            m_seq_shift2 = nanmean(seq_shift,2); % this is mean over trials;
            mm_seq_shift1(an,cn) = nanmean(m_seq_shift1); % this is mean over cells
            mm_seq_shift2(an,cn) = nanmean(m_seq_shift2); % this is mean over cells
            pLs = tr_seq;
            mshift = []; 
            for ii = 1:size(pLs,1)
                tcell = pLs(ii,:);
                ntcell = tcell(~isnan(tcell));
                mshift(ii,1) = mean(diff(ntcell));
            end
            mshifts(an,cn) = mean(mshift);
            mshiftsB(an,cn) = mean(mshift(mshift<0));
            mshiftsF(an,cn) = mean(mshift(mshift>0));
            fshiftB(an,cn) = sum(mshift<0)/length(mshift);
            fshiftF(an,cn) = sum(mshift>0)/length(mshift);
            fshiftZ(an,cn) = sum(mshift==0)/length(mshift);
    end
end

out.shift_trials_cells = [m_f_seq_shift_trials_LZ m_f_seq_shift_trials_GZ];%[f_seq_shift_GZ f_seq_shift_LZ];
% out.shift_trials_cells = [f_seq_shift_GZ];
out.shift_trials_cells = [mm_seq_shift1 mm_seq_shift2];%[f_seq_shift_GZ f_seq_shift_LZ];
% out.shift_trials_cells = [mm_seq_shift];%[f_seq_shift_GZ f_seq_shift_LZ];
out.shift_trials_cells = [mshiftsF abs(mshiftsB)];%[f_seq_shift_GZ f_seq_shift_LZ];
% out.shift_trials_cells = [fshiftF fshiftB fshiftZ];%[f_seq_shift_GZ f_seq_shift_LZ];
% out.shift_trials_cells = [mshifts];%[f_seq_shift_GZ f_seq_shift_LZ];
out.f_shifts = [fshiftF fshiftB fshiftZ];%[f_seq_shift_GZ f_seq_shift_LZ];