function out = find_population_vector_trial_to_trial_corr(Rs,resp)


for rr = 1:size(Rs,1)
    R = Rs{rr};
    rasters = permute(R.sp_rasters1,[2 1 3]);
    sum_rasters_bins = squeeze(nansum(rasters,1));
    t_resp = resp{rr};
    for cn = 1:length(t_resp)
        if ~t_resp(cn)
            continue;
        end
        tsum = sum_rasters_bins(:,cn);
        tsumL = tsum > 0;
        mRs = transpose(rasters(:,tsumL,cn));
    end
    
end


% find binwidths and column sizes of all population activity to later fill
% in missing values and bring all correlation matrices to equal size
for rr = 1:size(mRs,1)
    for cc1 = 1:size(mRs,2)
        cols(rr,cc1) = size(mRs{rr,cc1},2);
        binwidths(rr,cc1) = Rs{rr,cc1}.bin_width;
    end
end

maxcolsz = max(cols(:));
bw = unique(binwidths(:));
if isfield(Rs{1,1}.resp,'cis')
    cis = Rs{1,1}.resp.cis;
    xs = Rs{1,1}.xs;
    xs = xs - xs(cis(1,2));
    xticks = [cis(1,:) maxcolsz];
    xs = round(xs);
    xso.vals = xs;
    xso.ticks = xticks;
    xso.label = 'Time (secs)';
else
    xs = 0:bw:1000;
    xs = xs(1:maxcolsz);
    cols = maxcolsz;
    colsHalf = round(cols/2);
    xso.vals = round(xs);
    xso.ticks = [1 colsHalf cols];
    xso.label = 'Position (cm)';
end

% correlations with means e.g., of all trial responses with mean of raster
% (same thing for accross conditions) - HaoRan's idea of looking at how
% similar each trial is to the mean
for rr = 1:size(mRs,1) % for each animal
    this_mean = NaN(size(mRs{rr,1},1),maxcolsz,size(mRs,2));
    for cc1 = 1:size(mRs,2)
        this_mean(:,1:size(mRs{rr,cc1},2),cc1) = mRs{rr,cc1};
    end
    mean_mRs{rr,1} = nanmean(this_mean,3);
    for cc1 = 1:size(mRs,2)
        RV1 = mRs{rr,cc1}; RV2 = mean_mRs{rr,1}; % get rate vectors from 1 condition and the other one is mean
        RV1 = correctsz(RV1(resp{rr},:),maxcolsz); RV2 = correctsz(RV2(resp{rr},:),maxcolsz); % make sizes equal
        [RV1_ordered,~,cellnums] = findPopulationVectorPlot(RV1,[]); % order RV1 according to peak firing
        [RV2_ordered,~,~] = findPopulationVectorPlot(RV2,[],cellnums); % order RV2 the same as RV1
        [this_SP_corr,pSP] = corr(RV1_ordered',RV2_ordered');
        this_SP_corr = fillmissing(this_SP_corr,'linear',2,'EndValues','nearest');
        this_SP_corr = fillmissing(this_SP_corr,'linear',1,'EndValues','nearest');
        SP_corr_with_mean{cc1,1} = this_SP_corr;
        SP_corr_diag_with_mean{cc1,1} = diag(this_SP_corr); % length equal to number of neurons
    end
    all_SP_corr_with_mean{rr} = SP_corr_with_mean;
    all_SP_corr_diag_with_mean{rr} = SP_corr_diag_with_mean;
end

for rr = 1:size(mRs,1) % for each animal
    PV_corr = []; PV_corr_diag = [];
    SP_corr =[]; SP_corr_diag = [];
    RR = []; RR_SP = []; RR_PV = []; RR_score_PV_corr_diag = []; RR_score_PF_corr_diag = [];
    for cc1 = 1:size(mRs,2) % for each population activity (rate vectors) matrix
        for cc2 = 1:size(mRs,2) % for each population activity matrix --> comparison with itself and all others
            RV1 = mRs{rr,cc1}; RV2 = mRs{rr,cc2}; % get rate vectors from two conditions
            RV1 = correctsz(RV1(resp{rr},:),maxcolsz); RV2 = correctsz(RV2(resp{rr},:),maxcolsz); % make sizes equal
            [RV1_ordered,~,cellnums] = findPopulationVectorPlot(RV1,[]); % order RV1 according to peak firing
            [RV2_ordered,~,~] = findPopulationVectorPlot(RV2,[],cellnums); % order RV2 the same as RV1
            [this_PV_corr,pPV] = corr(RV1_ordered,RV2_ordered); % find correlation
            this_PV_corr = fillmissing(this_PV_corr,'linear',2,'EndValues','nearest'); 
            this_PV_corr = fillmissing(this_PV_corr,'linear',1,'EndValues','nearest');
            PV_corr{cc1,cc2} = this_PV_corr;
            PV_corr_diag{cc1,cc2} = diag(this_PV_corr);% length equal to number of bins (for autocorrelation all values would be 1) 
            %the diagonal indicates the corresponding bin comparison in two conditions
%             [param.FD(cc1,cc2,rr),param.all_bw{cc1,cc2}] = findFractalDim(this_PV_corr);
            
            % spatial correlation
            [this_SP_corr,pSP] = corr(RV1_ordered',RV2_ordered');
            this_SP_corr = fillmissing(this_SP_corr,'linear',2,'EndValues','nearest');
            this_SP_corr = fillmissing(this_SP_corr,'linear',1,'EndValues','nearest');
            SP_corr{cc1,cc2} = this_SP_corr;
            SP_corr_diag{cc1,cc2} = diag(this_SP_corr); % length equal to number of neurons
            % the diagnoal corresponds to the same cell correlation
            % comparison which means we are looking at how the same cell
            % behaved in two conditions
            
            
            % RR --> rate remapping
            RR_score = abs(RV1_ordered-RV2_ordered)./(RV1_ordered+RV2_ordered);
            RR{cc1,cc2} = RR_score;
            RR_score_PV_corr = corr(RR_score);
            RR_score_PF_corr = corr(RR_score');
            RR_SP{cc1,cc2} = nanmean(RR_score,2);
            RR_PV{cc1,cc2} = [nanmean(RR_score,1)]';
            RR_score_PV_corr_diag{cc1,cc2} = diag(RR_score_PV_corr);
            RR_score_PF_corr_diag{cc1,cc2} = diag(RR_score_PF_corr);
        end
   end
    all_PV_corr{rr} = PV_corr;
    all_PV_corr_diag{rr} = PV_corr_diag;
    all_SP_corr{rr} = SP_corr;
    all_SP_corr_diag{rr} = SP_corr_diag;
    all_RR{rr} = RR;
    all_RR_PV{rr} = RR_PV;
    all_RR_SP{rr} = RR_SP;
    all_RR_score_PV_corr_diag{rr} = RR_score_PV_corr_diag;
    all_RR_score_PF_corr_diag{rr} = RR_score_PF_corr_diag;
    
    ind = 1; auto_SP_corr = [];
    for cc1 = 1:size(mRs,2)
        for cc2 = 1:size(mRs,2)
            if cc1 == cc2
                auto_SP_corr{ind} = SP_corr{cc1,cc2}; ind = ind + 1;
            end
        end
    end
    auto_SP_corr_corr = [];
    mask = logical(triu(ones(size(auto_SP_corr{1})),1));
    for cc1 = 1:length(auto_SP_corr)
        for cc2 = 1:length(auto_SP_corr)
            set1 = auto_SP_corr{cc1}; set2 = auto_SP_corr{cc2};
            auto_SP_corr_corr(cc1,cc2) = corr(set1(mask),set2(mask));
        end
    end
    all_auto_SP_corr_corr(:,:,rr) = auto_SP_corr_corr;
end

N = size(mRs,2);
mean_corr = cell(N,N);
mask = ones(size(mean_corr)); mask = triu(mask,1) & ~triu(mask,2);
ind = 1;
for rr = 1:N
    for cc = 1:N
        for_mean_PV_corr = [];
        for an = 1:length(all_PV_corr)
            for_mean_PV_corr(:,:,an) = all_PV_corr{an}{rr,cc};
        end
        mean_PV_corr{rr,cc} = nanmean(for_mean_PV_corr,3);
        
        % for across adjacent conditions
        if mask(rr,cc)
            for an = 1:length(all_PV_corr)
                adj_PV_corr_diag{an,ind} = all_PV_corr_diag{an}{rr,cc};
                adj_SP_corr_diag{an,ind} = all_SP_corr_diag{an}{rr,cc};
                adj_RR_SP{an,ind} = all_RR_SP{an}{rr,cc};
                adj_RR_PV{an,ind} = all_RR_PV{an}{rr,cc};
            end
            ind = ind + 1;
        end
    end
end
out.mean_PV_corr = mean_PV_corr;
out.all_SP_corr = all_SP_corr;
out.all_SP_corr_diag = all_SP_corr_diag;
out.all_PV_corr = all_PV_corr;
out.all_PV_corr_diag = all_PV_corr_diag;
out.all_RR = all_RR;
out.all_RR_SP = all_RR_SP;
out.all_RR_PV = all_RR_PV;
out.adj_PV_corr_diag = adj_PV_corr_diag;
out.adj_SP_corr_diag = adj_SP_corr_diag;
out.adj_RR_SP = adj_RR_SP;
out.adj_RR_PV = adj_RR_PV;
out.all_auto_SP_corr_corr = all_auto_SP_corr_corr;
out.xs = xso;
out.all_SP_corr_with_mean = all_SP_corr_with_mean;
out.all_SP_corr_diag_with_mean = all_SP_corr_diag_with_mean;


n = 0;




function RV1 = correctsz(RV1,maxcolsz)
if size(RV1,2) < maxcolsz
    diffsz = maxcolsz - size(RV1,2);
    nanmat = NaN(size(RV1,1),diffsz);
    RV1 = [RV1 nanmat];
end

function [fd,bw] = findFractalDim(corrV)
bw = imbinarize(corrV,'adaptive','Sensitivity',0.75);
fd = BoxCountfracDim(bw);
