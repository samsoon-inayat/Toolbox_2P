function out = find_population_vector_corr_remap(Rs,mRs,resp)
% function [all_corr_an,all_corr_cell_an,mean_corr,mean_corr_popV,mean_corr_cell,xso,param] = find_population_vector_corr_remap(Rs,mRs,resp)

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
    PV_corr = []; PV_corr_diag = []; PV_corr_ind = []; PV_corr_diag_ind = [];
    SP_corr =[]; SP_corr_diag = []; SP_corr_ind =[]; SP_corr_diag_ind = [];
    RR = []; RR_SP = []; RR_PV = []; RR_score_PV_corr_diag = []; RR_score_PF_corr_diag = [];
    pop_sim = []; pop_sim_ind = [];
    for cc1 = 1:size(mRs,2) % for each population activity (rate vectors) matrix
        for cc2 = 1:size(mRs,2) % for each population activity matrix --> comparison with itself and all others
            RV1 = mRs{rr,cc1}; RV2 = mRs{rr,cc2}; % get rate vectors from two conditions
            RV1 = correctsz(RV1(resp{rr},:),maxcolsz); RV2 = correctsz(RV2(resp{rr},:),maxcolsz); % make sizes equal
            [RV1_ordered,~,cellnums] = findPopulationVectorPlot(RV1,[]); % order RV1 according to peak firing
            [RV2_ordered,~,~] = findPopulationVectorPlot(RV2,[],cellnums); % order RV2 the same as RV1
            [RV2_ordered_ind,~,cellnums_ind] = findPopulationVectorPlot(RV2,[]); % order RV2 the same as RV1
            [this_PV_corr,pPV] = corr(RV1_ordered,RV2_ordered); % find correlation
            [this_PV_corr_ind,pPV] = corr(RV1_ordered,RV2_ordered_ind); % find correlation
            
            R2d(cc1,cc2) = corr2(RV1_ordered,RV2_ordered); % find correlation
            R2d_ind(cc1,cc2) = corr2(RV1_ordered,RV2_ordered_ind); % find correlation
            cell_pos_shift(cc1,cc2) = mean(find_shifts(cellnums,cellnums_ind));
            
            this_PV_corr = fillmissing(this_PV_corr,'linear',2,'EndValues','nearest'); 
            this_PV_corr = fillmissing(this_PV_corr,'linear',1,'EndValues','nearest');
            
            this_PV_corr_ind = fillmissing(this_PV_corr_ind,'linear',2,'EndValues','nearest'); 
            this_PV_corr_ind = fillmissing(this_PV_corr_ind,'linear',1,'EndValues','nearest');
            
            PV_corr{cc1,cc2} = this_PV_corr;
            PV_corr_diag{cc1,cc2} = diag(this_PV_corr);% length equal to number of bins (for autocorrelation all values would be 1) 
            %the diagonal indicates the corresponding bin comparison in two conditions
%             [param.FD(cc1,cc2,rr),param.all_bw{cc1,cc2}] = findFractalDim(this_PV_corr);

            PV_corr_ind{cc1,cc2} = this_PV_corr_ind;
            PV_corr_diag_ind{cc1,cc2} = diag(this_PV_corr_ind);% length equal to number of bins (for autocorrelation all values would be 1) 
            
            % spatial correlation
            [this_SP_corr,pSP] = corr(RV1_ordered',RV2_ordered');
            this_SP_corr = fillmissing(this_SP_corr,'linear',2,'EndValues','nearest');
            this_SP_corr = fillmissing(this_SP_corr,'linear',1,'EndValues','nearest');
            SP_corr{cc1,cc2} = this_SP_corr;
            SP_corr_diag{cc1,cc2} = diag(this_SP_corr); % length equal to number of neurons
            % the diagnoal corresponds to the same cell correlation
            % comparison which means we are looking at how the same cell
            % behaved in two conditions
            
            [this_SP_corr_ind,pSP] = corr(RV1_ordered',RV2_ordered_ind');
            this_SP_corr_ind = fillmissing(this_SP_corr_ind,'linear',2,'EndValues','nearest');
            this_SP_corr_ind = fillmissing(this_SP_corr_ind,'linear',1,'EndValues','nearest');
            SP_corr_ind{cc1,cc2} = this_SP_corr_ind;
            SP_corr_diag_ind{cc1,cc2} = diag(this_SP_corr_ind); % length equal to number of neurons
            
            average_mR1 = mean(RV1_ordered,2); average_mR2 = mean(RV2_ordered,2); average_mR2_ind = mean(RV2_ordered_ind,2);
            
%             figure(10000);clf;plot(1:351,average_mR1);hold on;plot(1:351,average_mR2);plot(1:351,average_mR2_ind);
            
            pop_sim(cc1,cc2) = corr(average_mR1,average_mR2);
            pop_sim_ind(cc1,cc2) = corr(average_mR1,average_mR2_ind);
            
            
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
    all_R2d{rr} = R2d;
    all_R2d_ind{rr} = R2d_ind;
    all_cell_pos_shift{rr} = cell_pos_shift;
    
    all_PV_corr_ind{rr} = PV_corr_ind;
    all_PV_corr_diag_ind{rr} = PV_corr_diag_ind;
    all_SP_corr_ind{rr} = SP_corr_ind;
    all_SP_corr_diag_ind{rr} = SP_corr_diag_ind;
    
    all_pop_sim{rr} = pop_sim;
    all_pop_sim_ind{rr} = pop_sim_ind;
    
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
mean_corr_ind = mean_corr;
ind = 1;
for rr = 1:N
    for cc = 1:N
        for_mean_PV_corr = [];
        for_mean_PV_corr_ind = [];
        for an = 1:length(all_PV_corr)
            for_mean_PV_corr(:,:,an) = all_PV_corr{an}{rr,cc};
            for_mean_PV_corr_ind(:,:,an) = all_PV_corr_ind{an}{rr,cc};
        end
        mean_PV_corr{rr,cc} = nanmean(for_mean_PV_corr,3);
        mean_PV_corr_ind{rr,cc} = nanmean(for_mean_PV_corr_ind,3);
        
        % for across adjacent conditions
        if mask(rr,cc)
            for an = 1:length(all_PV_corr)
                adj_PV_corr_diag{an,ind} = all_PV_corr_diag{an}{rr,cc};
                adj_SP_corr_diag{an,ind} = all_SP_corr_diag{an}{rr,cc};
                adj_PV_corr_diag_ind{an,ind} = all_PV_corr_diag_ind{an}{rr,cc};
                adj_SP_corr_diag_ind{an,ind} = all_SP_corr_diag_ind{an}{rr,cc};
                adj_RR_SP{an,ind} = all_RR_SP{an}{rr,cc};
                adj_RR_PV{an,ind} = all_RR_PV{an}{rr,cc};
                adj_R2d{an,ind} = all_R2d{an}(rr,cc);
                adj_R2d_ind{an,ind} = all_R2d_ind{an}(rr,cc);
                adj_cell_pos_shift{an,ind} = all_cell_pos_shift{an}(rr,cc);
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
out.all_R2d = all_R2d;
out.all_R2d_ind = all_R2d_ind;
out.all_cell_pos_shift = all_cell_pos_shift;

out.mean_PV_corr_ind = mean_PV_corr_ind;
out.all_SP_corr_ind = all_SP_corr_ind;
out.all_SP_corr_diag_ind = all_SP_corr_diag_ind;
out.all_PV_corr_ind = all_PV_corr_ind;
out.all_PV_corr_diag_ind = all_PV_corr_diag_ind;

out.all_RR = all_RR;
out.all_RR_SP = all_RR_SP;
out.all_RR_PV = all_RR_PV;
out.adj_PV_corr_diag = adj_PV_corr_diag;
out.adj_SP_corr_diag = adj_SP_corr_diag;

out.adj_PV_corr_diag_ind = adj_PV_corr_diag_ind;
out.adj_SP_corr_diag_ind = adj_SP_corr_diag_ind;

out.all_pop_sim = all_pop_sim;
out.all_pop_sim_ind = all_pop_sim_ind;

out.adj_R2d = adj_R2d;
out.adj_R2d_ind = adj_R2d_ind;
out.adj_cell_pos_shift = adj_cell_pos_shift;

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

function sh = find_shifts(c1,c2)
sh = NaN(size(c1));
tot = length(c1);
parfor ii = 1:length(c1)
    ind = find(c2 == c1(ii));
    sh(ii) = 100*abs(ii - ind)/tot;
end
n=0;