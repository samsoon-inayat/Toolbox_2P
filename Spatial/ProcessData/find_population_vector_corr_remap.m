function out = find_population_vector_corr_remap(Rs,mRs,resp)
% function [all_corr_an,all_corr_cell_an,mean_corr,mean_corr_popV,mean_corr_cell,xso,param] = find_population_vector_corr_remap(Rs,mRs,resp)

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

for rr = 1:size(mRs,1)
    PV_corr = []; PV_corr_diag = [];
    SP_corr =[]; SP_corr_diag = [];
    RR = []; RR_SP = []; RR_PV = []; RR_score_PV_corr_diag = []; RR_score_PF_corr_diag = [];
    for cc1 = 1:size(mRs,2)
        for cc2 = 1:size(mRs,2)
            RV1 = mRs{rr,cc1}; RV2 = mRs{rr,cc2};
            RV1 = correctsz(RV1(resp{rr},:),maxcolsz); RV2 = correctsz(RV2(resp{rr},:),maxcolsz);
            [RV1_ordered,~,cellnums] = findPopulationVectorPlot(RV1,[]);
            [RV2_ordered,~,~] = findPopulationVectorPlot(RV2,[],cellnums);
            [this_PV_corr,pPV] = corr(RV1_ordered,RV2_ordered);
            this_PV_corr = fillmissing(this_PV_corr,'linear',2,'EndValues','nearest');
            this_PV_corr = fillmissing(this_PV_corr,'linear',1,'EndValues','nearest');
            PV_corr{cc1,cc2} = this_PV_corr;
            PV_corr_diag{cc1,cc2} = diag(this_PV_corr);
            
%             [param.FD(cc1,cc2,rr),param.all_bw{cc1,cc2}] = findFractalDim(this_PV_corr);
            
            % spatial correlation
            [this_SP_corr,pSP] = corr(RV1_ordered',RV2_ordered');
            this_SP_corr = fillmissing(this_SP_corr,'linear',2,'EndValues','nearest');
            this_SP_corr = fillmissing(this_SP_corr,'linear',1,'EndValues','nearest');
            SP_corr{cc1,cc2} = this_SP_corr;
            SP_corr_diag{cc1,cc2} = diag(this_SP_corr);
            
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
out.xs = xso;


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
