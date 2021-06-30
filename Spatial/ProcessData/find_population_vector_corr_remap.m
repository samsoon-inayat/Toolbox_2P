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
    all_corr = [];
    all_corr_cell =[];
    for cc1 = 1:size(mRs,2)
        for cc2 = 1:size(mRs,2)
            set1 = mRs{rr,cc1}; set2 = mRs{rr,cc2};
            set1 = correctsz(set1(resp{rr},:),maxcolsz); set2 = correctsz(set2(resp{rr},:),maxcolsz);
            [pv1,~,cellnums] = findPopulationVectorPlot(set1,[]);
            [pv2,~,~] = findPopulationVectorPlot(set2,[],cellnums);
            [corrV,pV] = corr(pv1,pv2);
            corrV = fillmissing(corrV,'linear',2,'EndValues','nearest');
            corrV = fillmissing(corrV,'linear',1,'EndValues','nearest');
            all_corr{cc1,cc2} = corrV;
            
            [param.FD(cc1,cc2,rr),param.all_bw{cc1,cc2}] = findFractalDim(corrV);
            
            % spatial correlation
            [corrCV,pCV] = corr(pv1',pv2');
            corrCV = fillmissing(corrCV,'linear',2,'EndValues','nearest');
            corrCV = fillmissing(corrCV,'linear',1,'EndValues','nearest');
            all_corr_cell{cc1,cc2} = corrCV;
            
            % rate remapping
            rate_remapping{cc1,cc2} = abs(pv1-pv2)./(pv1+pv2);
            rate_remapping_cells{cc1,cc2} = nanmean(abs(pv1-pv2)./(pv1+pv2),2);
            rate_remapping_pv{cc1,cc2} = nanmean(abs(pv1-pv2)./(pv1+pv2),1);
        end
   end
    all_corr_an{rr} = all_corr;
    all_corr_cell_an{rr} = all_corr_cell;
    rate_remapping_an{rr} = rate_remapping;
    rate_remapping_an_cell{rr} = rate_remapping_cells;
    rate_remapping_an_pv{rr} = rate_remapping_pv;
end

N = size(mRs,2);
mean_corr = cell(N,N);
for rr = 1:N
    for cc = 1:N
        thisCC = [];
        for an = 1:length(all_corr_an)
            temp_popV = all_corr_an{an}{rr,cc};
            thisCC(:,:,an) = temp_popV;
            mean_corr_popV(rr,cc,an) = nanmean(diag(temp_popV));
            corr_popV{rr,cc,an} = diag(temp_popV);
            mean_corr_cell(rr,cc,an) = mean(diag(all_corr_cell_an{an}{rr,cc}));
            corr_cell{rr,cc,an} = diag(all_corr_cell_an{an}{rr,cc});
            mean_rate_remap_cells(rr,cc,an) = nanmean(rate_remapping_an_cell{an}{rr,cc});
            mean_rate_remap_pv(rr,cc,an) = nanmean(rate_remapping_an_pv{an}{rr,cc});
        end
        mean_corr{rr,cc} = nanmean(thisCC,3);
    end
end
out.popV.animals = all_corr_an;
out.popV.mean_corr_animal_wise = mean_corr;
out.popV.mean_popV_wise = mean_corr_popV;
out.popV.corr = corr_popV;
out.popV.param = param;
out.spatial.animals = all_corr_cell_an;
out.spatial.mean_cell_wise = mean_corr_cell;
out.spatial.corr = corr_cell;
out.rate_remap.animals = rate_remapping_an;
out.rate_remap.animals_cells = rate_remapping_an_cell;
out.rate_remap.animals_pv = rate_remapping_an_pv;
out.rate_remap.mean_cells = mean_rate_remap_cells;
out.rate_remap.mean_pv = mean_rate_remap_pv;
out.xs = xso;


n = 0;




function set1 = correctsz(set1,maxcolsz)
if size(set1,2) < maxcolsz
    diffsz = maxcolsz - size(set1,2);
    nanmat = NaN(size(set1,1),diffsz);
    set1 = [set1 nanmat];
end

function [fd,bw] = findFractalDim(corrV)
bw = imbinarize(corrV,'adaptive','Sensitivity',0.75);
fd = BoxCountfracDim(bw);
