function [all_corr_an,all_corr_cell_an,mean_corr,mean_corr_cell,xs] = find_population_vector_corr_remap(Rs,mRs,resp)

for rr = 1:size(mRs,1)
    for cc1 = 1:size(mRs,2)
        cols(rr,cc1) = size(mRs{rr,cc1},2);
        binwidths(rr,cc1) = Rs{rr,cc1}.bin_width;
    end
end

maxcolsz = max(cols(:));
bw = unique(binwidths(:));
xs = 3:bw:1000;
xs = xs(1:maxcolsz);

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
            [corrCV,pCV] = corr(pv1',pv2');
            corrCV = fillmissing(corrCV,'linear',2,'EndValues','nearest');
            corrCV = fillmissing(corrCV,'linear',1,'EndValues','nearest');
            all_corr_cell{cc1,cc2} = corrCV;
        end
   end
    all_corr_an{rr} = all_corr;
    all_corr_cell_an{rr} = all_corr_cell;
end

N = size(mRs,2);
mean_corr = cell(N,N);
for rr = 1:N
    for cc = 1:N
        thisCC = [];
        for an = 1:length(all_corr_an)
            thisCC(:,:,an) = all_corr_an{an}{rr,cc};
            mean_corr_cell(rr,cc,an) = mean(diag(all_corr_cell_an{an}{rr,cc}));
        end
        mean_corr{rr,cc} = nanmean(thisCC,3);
    end
end


function set1 = correctsz(set1,maxcolsz)
if size(set1,2) < maxcolsz
    diffsz = maxcolsz - size(set1,2);
    nanmat = NaN(size(set1,1),diffsz);
    set1 = [set1 nanmat];
end