function [odata,mData] = getRasterDataT(iei,context,stimMarker)
mData = [];

for aa = 1:length(iei)
    ei = iei{aa};
    belt_length = ei.b.belt_length;
    datap = [];
    for pp = 1:length(ei.plane)
        [~,data] = getDataContexts(ei.plane{pp},context,stimMarker);
        for ii = 1:length(data)
            [aa pp ii]
            if isequal([aa pp ii],[4 1 1])
                n = 0;
            end
            mrfs = data{ii}.gauss_fit_on_mean;
            [rs,coeff] = getMRFS_vals(mrfs);
            as = coeff(:,1);
            bs = coeff(:,2);
            cs = coeff(:,3);
            PWs = 2.36*cs./sqrt(2)*(belt_length/50);
            data{ii}.rs = rs;
            data{ii}.coeff = coeff';
            data{ii}.pws = PWs';
            data{ii}.centers = (bs * belt_length/50)';
            data{ii}.peaks = as';
            data{ii}.sel = (rs > 0.5) & (bs > 1 & bs < 50) & (PWs > 5 & PWs < 120)';% & (data{ii}.SI > 4));
            data{ii}.placeCells5 = (rs > 0.5) & (bs > 1 & bs < 50)' & (PWs > 5 & PWs < 120)' & (data{ii}.SI > 5);
            data{ii}.placeCells4 = (rs > 0.5) & (bs > 1 & bs < 50)' & (PWs > 5 & PWs < 120)' & (data{ii}.SI > 4);
            data{ii}.placeCells3 = (rs > 0.5) & (bs > 1 & bs < 50)' & (PWs > 5 & PWs < 120)' & (data{ii}.SI > 3);
            data{ii}.percentPlaceCells5 = sum(data{ii}.placeCells5)/length(data{ii}.placeCells5);
            data{ii}.percentPlaceCells4 = sum(data{ii}.placeCells4)/length(data{ii}.placeCells4);
            data{ii}.percentPlaceCells3 = sum(data{ii}.placeCells3)/length(data{ii}.placeCells3);
            data{ii}.formula = mrfs.gauss1Formula;
            data{ii}.belt_length = belt_length;
            [h,pvals,mvals] = getResponsiveness(data{ii}.rasters);
            data{ii}.h = h';
            data{ii}.p = pvals';
            data{ii}.excR = (h & mvals(:,1) < mvals(:,2))';
            data{ii}.inhR = (h & mvals(:,1) > mvals(:,2))';
        end
        datap{pp} = data;
    end
    odata{aa} = datap;
end


function [h,pvals,mvals] = getResponsiveness(tempD)
for kk = 1:size(tempD,3)
    thisRaster = tempD(:,:,kk);
%     thisRaster1 = thisRaster(:,1:25); thisRaster1 = thisRaster1(~isnan(thisRaster1));
%     thisRaster2 = thisRaster(:,26:50); thisRaster2 = thisRaster2(~isnan(thisRaster2));
    thisRaster1 = nanmean(thisRaster(:,1:25),2);
    thisRaster2 = nanmean(thisRaster(:,26:50),2);
    try
        [~,pvals(kk,1),~,~] = ttest2(thisRaster1,thisRaster2);
%                 pvals(kk,1) = anova1(thisRaster,[ones(1,25) ones(1,25)*2],'off','alpha',0.01);
        mvals(kk,:) = [mean(thisRaster1) mean(thisRaster2)];
%         mvals(kk,:) = [mean(mean(thisRaster(:,1:25))) mean(mean(thisRaster(:,25:50)))];
    catch
        pvals(kk,1) = NaN;
        mvals(kk,:) = [NaN NaN];
    end
end
h(:,1) = pvals<0.05;
h(:,2) = pvals<0.01;
h(:,3) = pvals<0.001;

