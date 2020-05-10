function allA_P = getRasters_air_belt(ei,trials,ow)

if ~exist('ow','var')
    ow = 0;
end

if iscell(ei)
    for ii = 1:length(ei)
        temp = getRasters_air_belt(ei{ii},trials,ow);
        if ii == 1
            cellRaster = temp.cells;
            rasters = temp.rasters;
        else
            cellRaster = [cellRaster temp.cells];
            rasters = cat(3,rasters,temp.rasters);
        end
    end
    allA_P = temp;
    allA_P.cells = cellRaster;
    allA_P.rasters = rasters;
    return;
end

if ~exist('ei','var')
    ei = evalin('base','ei{1}');    
end

onsets = ei.b.air_puff_r(trials);
offsets = ei.b.air_puff_f(trials);
photo_sensor = ei.b.photo_sensor_f(ei.b.photo_sensor_f>onsets(1) & ei.b.photo_sensor_f<offsets(end));
dists = ei.b.dist(photo_sensor);
diff_dists = diff(dists);
inds = find(diff_dists < 100);
temp_photo_sensor = photo_sensor;
temp_photo_sensor(inds+1) = [];


a_p_onsets = onsets(2:end);
a_p_offsets = temp_photo_sensor(2:end);
% plotMarkers(ei,a_p_onsets,a_p_offsets,101);
allA_P = getDistRaster(ei,a_p_onsets,a_p_offsets,ow,0);

% 
% function rasters =  getDistRaster_2(ei,onsets,offsets,overwrite,plotFlag)
% 
% fileName = makeName(sprintf('raster_%d_%d_%d_%d.mat',onsets(1),onsets(end),offsets(1),offsets(end)),ei.folders.thispFolder);
% if exist(fileName,'file') && overwrite == 0
%     load(fileName);
%     rasters = findOccupancyNormalizedRasters(rasters,ei.deconv.spSigAll);
%     return;
% end
% 
% ccs = ei.areCells;
% b = ei.b;
% hWaitBar = waitbar(0,sprintf('Processed cell -'));
% for ii = 1:length(ccs)
%     waitbar(ii/length(ccs),hWaitBar,sprintf('Processing cell %d/%d',ii,length(ccs)));
%     tsp = ei.deconv.spSigAll{ii}';
% %     caSig = ei.signals(ii,:);
%     sCaSig = ei.deconv.caSigAll{ii}';
%     response_rasters(ii) = getDistRaster_1(b,sCaSig,tsp,onsets,offsets,plotFlag);
% end
% close(hWaitBar);
% 
% out = getDistRaster_forSpeed(b,onsets,offsets);
% rasters.cells = response_rasters;
% rasters.speed = out.distSpeedRaster;
% rasters.duration = out.distDurRaster;
% rasters.dist = out.dist;
% rasters.onsets = onsets;
% rasters.offsets = offsets;
% rasters = findOccupancyNormalizedRasters(rasters,ei.deconv.spSigAll);
% save(fileName,'rasters','-v7.3');
% 
% 
% function out =  getDistRaster_1(b,caSig,spSignal,onsets,offsets,plotFlag)
% 
% trialNumbers = 1:length(onsets);
% for iii = 1:length(trialNumbers)
%     ii = trialNumbers(iii);
%     st = onsets(ii);
%     se = offsets(ii);
%     frames = (find(b.frames_f >= st & b.frames_f <= se)); % find frame numbers
%     oSignals{iii} = caSig(frames);
%     spoSignals{iii} = spSignal(frames);
%     maxSignals(iii) = max(oSignals{iii});
%     minSignals(iii) = min(oSignals{iii});
%     ts{iii} = b.ts(b.frames_f(frames))-b.ts(st);
% %     temp = b.fSpeed(st:se);
% %     temp(temp<0) = 0;
% %     speed{iii} = temp;
%     tss{iii} = b.ts(st:se)-b.ts(st);
%     dist{iii} = b.dist(st:se)-b.dist(st);
%     distVal(iii) = dist{iii}(end);
% %     if isfield(b,'tone_light_stim')
% %         tempTL = (find(b.tone_light_stim_r >= st & b.tone_light_stim_r <= se)); % find frame numbers
% %         try
% %             tone_light{iii} = tempTL(3:4);
% %         catch
% %             tone_light{iii} = tempTL(1:2);
% %         end
% %         tstl{iii} = b.ts(b.tone_light_stim_r(tone_light{iii}))-b.ts(st);
% %         dstl(iii,:) = b.dist(b.tone_light_stim_r(tone_light{iii}))-b.dist(st);
% %     end
% end
% minDist = min(distVal);
% nbins = 50;
% bins = linspace(0,minDist,(nbins+1));
% sbins = bins(1:nbins);
% ebins = bins(2:end);
% 
% %%% Method %%%
% % 1) Divide the minDist into nbins
% 
% for iii = 1:length(trialNumbers)
%     thisDist = dist{iii};
%     thisT = tss{iii};
%     thisTSig = ts{iii};
%     thisspSign = spoSignals{iii};
%     thiscaSign = oSignals{iii};
% %     thisSpeed = speed{iii};
%     for jj = 1:nbins
%         st = sbins(jj);  se = ebins(jj);
%         dVals = find(thisDist >= st & thisDist <se);
%         tVals = thisT(dVals);
%         stt = tVals(1);
%         sett = tVals(end);
%         tValsInd = find(thisTSig >= stt & thisTSig < sett);
% %         if isnan(mean(thisspSign(tValsInd))/(sett-stt))
% %             error('NaN in raster calculation');
% %         end
%         distSigRaster(iii,jj) = mean(thisspSign(tValsInd));
%         distSpikesRaster(iii,jj) = sum(thisspSign(tValsInd))*(sett-stt);
%         distCaRaster(iii,jj) = mean(thiscaSign(tValsInd));
% %         distSigRaster(iii,jj) = mean(thisspSign(tValsInd))/(sett-stt);
% %         distDurRaster(iii,jj) = sett-stt;
% %         distSpeedRaster(iii,jj) = mean(thisSpeed(dVals));
%     end
% %     dSigF(iii,:) = applyGaussFilt(distSigRaster(iii,:),3);
% end
% 
% out.raster_meanFiringRate = distSigRaster;
% out.raster_spikeCount = distSpikesRaster;
% out.raster_df = distCaRaster;
% % out.distSigRaster = dSigF;
% % out.raster_duration = distDurRaster;
% % out.dist = sbins;
% % out.distSpeedRaster = distSpeedRaster;
% % out.distSpikesRaster = distSpikesRaster;
% % out.raster = distSigRaster./distDurRaster;
% % out.SI = spatial_information(b.ts(b.frames_f),spSignal,distSpikesRaster,distSigRaster,distDurRaster,'MI_formula');
% % out.MI = mutualinformationx(nanmean(distSigRaster),nanmean(distDurRaster));
% 
% if plotFlag
% n=0;
% mmSig = max(maxSignals);
% mnSig = min(minSignals);
% ff = makeFigureRowsCols(555,[0.1 -0.5 6 9],'RowsCols',[length(trialNumbers) 1],'spaceRowsCols',[0.021 0.0009],'rightUpShifts',[0.05 0.05],'widthHeightAdjustment',[-70 -25]);
% for ii = 1:length(trialNumbers)
%     axes(ff.h_axes(ii,1));
%     plot(ts{ii},oSignals{ii});hold on;
%     plot(ts{ii},spoSignals{ii});
%     tempSpeed = speed{ii}*0.45;
%     plot(tss{ii},tempSpeed);
%     if isfield(b,'tone_light_stim')
%         plot(tstl{ii},ones(size(tstl{ii})) * mmSig/2,'.');
%     end
%     ylim([mnSig mmSig]);
%     box off;
%     ylabel(ii);
%     set(ff.h_axes(ii,1),'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
% end
% 
% % mdstl = mean(dstl);
% % toneDist = mdstl(1);
% % lightDist = mdstl(2);
% % binI = bins-toneDist;
% % tBinI = find(binI>0,1,'first');
% % binI = bins-lightDist;
% % lBinI = find(binI>0,1,'first');
% 
% figure(666);clf;
% imagesc(distSigRaster./distDurRaster);
% colorbar;
% xVs = get(gca,'XTick');
% set(gca,'XTickLabel',round(bins(xVs+1)));
% text(tBinI,1,'T','color','r');
% text(lBinI,1,'L','color','r');
% figure(667);clf;
% mdSig = nanmean(distSigRaster);
% dSigT = applyGaussFilt(mdSig,3);
% plot(dSigT);
% end
% 
% function rasters = findOccupancyNormalizedRasters(rasters,spSignals)
% % for ii = 1:length(rasters.cells)
% %     spSigs(ii,:) = spSignals{ii}';
% % end
% for ii = 1:length(rasters.cells)
% %     thisSR = spSignals{ii};
% %     clus = kmeans(thisSR,2);
% %     mFRs = [mean(thisSR(clus==1)) mean(thisSR(clus==2))];
% %     thresh = max(mFRs);
% %     thresh = mean(thisSR)+3*std(thisSR);
% %     thresh = 0;
% %     rasterMFR = smoothts(rasters.cells(ii).raster_meanFiringRate,'g',3);%>thresh;
% %     thresh = 0;%mean(rasterMFR(:))+std(rasterMFR(:));
% %     rasterMFR1 = rasterMFR.*(rasterMFR>thresh);
% %     rasters.rasters(:,:,ii) = smoothts(rasterMFR./rasters.duration,'g',3);
% %     figure(100);clf;
% %     imagesc(rasters.rasters(:,:,ii));
% %     n=0;
%     rasters.rasters(:,:,ii) = rasters.cells(ii).raster_meanFiringRate./rasters.duration;
%     rasters.caRasters(:,:,ii) = rasters.cells(ii).raster_df./rasters.duration;
% end
% 
% 
% 
% 
% 
