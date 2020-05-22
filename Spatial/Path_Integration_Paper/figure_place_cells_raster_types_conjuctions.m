function figure_place_cells_vs_other_cells(fn,allRs,ccs)

ei = evalin('base','ei10');
mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;
% selAnimals = 1;
% selAnimals = 5:8;
selAnimals = [1:4 9];
mData.belt_length = ei{selAnimals(1)}.b.belt_length;
n = 0;

%%
selCells = 'areCells';
varName = 'rasters';
planeNumbers = 'All';
maxDistTime = [142 15];
contextNumbers = ones(1,4)*1;
stimMarkers = {'air','air','belt','airI'};
rasterTypes = {'dist','time','dist','time'};
% contextNumbers = [1 2 3 4];
% stimMarkers = {'air','air','air','air'};
% rasterTypes = {'dist','dist','dist','dist'};
trials = 3:10;
trials10 = 3:9;
% align cells
CNi = 1;
% select spatial cells
sNi = 1;
% cNi = sNi;
for ii = 1:length(contextNumbers)
    contextNumber = contextNumbers(ii);
    mRsi = []; distDi = [];
    for jj = 1:length(selAnimals)
        disp([ii jj]);
        if isequal([ii jj],[1 5])
            n = 0;
        end
%         [pcs_d cns areCells] = getParamValues('placeCells5',ei(selAnimals(jj)),planeNumbers,contextNumber,'airI','dist',selCells,maxDistTime);
        [pcs cns areCells] = getParamValues('placeCells3',ei(selAnimals(jj)),planeNumbers,contextNumber,stimMarkers{ii},rasterTypes{ii},selCells,maxDistTime);
%         [clus] = getParamValues('cluster3',ei(selAnimals(jj)),planeNumbers,contextNumber,stimMarkers{sNi},rasterTypes{sNi},selCells,maxDistTime);
%         pcs = pcs_d | pcs_t;
%         pcs = logical(ones(size(pcs)));
        [tempD cns] = getParamValues(varName,ei(selAnimals(jj)),planeNumbers,contextNumber,stimMarkers{ii},rasterTypes{ii},selCells,maxDistTime);
        [zMIs cns] = getParamValues('info_metrics.ShannonMI_Zsh',ei(selAnimals(jj)),planeNumbers,contextNumber,stimMarkers{ii},rasterTypes{ii},selCells,maxDistTime);
        [data cns areCells] = getParamValues('',ei(selAnimals(jj)),planeNumbers,contextNumber,stimMarkers{ii},rasterTypes{ii},selCells,maxDistTime);
%         cellSel = clus(areCells) == 1;
        cellSel = pcs;
%         cellSel = 1:sum(areCells)';
        all_pcs{jj,ii} = pcs;
        npcs(jj,ii) = 100*sum(pcs)/length(pcs);
        distDi = [distDi;cellSel];
        try
             mR = findMeanRasters(tempD,trials);
         catch
             mR = findMeanRasters(tempD,trials10);
         end
         mRsi = [mRsi;mR];
    end
    distD{ii} = distDi;
     allRs{ii} = mRsi;
    dxs = diff(data.xs); bin_width = dxs(1); xs = 0:bin_width:1000;
    time_xs{ii} = xs(1:size(mRsi,2));
end

npcs(:,1)
combs = nchoosek(1:size(all_pcs,2),2);
for ii = 1:size(all_pcs,1)
    thisA = all_pcs(ii,:);
    for jj = 1:size(combs,1)
        pcs1 = thisA{combs(jj,1)}; pcs2 = thisA{combs(jj,2)};
        conjuctive(ii,jj) = 100 * sum(pcs1 & pcs2)/length(pcs1);
    end
end
n = 0;
figure(100);clf;imagesc(conjuctive);colorbar;
set(gca,'XTickLabels',{'1-2','1-3','1-4','2-3','2-4','3-4'},'YTick',[1 2 3 4 5]);
ylabel('Animal Number');
xtickangle(20);title('Percent overlap');
%%
