function figure1_Distributions_df_by_fo

% dataAir = evalin('base','data');
% dataBelt = evalin('base','datab');
% dataAOn = evalin('base','dataAOn010');
% dataAOff = evalin('base','dataAOff010');
ei = evalin('base','ei10');
mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;
axes_font_size = mData.axes_font_size;
% allCells = mData.allCells;
selAnimals = 1:4;
% selAnimals = 5:8;
selAnimals = 1:8;
n = 0;

%%
runthis = 1;
if runthis
    planeNumbers = 'All';
    maxDistTime = [140 5];
    contextNumber = ones(1,4).*1;%[1 2 3 4];
    stimMarkers = {'air','belt','air','airI'};
    rasterTypes = {'dist','dist','time','time'};
%     stimMarkers = {'air','air','air','air'};
%     rasterTypes = {'dist','dist','dist','dist'};
    selCells = 'areCells';
    allVals = []; gAllVals = [];
    for ss = 1:length(stimMarkers)
        distD = [];
        for jj = 1:length(selAnimals)
            [pcs cns areCells] = getParamValues('placeCells1',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),'air','dist',selCells,maxDistTime);
            [pcs cns areCells] = getParamValues('placeCells1',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),'air','time',selCells,maxDistTime);
            pcs = logical(ones(size(pcs)));
            [tempVals cns ACs] = getParamValues('info_metrics.ShannonMI_Zsh',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
%             [tempVals cns ACs] = getParamValues('gauss_fit_on_mean.rs',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
%             [tempVals cns ACs] = getParamValues('data.rs',ei(selAnimals(jj)),planeNumbers,contextNumber,stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
%             [temppws cns ACs] = getParamValues('data.mean_trial_corr',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
%             tempVals = temppws;
%             tempVals = tempVals(logical(pcs));
            distD = [distD;tempVals(pcs)];
            zMI_mean(jj,ss) = nanmean(tempVals(pcs));
            data{jj,ss} = tempVals(pcs);
        end
        allVals{ss} = distD;
        gAllVals = [gAllVals;distD];
    end
    minBin = min(gAllVals);
    maxBin = 15;%max(gAllVals);
    incr = (maxBin-minBin)/100;
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 6.9 2.5],'color','w');
    hold on;
    [ha,hb,hca,sigR] = plotDistributions(data,'colors',colors,'maxY',90,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    hold on;
    legs = [];
    for ii = 1:length(stimMarkers)
        if ii == 1
            legs{ii} = sprintf('%s-%s (N = %d)',stimMarkers{ii},rasterTypes{ii},size(data,1));
        else
            legs{ii} = sprintf('%s-%s',stimMarkers{ii},rasterTypes{ii});
        end
    end
    ylim([0 120]);
    xlims = xlim; dx = xlims(2) - xlims(1); ylims = ylim; dy = ylims(2) - ylims(1);
    legs{ii+1} = [xlims(1)+dx/2 dx/50 ylims(1)+dy/2 dy/15];
    if size(data,2) > 2
        putLegend(ha,legs,'colors',colors,'sigR',{sigR,'anova',sigColor,8});
    else
        putLegend(ha,legs,'colors',colors,'sigR',{sigR,'ttest',sigColor,10});
    end
    axes(ha);
    h = xlabel('Mutual Information Z-Score');%changePosition(h,[0 -dy/3 0]);
    h = ylabel('Percentage');changePosition(h,[0 0 0]);
    set(gca,'FontSize',axes_font_size,'FontWeight','Bold');changePosition(ha,[-0.04 0.09 0.13 -0.05]);
    save_pdf(hf,mData.pdf_folder,'Distribution Of dfbyf',600);
    if size(data,2) > 2
        sigR.anova
    else
        sigR.ttest
    end
return;
end

