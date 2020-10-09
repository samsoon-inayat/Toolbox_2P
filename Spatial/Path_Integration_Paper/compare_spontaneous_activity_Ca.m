function compare_spontaneous_activity

protocol_C = '10_C';
protocol_A = '10_A';
ei_C = evalin('base','ei10_C');
ei_A = evalin('base','ei10_A');
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ET_C = evalin('base',sprintf('ET%s',protocol_C));
ET_A = evalin('base',sprintf('ET%s',protocol_A));
selAnimals_C = 1:length(ei_C);
selAnimals_A = 1:length(ei_A);
n = 0;
%%
out_C = get_spike_rate(ei_C);
out_A = get_spike_rate(ei_A);

n=0;
%%
if 1
    [mVar_C,semVar_C] = findMeanAndStandardError(out_C.count_prom_animal);
    [mVar_A,semVar_A] = findMeanAndStandardError(out_A.count_prom_animal);
    mVar = [mVar_C,mVar_A];
    semVar = [semVar_C,semVar_A];
    [h,p,ci,stats] = ttest2(out_C.count_prom_animal,out_A.count_prom_animal)
    combs = [1 2]; %p = ones(size(combs,1),1); h = logical(zeros(size(combs,1),1));
    xdata = [1 2]; maxY = 4;
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 2.25 1],'color','w');
    hold on;
    tcolors = {colors{1};colors{1}};
    hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'maxY',maxY,'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 2:2:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+1],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'Ctrl','AD'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.02 0.03 0.02 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{'Tra./min',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bargraph',mfilename),600);
end



function out = get_spike_rate(ei_C)
allVals = []; allVals_motion = []; allVals_rest = []; allVals_th = [];
for ii = 1:length(ei_C)
    ii
    ei = ei_C{ii};
    spSigAll = nan(length(ei.plane{1}.tP.deconv.caSigAll),length(ei.plane{1}.tP.deconv.caSigAll{1}));
    iscell = logical(ei.plane{1}.tP.iscell(:,1));
    count_prom = [];
    wp_prom = [];
    for cn = 1:length(ei.plane{1}.tP.deconv.caSigAll)
        if ~iscell(cn)
            continue;
        end
        spSigAll(cn,:) = ei.plane{1}.tP.deconv.caSigAll{cn}';
        thisSig = ei.plane{1}.tP.deconv.caSigAll{cn}';
        threshold = 3*std(thisSig);
        thisSig_th = thisSig > threshold;
        [pks,locs,widths,prom] = findpeaks(thisSig,ei.thorExp.frameRate);
        [pks,locs,widths,prom] = findpeaks(thisSig);
        mean_proms(cn) = mean(prom);
        std_proms(cn) = std(prom);
        prom_th = 35;
        indsprom = prom > prom_th;
        wp_prom(cn) = sum(pks(indsprom));
        count_prom(cn) = sum(indsprom);
%         figure(100);clf;
%         plot(thisSig);hold on;
%         plot(locs(indsprom),pks(indsprom),'*')
% %         plot(ones(size(thisSig))*threshold)
% %         plot(thisSig_th * 150)
%         n = 0;
    end
    ts = ei.b.ts(ei.plane{1}.b.frames_f);
    wp_prom = wp_prom(iscell);
    count_prom = count_prom(iscell);
    count_prom_animal(ii) = mean(count_prom)/(ts(end)/60);
    wp_prom_animal(ii) = mean(wp_prom)/(ts(end)/60);
end
out.count_prom_animal = count_prom_animal;
out.wp_prom_animal = wp_prom_animal;