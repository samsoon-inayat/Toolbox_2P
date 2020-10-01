function figure_speed_AD_Ctrl(fn,allRs,ccs)

protocol_C = '15_C';
protocol_A = '15_A';
ei_C = evalin('base','ei15_C');
ei_A = evalin('base','ei15_A');
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ET_C = evalin('base',sprintf('ET%s',protocol_C));
ET_A = evalin('base',sprintf('ET%s',protocol_A));
selAnimals_C = 1:length(ei_C)
selAnimals_A = 1:length(ei_A)
%%
[mean_speed_middle_C,mean_speed_end_C] = get_speeds(ei_C);
[mean_speed_middle_A,mean_speed_end_A] = get_speeds(ei_A);

data_CT = [ones(size(mean_speed_middle_C,1),1) mean_speed_middle_C mean_speed_end_C];
data_CT(:,2:2:end) = mean_speed_middle_C; data_CT(:,3:2:end) = mean_speed_end_C;
data_AT = [2*ones(size(mean_speed_middle_A,1),1) mean_speed_middle_A mean_speed_end_A];
data_AT(:,2:2:end) = mean_speed_middle_A; data_AT(:,3:2:end) = mean_speed_end_A;

w2 = {'Mid','End'};
varNames = []; ind = 1;
temp_tcolors = [];
for ii = 1:3
    for jj = 1:2
        varNames{ind} = sprintf('C%d%s',ii,w2{jj});
        temp_tcolors{ind} = colors{ii};
        ind = ind + 1;
    end
end

dataT = array2table([data_CT;data_AT]);
dataT.Properties.VariableNames = {'Group',varNames{:}};
numCols = size(dataT,2)-1;
dataT.Group = categorical(int32(dataT.Group))
colVar1 = [1 1 2 2 3 3];    colVar2 = [1 2 1 2 1 2];
within = table(colVar1',colVar2');
within.Properties.VariableNames = {'Condition','Bin'};
within.Condition = categorical(within.Condition);
within.Bin = categorical(within.Bin);

rm = fitrm(dataT,sprintf('%s,%s,%s,%s,%s,%s~Group',varNames{1},varNames{2},varNames{3},varNames{4},...
    varNames{5},varNames{6}),'WithinDesign',within,'WithinModel','Condition*Bin');
rtable = ranova(rm,'WithinModel',rm.WithinModel);
file_name = fullfile(mData.pdf_folder,sprintf('%s_Data.xlsx',mfilename));
    writetable(dataT,file_name,'WriteRowNames',true)
file_name = fullfile(mData.pdf_folder,sprintf('%s_RANOVA.xlsx',mfilename));
writetable(rtable,file_name,'WriteRowNames',true)
mauchlytbl = mauchly(rm);
% multcompare(rm,'Day','ComparisonType','bonferroni')
mcGroup = find_sig_mctbl(multcompare(rm,'Group','By','Condition','ComparisonType','bonferroni'),6);
mcGroup = multcompare(rm,'Group','By','Condition','ComparisonType','bonferroni');
file_name = fullfile(mData.pdf_folder,sprintf('%s_RANOVA_multcompare.xlsx',mfilename));
writetable(mcGroup,file_name,'WriteRowNames',true)

    [mVarT,semVarT] = findMeanAndStandardError(data_CT(:,2:end));
    [mVarT_A,semVarT_A] = findMeanAndStandardError(data_AT(:,2:end));
    mVar = [mVarT mVarT_A];
    semVar = [semVarT semVarT_A];
    tcolors = num2cell(mVar);
%     mVar = mVarT;semVar(1:2:TL) = semVarT;
%     mVar(2:2:TL) = mVarT_A;semVar(2:2:TL) = semVarT_A;
    TL = 12;
    tcolors(1:6) = temp_tcolors; tcolors(7:TL) = temp_tcolors;
    combs = nchoosek(1:TL,2); p = ones(size(combs,1),1); h = logical(zeros(size(combs,1),1));
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
    % row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
    % row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

    xdata = [1:6 (7:TL)+1]; maxY = 30;
    colors = mData.colors;
%     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 6.9 2],'color','w');
    hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 6.9 2],'color','w');
    hold on;
    ind = 1

    hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'maxY',maxY,'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',-0.3);
    for ii = 7:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+1],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:2:end)+0.5; 
    
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.02 0.03 0.02 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{'Speed (cm/s)',[0 0 0]});
%     rectangle(gca,'Position',[0.75 9 1 2],'edgecolor','k','facecolor','k');
%     text(1.85,10,'Control','FontSize',5);
%     rectangle(gca,'Position',[6 9 1 2],'edgecolor','k');
%     text(7.2,10,'APP','FontSize',5);
%     text(3.5,18,'Control','FontSize',7);
%     text(18.5,18,'APP','FontSize',7);
    % applyhatch_plusC(gcf

save_pdf(hf,mData.pdf_folder,sprintf('%s_speed',mfilename),600);

    
    function [mean_speed_middle,mean_speed_end] = get_speeds(ei_C)
for an = 1:length(ei_C)
    temp_bl = [];
    for cn = 1:3
        try
            onsets = ei_C{an}.plane{1}.contexts(cn).markers.air_onsets;
            offsets = ei_C{an}.plane{1}.contexts(cn).markers.air_offsets;
        catch
            onsets = ei_C{an}.plane{1}.contexts(cn+2).markers.air_onsets;
            offsets = ei_C{an}.plane{1}.contexts(cn+2).markers.air_offsets;
        end
        b = ei_C{an}.b;
        temp_bl(cn,:) = b.dist(offsets)-b.dist(onsets);
        m_belt_length(an,cn) = mean(b.dist(offsets)-b.dist(onsets));
        try
            sz(an,cn) = size(ei_C{an}.plane{1}.contexts(cn).rasters.airD.sp_rasters_nan_corrected,2);
            speeds_C{an,cn} = ei_C{an}.plane{1}.contexts(cn).rasters.airD.speed(:,1:sz(an,cn));
        catch
            sz(an,cn) = size(ei_C{an}.plane{1}.contexts(cn+2).rasters.airD.sp_rasters_nan_corrected,2);
            speeds_C{an,cn} = ei_C{an}.plane{1}.contexts(cn+2).rasters.airD.speed(:,1:sz(an,cn));
        end
        this_speed_raster = speeds_C{an,cn};
        middle_bin = floor(sz(an,cn)/2);
        start_bin = middle_bin - 2;
        end_bin = middle_bin + 2;
        middle_speeds = this_speed_raster(:,start_bin:end_bin);
        mean_speed_middle(an,cn) = nanmean(middle_speeds(:));
        end_speeds = this_speed_raster(:,(sz(an,cn)-4):end);
        mean_speed_end(an,cn) = nanmean(end_speeds(:));
        
    end
    belt_lengths{an} = temp_bl;
end