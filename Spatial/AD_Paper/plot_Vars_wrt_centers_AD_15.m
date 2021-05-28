function plot_Vars_wrt_centers_AD_15

protocol_C = '15_C';
protocol_A = '15_A';
ei_C = evalin('base','ei15_C');
ei_A = evalin('base','ei15_A');
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ET_C = evalin('base',sprintf('ET%s',protocol_C));
ET_A = evalin('base',sprintf('ET%s',protocol_A));
selAnimals_C = 1:length(ei_C)
selAnimals_A = 1:length(ei_A)

% in the following variable all the measurements are in the matrices form
% for each variable colums indicate raster and stim marker types specified 
% the rows indicate condition numbers.
paramMs_C = parameter_matrices('get',protocol_C);
paramMs_A = parameter_matrices('get',protocol_A);
paramMs_C.belt_lengths = get_mean_belt_length(ei_C,protocol_C)
paramMs_A.belt_lengths = get_mean_belt_length(ei_A,protocol_A)
% after getting all matrics, we can apply selection criteria to select a
% subgroup of cells
% here is the selection criteria in make_selC_structure function
all_variables = {'all_zMIs','all_fFR','all_fwidths','all_frs',''};
ylabels = {'zMIs','Firing Rate','PF Widths','RS','Percent PCs'};
maxYs = [10,40,50,0.7,100];
svn = 5; %gcn = 3
if svn == 5
    selected_variable = all_variables{1};
    selected_variable_f = 'Percent_PCs';
    number_of_bins = 4;
else
    selected_variable = all_variables{svn};
    selected_variable_f = selected_variable;
    number_of_bins = 1;
end
cellsOrNot = 1; planeNumber = NaN; zMI_Th = 3; fwids = [1 120]; fcens = [0 140]; rs_th = 0.4;
% cellsOrNot = 1; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
conditionsAndRasterTypes = [31 41 51];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN);
[cpMs_C,pMs_C] = parameter_matrices('select',protocol_C,{paramMs_C,selC});
[cpMs_A,pMs_A] = parameter_matrices('select',protocol_A,{paramMs_A,selC});
% perc_cells_C = parameter_matrices('print','10_C',{cpMs_C,pMs_C,ET_C,selAnimals_C});
% perc_cells_A = parameter_matrices('print','10_A',{cpMs_A,pMs_A,ET_A,selAnimals_A});

out_C = get_values(paramMs_C,pMs_C,selAnimals_C,selected_variable,conditionsAndRasterTypes,number_of_bins);
out_A = get_values(paramMs_A,pMs_A,selAnimals_A,selected_variable,conditionsAndRasterTypes,number_of_bins);

all_conds = unique(out_C.all_conds); all_rts = unique(out_C.all_rts);
var_oi_A = squeeze(out_A.sV);
var_oi_C = squeeze(out_C.sV);
dist_oi_A = squeeze(out_A.a_sV);
dist_oi_C = squeeze(out_C.a_sV);
varNames = out_C.varNames;
%% for making bar graphs from means over cells (for individual animals)
runthis = 0;
if runthis
    numCols = length(all_rts);
    data = var_oi_C;
    cmdTxt = sprintf('dataT_C = table(');
    for ii = 1:(size(data,2)-1)
        cmdTxt = sprintf('%sdata(:,%d),',cmdTxt,ii);
    end
    cmdTxt = sprintf('%sdata(:,size(data,2)));',cmdTxt);
    eval(cmdTxt);
    data = var_oi_A;
    cmdTxt = sprintf('dataT_A = table(');
    for ii = 1:(size(data,2)-1)
        cmdTxt = sprintf('%sdata(:,%d),',cmdTxt,ii);
    end
    cmdTxt = sprintf('%sdata(:,size(data,2)));',cmdTxt);
    eval(cmdTxt);
    dataT = [dataT_C;dataT_A]
    dataT.Properties.VariableNames = varNames;
    dataT = [table([ones(length(ei_C),1);2*ones(length(ei_A),1)]) dataT];
    dataT.Properties.VariableNames{1} = 'Group';
    dataT.Group = categorical(dataT.Group)
    colVar1 = [ones(1,numCols) 2*ones(1,numCols) 3*ones(1,numCols) 4*ones(1,numCols)];    colVar2 = [1:numCols 1:numCols 1:numCols 1:numCols];
    within = table(colVar1');
    within.Properties.VariableNames = {'Condition'};
    within.Condition = categorical(within.Condition);
%     within.Raster_Type = categorical(within.Raster_Type);
    file_name = fullfile(mData.pdf_folder,sprintf('%s_Data.xlsx',mfilename));
    writetable(dataT,file_name,'WriteRowNames',true)

    rm = fitrm(dataT,sprintf('%s,%s,%s,%s~Group',varNames{1},varNames{2},varNames{3},varNames{4}),'WithinDesign',within,'WithinModel','Condition');
    rtable = ranova(rm,'WithinModel',rm.WithinModel);
%     writetable(dataT,'Percentage_of_PCs_10.xlsx','WriteRowNames',true)
    file_name = fullfile(mData.pdf_folder,sprintf('%s_RANOVA.xlsx',mfilename));
    writetable(rtable,file_name,'WriteRowNames',true)
    mauchlytbl = mauchly(rm);
    % multcompare(rm,'Day','ComparisonType','bonferroni')
    mcGroup = find_sig_mctbl(multcompare(rm,'Group','By','Condition','ComparisonType','bonferroni'),6);
    mcGroup = multcompare(rm,'Group','By','Condition','ComparisonType','bonferroni');
    file_name = fullfile(mData.pdf_folder,sprintf('%s_RANOVA_multcompare.xlsx',mfilename));
    writetable(mcGroup,file_name,'WriteRowNames',true)
%     mcTI = find_sig_mctbl(multcompare(rm,'TI','By','Day','ComparisonType','bonferroni'),6);
%     mcDays = find_sig_mctbl(multcompare(rm,'Day','By','TI','ComparisonType','bonferroni'),6);
    

    [mVarT,semVarT] = findMeanAndStandardError(var_oi_C);
    [mVarT_A,semVarT_A] = findMeanAndStandardError(var_oi_A);

    mVar = NaN(1,2*(length(mVarT)));
    semVar = mVar;
    mVar(1:2:8) = mVarT;semVar(1:2:8) = semVarT;
    mVar(2:2:8) = mVarT_A;semVar(2:2:8) = semVarT_A;
    combs = nchoosek(1:8,2); p = ones(size(combs,1),1); h = logical(zeros(size(combs,1),1));
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
    % row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
    % row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

    xdata = [1 2 4 5 7 8 10 11]; maxY = maxYs(svn);
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 2.25 1],'color','w');
    hold on;
    tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
    hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'maxY',maxY,'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',-0.3);
    for ii = 2:2:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+1],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:2:end)+0.5; xticklabels = {'C1','C2','C3','C4'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.02 0.03 0.02 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{ylabels{svn},[0 0 0]});
    rectangle(gca,'Position',[0.75 9 1 2],'edgecolor','k','facecolor','k');
    text(1.85,10,'Control','FontSize',5);
    rectangle(gca,'Position',[6 9 1 2],'edgecolor','k');
    text(7.2,10,'APP','FontSize',5);
%     text(3.5,18,'Control','FontSize',7);
%     text(18.5,18,'APP','FontSize',7);
    % applyhatch_plusC(gcf
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bargraph',mfilename),600);
return;
end

%% overall average distributions for the two groups  (not w.r.t centers)
runthis = 0;
if runthis
    data_C = squeeze(dist_oi_C);
    data_A = squeeze(dist_oi_A);
    data = [data_C(:,gcn) data_A(:,gcn)];
    minBin = min([out_C.gAllVals;out_A.gAllVals]);
    maxBin = 15;%max(gAllVals);
    incr = (maxBin-minBin)/100;
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 6 2 1.5],'color','w');
    hold on;
    [ha,hb,hca,sigR] = plotDistributions(data,'colors',colors,'maxY',90,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    hold on;
    legs = {sprintf('%s-C',varNames{ii}),sprintf('%s-A',varNames{ii})};
%     paramMs = paramMs_A;
%     for ii = 1:length(varNames)
%         if ii == 1
%             legs{ii} = sprintf('%s-%s (N = %d)',paramMs.stimMarkers{ii},paramMs.rasterTypes{ii}(1),size(data,1));
%         else
%             legs{ii} = sprintf('%s-%s',paramMs.stimMarkers{ii},paramMs.rasterTypes{ii}(1));
%         end
%     end
    ylim([0 110]);
    xlim([-3 12])
    xlims = xlim; dx = xlims(2) - xlims(1); ylims = ylim; dy = ylims(2) - ylims(1);
    legs{length(legs)+1} = [xlims(1)+dx/1.5 dx/30 ylims(1)+dy/3 dy/15];
    putLegend(ha,legs,'colors',colors);
    axes(ha);
    h = xlabel(ylabels{svn});%changePosition(h,[0 -dy/3 0]);
    h = ylabel('Percentage');changePosition(h,[-0.1 0 0]);
    set(gca,'FontSize',axes_font_size,'FontWeight','Bold');changePosition(ha,[0.07 0.01 -0.05 -0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution Of zMI %d',ii),600);
return;
end

%% average distributions w.r.t centers for the two groups
runthis = 1;
if runthis
%     cN = 1;
    data_C = [];
    data_A = [];
    ind = 1;
    for gcni = 1:3
        if svn == 5
            tempC = squeeze(out_C.all_data_N(gcni,:,:))';
            tempA = squeeze(out_A.all_data_N(gcni,:,:))';
            if size(tempC,1) == 1
                tempC = tempC'; tempA = tempA';
            end
            sum_C = sum(tempC,2);
            sum_A = sum(tempA,2);
            for ii = 1:size(tempC,1)
                tempC(ii,:) = 100*tempC(ii,:)/sum_C(ii);
            end
            for ii = 1:size(tempA,1)
                tempA(ii,:) = 100*tempA(ii,:)/sum_A(ii);
            end
        else
            tempC = squeeze(out_C.all_data(gcni,:,:))';
            tempA = squeeze(out_A.all_data(gcni,:,:))';
            if size(tempC,1) == 1
                tempC = tempC'; tempA = tempA';
            end
        end
        
        data_C = [data_C tempC];
        data_A = [data_A tempA];
        if size(tempC,2) > 1
            for ii = 1:size(tempC,2)
                varNames{ind} = sprintf('C%dB%d',gcni,ii);
                xticklabels{ind} = sprintf('B%d',ii);
                temp_tcolors{ind} = colors{gcni};
                ind = ind + 1;
            end
        else
            varNames{ind} = sprintf('C%d',gcni); xticklabels{ind} = varNames{ind};
            temp_tcolors{ind} = colors{gcni};ind = ind + 1;
        end
    end
    data_CT = [ones(size(data_C,1),1) data_C];
    data_AT = [2*ones(size(data_A,1),1) data_A];
    
    dataT = array2table([data_CT;data_AT]);
    dataT.Properties.VariableNames = {'Group',varNames{:}};
    numCols = size(tempC,2);
    dataT.Group = categorical(int32(dataT.Group))
    colVar1 = [ones(1,numCols) 2*ones(1,numCols) 3*ones(1,numCols)];    colVar2 = [1:numCols 1:numCols 1:numCols];
    within = table(colVar1',colVar2');
    within.Properties.VariableNames = {'Condition','Bin'};
    within.Condition = categorical(within.Condition);
    within.Bin = categorical(within.Bin);
%     rmaR = repeatedMeasuresAnova(dataT,within);
%     rm = fitrm(dataT,sprintf('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s~Group',...
%         varNames{1},varNames{2},varNames{3},varNames{4},varNames{5},varNames{6},varNames{7},varNames{8},...
%         varNames{9},varNames{10},varNames{11},varNames{12},varNames{13},varNames{14},varNames{15},varNames{16}),'WithinDesign',within,'WithinModel','Condition*Bin');
%     rtable = ranova(rm,'WithinModel',rm.WithinModel);
%     file_name = fullfile(mData.pdf_folder,sprintf('%s_RANOVA_%s.xlsx',mfilename,selected_variable_f));
%     writetable(rtable,file_name,'WriteRowNames',true);
    file_name = fullfile(mData.pdf_folder,sprintf('%s_Data_%s.xlsx',mfilename,selected_variable_f));
    writetable(dataT,file_name,'WriteRowNames',true)
%     mauchlytbl = mauchly(rm);
%     % multcompare(rm,'Day','ComparisonType','bonferroni')
%     mcGbyC = find_sig_mctbl(multcompare(rm,'Group','By','Condition','ComparisonType','bonferroni'),6);
%     mcGbyB = find_sig_mctbl(multcompare(rm,'Group','By','Bin','ComparisonType','bonferroni'));
% %     writetable(mcGroup,'RANOVA_Number_of_PCs_multcompare_10.xlsx')
    
    
    [mVarT,semVarT] = findMeanAndStandardError(data_C);
    [mVarT_A,semVarT_A] = findMeanAndStandardError(data_A);
    mVar = NaN(1,2*(length(mVarT)));
    semVar = mVar; TL = 2*length(mVarT);
    tcolors = num2cell(mVar);
    mVar(1:2:TL) = mVarT;semVar(1:2:TL) = semVarT;
    mVar(2:2:TL) = mVarT_A;semVar(2:2:TL) = semVarT_A;
    tcolors(1:2:TL) = temp_tcolors; tcolors(2:2:TL) = temp_tcolors;
    combs = nchoosek(1:TL,2); p = ones(size(combs,1),1); h = logical(zeros(size(combs,1),1));
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
    % row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
    % row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

    xdata = 1:TL; maxY = maxYs(svn);
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 6.9 2],'color','w');
%     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 2.25 1],'color','w');
    hold on;
    ind = 1

    hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'maxY',maxY,'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',-0.3);
    for ii = 2:2:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+1],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:2:end)+0.5; 
    
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.02 0.03 0.02 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{ylabels{svn},[0 0 0]});
%     rectangle(gca,'Position',[0.75 9 1 2],'edgecolor','k','facecolor','k');
%     text(1.85,10,'Control','FontSize',5);
%     rectangle(gca,'Position',[6 9 1 2],'edgecolor','k');
%     text(7.2,10,'APP','FontSize',5);
%     text(3.5,18,'Control','FontSize',7);
%     text(18.5,18,'APP','FontSize',7);
    % applyhatch_plusC(gcf
    if svn == 5
        save_pdf(hf,mData.pdf_folder,sprintf('%s_distributions_PercentPCs_centers',mfilename),600);
    else
        save_pdf(hf,mData.pdf_folder,sprintf('%s_distributions_%s_centers',mfilename,selected_variable_f),600);
    end
return;
end


function out = get_values(paramMs_C,pMs_C,selAnimals_C,selected_variable,conditionsAndRasterTypes,number_of_bins)
all_conds = []; all_rts = []; gAllVals_C = [];
for rr = 1:size(pMs_C,1)
    for cc = 1:size(pMs_C,2)
        tcond = conditionsAndRasterTypes(rr,cc);
        nds = dec2base(tcond,10) - '0';
        varNames{rr,cc} = sprintf('C%dR%d',nds(1),nds(2));
        all_conds = [all_conds nds(1)]; all_rts = [all_rts nds(2)];
        xticklabels{cc,rr} = sprintf('%s-%s',paramMs_C.stimMarkers{nds(2)},paramMs_C.rasterTypes{nds(2)}(1));
        for an = 1:length(selAnimals_C)
            f_centers_C{an,rr,cc} = squeeze(pMs_C{rr,cc}.all_fcenters{selAnimals_C(an)}(nds(1)-2,nds(2),:));
            cmdTxt = sprintf('sV_C(an,rr,cc) = nanmean(squeeze(pMs_C{rr,cc}.%s{selAnimals_C(an)}(nds(1)-2,nds(2),:)));',selected_variable);
            eval(cmdTxt)
            cmdTxt = sprintf('a_sV_C{an,rr,cc} = squeeze(pMs_C{rr,cc}.%s{selAnimals_C(an)}(nds(1)-2,nds(2),:));',selected_variable);
            eval(cmdTxt)
            gAllVals_C = [gAllVals_C;a_sV_C{an,rr,cc}];
        end
    end
end


bins = 0:number_of_bins:155;
for an = 1:length(selAnimals_C)
    for cc = 1:length(conditionsAndRasterTypes)
        theseCenters = f_centers_C{an,1,cc};
        these_sV_Vals = a_sV_C{an,1,cc};
         mbl = paramMs_C.belt_lengths{an}(cc);
        binSize = mbl/number_of_bins;
        bins = 0:binSize:mbl;
        [N,E,Bi] = histcounts(theseCenters,bins);
        mean_sv_Vals = [];
        for bb = 1:length(N)
            mean_sv_Vals(bb) = nanmean(these_sV_Vals(Bi == bb));
%             mean_sv_Vals(bb) = nanmedian(these_sV_Vals(Bi == bb));
        end
%         mean_sv_Vals(isnan(mean_sv_Vals)) = 0;
        all_data_C(cc,:,an) = mean_sv_Vals;
        all_data_N(cc,:,an) = N;
    end
end
out.varNames = varNames;
out.all_conds = all_conds;
out.all_rts = all_rts;
out.gAllVals = gAllVals_C;
out.sV = sV_C;
out.a_sV = a_sV_C;
out.all_data = all_data_C;
out.all_data_N = all_data_N;
