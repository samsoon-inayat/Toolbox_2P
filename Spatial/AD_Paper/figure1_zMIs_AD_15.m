function figure1_number_of_PCs

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
paramMs_C = parameter_matrices('get','15_C');
paramMs_A = parameter_matrices('get','15_A');
% after getting all matrics, we can apply selection criteria to select a
% subgroup of cells
% here is the selection criteria in make_selC_structure function
% cellsOrNot = 1; planeNumber = NaN; zMI_Th = 3; fwids = [1 120]; fcens = [0 140]; rs_th = 0.4;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
conditionsAndRasterTypes = [31 41 51];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN);
[cpMs_C,pMs_C] = parameter_matrices('select','15_C',{paramMs_C,selC});
[cpMs_A,pMs_A] = parameter_matrices('select','15_A',{paramMs_A,selC});
perc_cells_C = parameter_matrices('print','15_C',{cpMs_C,pMs_C,ET_C,selAnimals_C});
perc_cells_A = parameter_matrices('print','15_A',{cpMs_A,pMs_A,ET_A,selAnimals_A});

all_conds = []; all_rts = [];
gAllVals_C = []; gAllVals_A = [];
for rr = 1:size(pMs_C,1)
    for cc = 1:size(pMs_C,2)
        tcond = conditionsAndRasterTypes(rr,cc);
        nds = dec2base(tcond,10) - '0';
        varNames{rr,cc} = sprintf('C%dR%d',nds(1),nds(2));
        all_conds = [all_conds nds(1)]; all_rts = [all_rts nds(2)];
        xticklabels{cc,rr} = sprintf('%s-%s',paramMs_C.stimMarkers{nds(2)},paramMs_C.rasterTypes{nds(2)}(1));
        for an = 1:length(selAnimals_C)
            zMIs_C(an,rr,cc) = nanmean(squeeze(pMs_C{rr,cc}.all_zMIs{selAnimals_C(an)}(nds(1)-2,nds(2),:)));
            a_zMIs_C{an,rr,cc} = squeeze(pMs_C{rr,cc}.all_zMIs{selAnimals_C(an)}(nds(1)-2,nds(2),:));
            gAllVals_C = [gAllVals_C;a_zMIs_C{an,rr,cc}];
        end
    end
end

for rr = 1:size(pMs_A,1)
    for cc = 1:size(pMs_A,2)
        tcond = conditionsAndRasterTypes(rr,cc);
        nds = dec2base(tcond,10) - '0';
        varNames{rr,cc} = sprintf('C%dR%d',nds(1),nds(2));
        all_conds = [all_conds nds(1)]; all_rts = [all_rts nds(2)];
        xticklabels{cc,rr} = sprintf('%s-%s',paramMs_A.stimMarkers{nds(2)},paramMs_A.rasterTypes{nds(2)}(1));
        for an = 1:length(selAnimals_A)
            zMIs_A(an,rr,cc) = nanmean(squeeze(pMs_A{rr,cc}.all_zMIs{selAnimals_A(an)}(nds(1)-2,nds(2),:)));
            a_zMIs_A{an,rr,cc} = squeeze(pMs_A{rr,cc}.all_zMIs{selAnimals_A(an)}(nds(1)-2,nds(2),:));
            gAllVals_A = [gAllVals_A;a_zMIs_A{an,rr,cc}];
        end
    end
end
all_conds = unique(all_conds); all_rts = unique(all_rts);
var_oi_A = squeeze(zMIs_A);
var_oi_C = squeeze(zMIs_C);
%%
runthis = 1;
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
    colVar1 = [ones(1,numCols) 2*ones(1,numCols) 3*ones(1,numCols)];    colVar2 = [1:numCols 1:numCols 1:numCols 1:numCols];
    within = table(colVar1');
    within.Properties.VariableNames = {'Condition'};
    within.Condition = categorical(within.Condition);
%     within.Raster_Type = categorical(within.Raster_Type);

    rm = fitrm(dataT,sprintf('%s,%s,%s~Group',varNames{1},varNames{2},varNames{3}),'WithinDesign',within,'WithinModel','Condition');
    rtable = ranova(rm,'WithinModel',rm.WithinModel);
    file_name = fullfile(mData.pdf_folder,sprintf('%s_RANOVA.xlsx',mfilename));
    writetable(rtable,file_name,'WriteRowNames',true)
    file_name = fullfile(mData.pdf_folder,sprintf('%s_dataT.xlsx',mfilename));
    writetable(dataT,file_name,'WriteRowNames',true)
    mauchlytbl = mauchly(rm);
    % multcompare(rm,'Day','ComparisonType','bonferroni')
    mcGroup = find_sig_mctbl(multcompare(rm,'Group','By','Condition','ComparisonType','bonferroni'),6);
    mcGroup = multcompare(rm,'Group','By','Condition','ComparisonType','bonferroni');
    file_name = fullfile(mData.pdf_folder,sprintf('%s_RANOVA_multcompare.xlsx',mfilename));
    writetable(mcGroup,file_name,'WriteRowNames',true);
    

    [mVarT,semVarT] = findMeanAndStandardError(var_oi_C);
    [mVarT_A,semVarT_A] = findMeanAndStandardError(var_oi_A);

    mVar = NaN(1,2*(length(mVarT)));
    semVar = mVar;
    mVar(1:2:6) = mVarT;semVar(1:2:6) = semVarT;
    mVar(2:2:6) = mVarT_A;semVar(2:2:6) = semVarT_A;
    combs = nchoosek(1:6,2); p = ones(size(combs,1),1); h = logical(zeros(size(combs,1),1));
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
    % row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
    % row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

    xdata = [1 2 4 5 7 8]; maxY = 3;
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 2.25 1],'color','w');
    hold on;
    tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3}};
    hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'maxY',maxY,'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',-0.3);
    for ii = 2:2:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+1],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:2:end)+0.5; xticklabels = {'C1','C2','C3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.02 0.03 0.02 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{'Average zMI',[0 0 0]});
    rectangle(gca,'Position',[0.75 3.5 1 0.5],'edgecolor','k','facecolor','k');
    text(1.85,3.5,'Control','FontSize',5);
    rectangle(gca,'Position',[6 3.5 1 0.5],'edgecolor','k');
    text(7.2,3.5,'APP','FontSize',5);
%     text(3.5,18,'Control','FontSize',7);
%     text(18.5,18,'APP','FontSize',7);
    % applyhatch_plusC(gcf
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bargraph',mfilename),600);
return;
end

%%
runthis = 1;
if runthis
    data_C = squeeze(a_zMIs_C);
    data_A = squeeze(a_zMIs_A);
%     data_A
    ii = 1; data = [data_C(:,ii) [data_A(:,ii);{NaN}]];
    minBin = min([gAllVals_C;gAllVals_A]);
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
    h = xlabel('Mutual Information Z-Score');%changePosition(h,[0 -dy/3 0]);
    h = ylabel('Percentage');changePosition(h,[-0.1 0 0]);
    set(gca,'FontSize',axes_font_size,'FontWeight','Bold');changePosition(ha,[0.07 0.01 -0.05 -0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution Of zMI %d _15',ii),600);
return;
end