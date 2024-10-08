function figure_speed_AD_Ctrl(fn,allRs,ccs)

ei_C = evalin('base','ei10_C');
ei_A = evalin('base','ei10_A');
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;

selContexts1 = [1 2 3 4 1 2 3 4 ];
rasterNames1 = {'airD','airD','airD','airD','airIT','airIT','airIT','airIT'};
[~,ei_C] = get_rasters_data(ei_C,selContexts1,rasterNames1); 
[~,ei_A] = get_rasters_data(ei_A,selContexts1,rasterNames1); 

out_C = get_speeds(ei_C);
out_A = get_speeds(ei_A);
n = 0;
%% R1 (during revision)
ms_C = []; ms_A = [];
for an = 1:5
    mspeed = out_C.mspeed(:,:,an);
    mspeed = reshape(mspeed,1,40);
    ms_C = [ms_C;mspeed];
    
    mspeed = out_A.mspeed(:,:,an);
    mspeed = reshape(mspeed,1,40);
    ms_A = [ms_A;mspeed];
end

var_C = ms_C;
var_A = ms_A;
[within,dvn,xlabels] = make_within_table({'Cond','Trials'},[4,10]);
dataT = make_between_table({var_C;var_A},dvn);
ra1 = RMA(dataT,within,{0.05,{}});
% ra = RMA(dataT,within,{0.05,{'hsd'}});
ra1.ranova

print_for_manuscript(ra1) % same result as below ... confirmation with a different variable


%%
tic
for an = 1:5
    for cn = 1:4
        [an cn]
        for tn = 1:10
            mspeed_C(tn,cn,an) = mean(out_C.fspeeds{an,cn,tn});
            mspeed_A(tn,cn,an) = mean(out_A.fspeeds{an,cn,tn});
            mspeed_CI(tn,cn,an) = mean(out_C.fspeedsI{an,cn,tn});
            mspeed_AI(tn,cn,an) = mean(out_A.fspeedsI{an,cn,tn});
        end
    end
end
toc

rmspeed_C = mspeed_C;
rmspeed_A = mspeed_A;

var_C = reshape(rmspeed_C,40,5)';
var_A = reshape(rmspeed_A,40,5)';
[within,dvn,xlabels] = make_within_table({'Cond','Trials'},[4,10]);
dataT = make_between_table({var_C;var_A},dvn);
ra = RMA(dataT,within,{0.05,{}});
% ra = RMA(dataT,within,{0.05,{'hsd'}});
ra.ranova

print_for_manuscript(ra)
%%
Kmax = 20;
tic
for an = 1:5
    for cn = 1:4
        [an cn]
        parfor tn = 1:10
            t_fspeed = out_C.fspeeds{an,cn,tn};
            hifd_Ct(an,cn,tn) = Higuchi_FD(zscore(t_fspeed),Kmax);
            t_fspeed = out_A.fspeeds{an,cn,tn};
            hifd_At(an,cn,tn) = Higuchi_FD(zscore(t_fspeed),Kmax);
        end
    end
end
toc
n = 0;
%%
var_C = reshape(permute(hifd_Ct,[3 2 1]),40,5)';
var_A = reshape(permute(hifd_At,[3 2 1]),40,5)';
[within,dvn,xlabels] = make_within_table({'Cond','Trials'},[4,10]);
dataT = make_between_table({var_C;var_A},dvn);
ra = RMA(dataT,within,{0.05,{}});
ra.ranova
%%
Kmax = 20;
tic
for an = 1:5
    for cn = 1:4
        [an cn]
        t_fspeed = out_C.fspeeds_cn{an,cn};
        hifd_C(an,cn) = Higuchi_FD(zscore(t_fspeed),Kmax);
        t_fspeed = out_A.fspeeds_cn{an,cn};
        hifd_A(an,cn) = Higuchi_FD(zscore(t_fspeed),Kmax);
    end
end
toc
n = 0;
[within,dvn,xlabels] = make_within_table({'Cond'},[4]);
dataT = make_between_table({hifd_C;hifd_A},dvn);
ra = RMA(dataT,within,{0.05,{}});
ra.ranova

%%
Kmax = 20;
tic
for an = 1:5
        [an]
        t_fspeed = out_C.fspeeds_an{an};
        hifd_Can(an) = Higuchi_FD(zscore(t_fspeed),Kmax);
        t_fspeed = out_A.fspeeds_an{an};
        hifd_Aan(an) = Higuchi_FD(zscore(t_fspeed),Kmax);
end
toc
n = 0;
[within,dvn,xlabels] = make_within_table({'Cond'},[4]);
dataT = make_between_table({hifd_C;hifd_A},dvn);
ra = RMA(dataT,within,{0.05,{}});
ra.ranova


%%
Kmax = 20;
tic
for an = 1:5
    for cn = 1:4
        [an cn]
        t_fspeed = out_C.fspeeds_cn{an,cn};
        hifd_C(an,cn) = mean(t_fspeed);
        t_fspeed = out_A.fspeeds_cn{an,cn};
        hifd_A(an,cn) = mean(t_fspeed);
    end
end
toc
n = 0;
%%
figure(1000);clf;
for scn = 1:4
    for an = 1:5
        for cn = scn
            t_fspeed_C(an,:) = (out_C.fspeeds_cn{an,cn});
            t_fspeed_A(an,:) = (out_A.fspeeds_cn{an,cn});
        end
    end
    ax(scn) = subplot(4,1,scn);
    mspeed_C = mean(t_fspeed_C);
    mspeed_A = mean(t_fspeed_A);
    semspeed_C = std(t_fspeed_C)/sqrt(5);
    semspeed_A = std(t_fspeed_A)/sqrt(5);
    shadedErrorBar(1:length(mspeed_C),mspeed_C,semspeed_C,'b'); hold on;
    shadedErrorBar(1:length(mspeed_A),mspeed_A,semspeed_A,'r');
end
linkaxes(ax);
%%
while 0
    msr_C = mean(out_C.mspeed_raster,3);
    semsr_C = std(out_C.mspeed_raster,[],3)/sqrt(5);

    msr_A = mean(out_A.mspeed_raster,3);
    semsr_A = std(out_A.mspeed_raster,[],3)/sqrt(5);

    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[10 3 6.9 1.25],'RowsCols',[1 4],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
    varT = 1;
    switch varT
        case 1 % responsive cells 
            MY = 1.6; ysp = 3; mY = 1.1; titletxt = 'Speed Correlation'; ylabeltxt = {'Corr'};
        
    end
    stp = 0.25*magfac; widths = ([1.2 1.3 1.3 1.3 1.3 0.5 0.5 0.5]+0.25)*magfac; gap = 0.16*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{''});
    xs = 1:size(msr_C,1);
    axes(ff.h_axes(1,1));
    shadedErrorBar(xs,msr_C(:,1),semsr_C(:,1)); hold on;
    plot(xs,msr_A(:,1),'r'); shadedErrorBar(xs,msr_A(:,1),semsr_A(:,1));
    %%
    break;
end
%%
while 1
    var_C = reshape(out_C.hifd,40,5)';
    var_A = reshape(out_A.hifd,40,5)';
    
    [within,dvn,xlabels,withinD] = make_within_table({'Cond','Trials'},[4,10]); awithinD = withinD;
    dataT = make_between_table({var_C;var_A},dvn);
    ra = RMA(dataT,within,{0.05,{}});
    ra.ranova
    n = 0;
    %%
    cn = 1;
    [within,dvn,xlabels,withinD] = make_within_table({'Trials'},[10]);
    dataT = make_between_table({var_C(:,awithinD(:,1)==cn);var_A(:,awithinD(:,1)==cn)},dvn);
    ra1 = RMA(dataT,within,{0.05,{'hsd'}});
    ra1.ranova
    n = 0;
    %%
    var_C = squeeze(mean(out_C.mspeedI,1))';
    var_A = squeeze(mean(out_A.mspeedI,1))';
    [within,dvn,xlabels,withinD] = make_within_table({'Cond'},[4]);
    dataT = make_between_table({var_C;var_A},dvn);
    ra1 = RMA(dataT,within,{0.05,{'hsd'}});
    ra1.ranova
    n = 0;
    %%
    break;
end

%% one graph
while 1
    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[10 3 2.5 1.25],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
    varT = 1;
    switch varT
        case 1 % responsive cells 
            MY = 1.6; ysp = 3; mY = 1.1; titletxt = 'Speed Correlation'; ylabeltxt = {'Corr'};
        
    end
    stp = 0.25*magfac; widths = ([1.95 1.3 1.3 1.3 1.3 0.5 0.5 0.5]+0.25)*magfac; gap = 0.16*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{''});
    tcolors = repmat(dcolors,2,1);

    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra1,{'Group_by_Trials','hsd'},[1.5 1 1]);
        xdata = make_xdata([10 10],[1 1.5]);   
%         combs = [[1:2:12]' [2:2:12]']; p = ra.MC.bonferroni.Group_by_Cond{1:2:12,6}; h = p<0.05;
    h(h==1) = 0;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
    xticklabels = {'C1','C2','C3','C4'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(0);
    make_bars_hollow(hbs(11:end));
    put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
%     set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Control','APP'});
    ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.051 0 0]);
    format_axes_b(gca);
    set(ht,'FontWeight','Bold');
    save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
    %%
    break;
end


%%
data_CT = [ones(size(mean_speed_middle_C,1),1) mean_speed_middle_C mean_speed_end_C];
data_CT(:,2:2:end) = mean_speed_middle_C; data_CT(:,3:2:end) = mean_speed_end_C;
data_AT = [2*ones(size(mean_speed_middle_A,1),1) mean_speed_middle_A mean_speed_end_A];
data_AT(:,2:2:end) = mean_speed_middle_A; data_AT(:,3:2:end) = mean_speed_end_A;

w2 = {'Mid','End'};
varNames = []; ind = 1;
temp_tcolors = [];
for ii = 1:4
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
colVar1 = [1 1 2 2 3 3 4 4];    colVar2 = [1 2 1 2 1 2 1 2];
within = table(colVar1',colVar2');
within.Properties.VariableNames = {'Condition','Bin'};
within.Condition = categorical(within.Condition);
within.Bin = categorical(within.Bin);

rm = fitrm(dataT,sprintf('%s,%s,%s,%s,%s,%s,%s,%s~Group',varNames{1},varNames{2},varNames{3},varNames{4},...
    varNames{5},varNames{6},varNames{7},varNames{8}),'WithinDesign',within,'WithinModel','Condition*Bin');
rtable = ranova(rm,'WithinModel',rm.WithinModel);
file_name = fullfile(mData.pdf_folder,sprintf('%s_Data.xlsx',mfilename));
%     writetable(dataT,file_name,'WriteRowNames',true)
file_name = fullfile(mData.pdf_folder,sprintf('%s_RANOVA.xlsx',mfilename));
% writetable(rtable,file_name,'WriteRowNames',true)
mauchlytbl = mauchly(rm);
% multcompare(rm,'Day','ComparisonType','bonferroni')
mcGroup = find_sig_mctbl(multcompare(rm,'Group','By','Condition','ComparisonType','bonferroni'),6);
mcGroup = multcompare(rm,'Group','By','Condition','ComparisonType','bonferroni');
file_name = fullfile(mData.pdf_folder,sprintf('%s_RANOVA_multcompare.xlsx',mfilename));
% writetable(mcGroup,file_name,'WriteRowNames',true)

    [mVarT,semVarT] = findMeanAndStandardError(data_CT(:,2:end));
    [mVarT_A,semVarT_A] = findMeanAndStandardError(data_AT(:,2:end));
    mVar = [mVarT mVarT_A];
    semVar = [semVarT semVarT_A];
    tcolors = num2cell(mVar);
%     mVar = mVarT;semVar(1:2:TL) = semVarT;
%     mVar(2:2:TL) = mVarT_A;semVar(2:2:TL) = semVarT_A;
    TL = 16;
    tcolors(1:8) = temp_tcolors; tcolors(9:TL) = temp_tcolors;
    combs = nchoosek(1:TL,2); p = ones(size(combs,1),1); h = logical(zeros(size(combs,1),1));
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
    % row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
    % row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

    xdata = [1:8 (9:TL)+1]; maxY = 30;
    colors = mData.colors;
%     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 6.9 2],'color','w');
    hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 6.9 2],'color','w');
    hold on;
    ind = 1

    hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'maxY',maxY,'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',-0.3);
    for ii = 9:length(hbs)
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

   
function out = get_speeds(ei_C)
    for an = 1:length(ei_C)
        temp_bl = [];
        b = ei_C{an}.b;
        for cn = 1:4
%             [an cn]
            onsets = ei_C{an}.plane{1}.contexts(cn).markers.air_onsets;
            offsets = ei_C{an}.plane{1}.contexts(cn).markers.air_offsets;
            onsetsI = ei_C{an}.plane{1}.contexts(cn).markers.airI_onsets;
            offsetsI = ei_C{an}.plane{1}.contexts(cn).markers.airI_offsets;
            temp_bl(cn,:) = b.dist(offsets)-b.dist(onsets);
            m_belt_length(an,cn) = mean(b.dist(offsets)-b.dist(onsets));
            sz(an,cn) = size(ei_C{an}.plane{1}.contexts(cn).rasters.airD.sp_rasters_nan_corrected,2);
            speeds_C{an,cn} = ei_C{an}.plane{1}.contexts(cn).rasters.airD.speed(:,1:sz(an,cn));
            this_speed_raster = speeds_C{an,cn};
            mspeed_raster(:,cn,an) = mean(this_speed_raster(:,1:35));
            c_s = corr(this_speed_raster');
            mask = ones(size(c_s)); mask = triu(mask,1) & ~triu(mask,2);
            c_s_tp(:,cn,an) = c_s(mask);
            middle_bin = floor(sz(an,cn)/2);
            start_bin = middle_bin - 2;
            end_bin = middle_bin + 2;
            middle_speeds = this_speed_raster(:,start_bin:end_bin);
            mean_speed_middle(an,cn) = nanmean(middle_speeds(:));
            end_speeds = this_speed_raster(:,(sz(an,cn)-4):end);
            mean_speed_end(an,cn) = nanmean(end_speeds(:));
            seconds_before = 3; seconds_after = 5; 
            tfspeed = [];
            for tn = 1:10
                st = onsets(tn) - round(1e6 * seconds_before/b.si); se = onsets(tn) - round(1e6 * 0.01/b.si);
                speed_sig = b.fSpeed(st:se); ms_b(tn,cn,an) = mean(speed_sig);
                se = onsets(tn) + round(1e6 * seconds_after/b.si); st = onsets(tn) + round(1e6 * 0.01/b.si);
                speed_sig = b.fSpeed(st:se); ms_a(tn,cn,an) = mean(speed_sig);
                mspeed(tn,cn,an) = mean(b.fSpeed(onsets(tn):offsets(tn)));
                mspeedI(tn,cn,an) = mean(b.fSpeed(onsetsI(tn):offsetsI(tn)));
                son(tn,cn,an) = b.speed(onsets(tn));
                soff(tn,cn,an) = b.speed(offsets(tn));
                movL(tn,cn,an) = find(b.fSpeed(onsets(tn):offsets(tn)) > 0.1,1,'first')*b.si*1e-6;
                durT(tn,cn,an) = b.ts(offsets(tn))-b.ts(onsets(tn));
                durIT(tn,cn,an) = b.ts(offsetsI(tn))-b.ts(onsetsI(tn));
                below7(tn,cn,an) = sum(b.fSpeed(onsets(tn):offsets(tn)) < 7)*b.si*1e-6;
                above7(tn,cn,an) = sum(b.fSpeed(onsets(tn):offsets(tn)) > 7)*b.si*1e-6;
                t_speed = b.fSpeed(onsets(tn):offsets(tn));
%                 ts = b.ts(onsets(tn):offsets(tn))-b.ts(onsets(tn));
%                 figure(100);clf;plot(ts,t_speed);
                st = onsets(tn) - round(1e6 * 1/b.si); se = onsets(tn) + round(1e6 * 5/b.si);
                t_speed1 = b.fSpeed(st:100:se);
                tfspeed = [tfspeed t_speed1];
                fspeeds{an,cn,tn} = t_speed;
                fspeeds1{an,cn,tn} = t_speed1;
                t_speed = b.fSpeed(onsetsI(tn):offsetsI(tn));
                fspeedsI{an,cn,tn} = t_speed;
                n = 0;
            end
            fspeeds_cn{an,cn} = tfspeed; %b.fSpeed(onsets(1):offsets(end));
            n = 0;
        end
        fspeeds_an{an} = b.fSpeed;
        belt_lengths{an} = temp_bl;
    end
    out.diff_ms = ms_a - ms_b;
    out.mean_speed_middle = mean_speed_middle;
    out.mean_speed_end = mean_speed_end;
    out.corr_s_tp = c_s_tp;
    out.mspeed = mspeed;
    out.mspeedI = mspeedI;
    out.son = son;     out.soff = soff;
    out.movL = movL; out.durT = durT; out.durIT = durIT;
    out.below7 = below7; out.above7 = above7;
    out.mspeed_raster = mspeed_raster;
    out.fspeeds = fspeeds;out.fspeedsI = fspeedsI;
    out.fspeeds1 = fspeeds1;
    out.fspeeds_cn = fspeeds_cn;
    out.fspeeds_an = fspeeds_an;
    n = 0;
    
function out = get_higuchi(ei_C)
    for an = 1:length(ei_C)
        temp_bl = [];
        for cn = 1:4
            [an cn]
            onsets = ei_C{an}.plane{1}.contexts(cn).markers.air_onsets;
            offsets = ei_C{an}.plane{1}.contexts(cn).markers.air_offsets;
            onsetsI = ei_C{an}.plane{1}.contexts(cn).markers.airI_onsets;
            offsetsI = ei_C{an}.plane{1}.contexts(cn).markers.airI_offsets;
            b = ei_C{an}.b;
            temp_bl(cn,:) = b.dist(offsets)-b.dist(onsets);
            m_belt_length(an,cn) = mean(b.dist(offsets)-b.dist(onsets));
            sz(an,cn) = size(ei_C{an}.plane{1}.contexts(cn).rasters.airD.sp_rasters_nan_corrected,2);
            speeds_C{an,cn} = ei_C{an}.plane{1}.contexts(cn).rasters.airD.speed(:,1:sz(an,cn));
            this_speed_raster = speeds_C{an,cn};
            mspeed_raster(:,cn,an) = mean(this_speed_raster(:,1:35));
            c_s = corr(this_speed_raster');
            mask = ones(size(c_s)); mask = triu(mask,1) & ~triu(mask,2);
            c_s_tp(:,cn,an) = c_s(mask);
            middle_bin = floor(sz(an,cn)/2);
            start_bin = middle_bin - 2;
            end_bin = middle_bin + 2;
            middle_speeds = this_speed_raster(:,start_bin:end_bin);
            mean_speed_middle(an,cn) = nanmean(middle_speeds(:));
            end_speeds = this_speed_raster(:,(sz(an,cn)-4):end);
            mean_speed_end(an,cn) = nanmean(end_speeds(:));
            seconds_before = 3; seconds_after = 5; 
            parfor tn = 1:10
                hifd(tn,cn,an) = Higuchi_FD(b.fSpeed(onsets(tn):offsets(tn)),20);
            end
            
        end
        belt_lengths{an} = temp_bl;
    end
    out.hifd = hifd;
    n = 0;