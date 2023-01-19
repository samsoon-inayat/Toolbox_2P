%% Three-Way RM-ANOVA Most Recent compare all variables pooled air on or off brake vs no-brake exc and inh (and controls) big ANOVA test
    ntrials = 50; 
    event_type = {'B-AOn-Exc','B-AOn-Inh','B-AOff-Exc','B-AOff-Inh','B-Arb-Exc','B-Arb-Inh','NB-AOn-Exc','NB-AOn-Inh','NB-AOff-Exc','NB-AOff-Inh','NB-Arb-Exc','NB-Arb-Inh'};
    sic = {[Lb];[Lbs];[Ab_On Ab_Off];[Abs_On Abs_Off]};
    clear all_gFR all_exc all_inh all_gV all_exc_inh
    prop_names = {'resp','N_Resp_Trials','zMI','zMINaN','HaFD','HiFD','cells_pooled'};
    cell_sel = {'good_FR','vals','exc','inh'};
    varName = {'all_gFR','all_gV','all_exc','all_inh'};
    for cii = 1:length(cell_sel);
        cmdTxt = sprintf('clear %s',varName{cii});eval(cmdTxt);
    end
    pni = 7;
    all_exc_inh = [];
    for ii = 1:length(sic)
        sit = sic{ii};
        tRs = o.Rs(:,sit);
        props1 = get_props_Rs(tRs,ntrials);
        for cii = 1:length(cell_sel)
            cmdTxt = sprintf('gFR = props1.%s;',cell_sel{cii});eval(cmdTxt);
            if pni == 1
                rf = find_percent(gFR);
                cmdTxt = sprintf('%s(:,ii) = mean(rf,2);',varName{cii}); eval(cmdTxt)
            else
                if pni == 7
                    cmdTxt = sprintf('%s(:,ii) = cell_list_op(gFR,[],''or'',1);',varName{cii});eval(cmdTxt);
                else
                    if pni == 3 || pni == 4
                        cmdtxt = sprintf('rf = props1.%s;',prop_names{3}); eval(cmdtxt);
                        if pni == 3
                            cmdTxt = sprintf('%s(:,ii) = mean(exec_fun_on_cell_mat(rf,''nanmean'',gFR),2);',varName{cii});eval(cmdTxt);
                        else
                            for rrr = 1:size(rf,1)
                                for ccc = 1:size(rf,2)
                                    temp = isnan(rf{rrr,ccc});
                                    rf1(rrr,ccc) = 100*sum(temp(gFR{rrr,ccc}))/size(rf{rrr,ccc},1);
                                end
                            end
                            cmdTxt = sprintf('%s(:,ii) = mean(rf1,2);',varName{cii}); eval(cmdTxt)
                        end
                    else
                        cmdtxt = sprintf('rf = props1.%s;',prop_names{pni}); eval(cmdtxt);
                        cmdTxt = sprintf('%s(:,ii) = mean(exec_fun_on_cell_mat(rf,''mean'',gFR),2);',varName{cii});eval(cmdTxt);
                    end
                end
            end
        end
        all_exc_inh = [all_exc_inh all_exc(:,ii) all_inh(:,ii)];
    end
%%
 %
avar = find_percent(all_gV);
[within,dvn,xlabels,withinD] = make_within_table({'ST','EL'},[2,2]); withinD3 = withinD;
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{'hsd'}});
ra.ranova
print_for_manuscript(ra)

%%
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[5 5 1.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.01 0.25],'widthHeightAdjustment',[10 -320]);
MY = 20; ysp = 3; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'Percent of Cells'};

stp = 0.3;magfac; widths = ([0.75 ])*magfac; gap = 0.01*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'ST','hsd'},[1.5 1 1]);
xdata = make_xdata([2],[1 1.5]);   
tcolors = mData.dcolors(4:end);
axes(ff.h_axes(1,1))
[hbs,maxY]  = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',4,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
% make_bars_hollow(hbs(13:end));
set_axes_limits(gca,[0.25 xdata(end)+.75],[mY MY]); format_axes(gca); xticks = xdata; 
set(gca,'xtick',xticks,'xticklabels',{'Light','Air'});xtickangle(45)
put_axes_labels(gca,{[],[0 0 0]},{ylabeltxt,[0 0 0]});
ht = set_axes_top_text_no_line(gcf,gca,'',[0 -0.051 0 0]);
save_pdf(ff.hf,mData.pdf_folder,sprintf('light_vs_air.pdf'),600);

%% Three-Way RM-ANOVA Most Recent compare all variables pooled air on or off brake vs no-brake exc and inh (and controls) big ANOVA test
    ntrials = 50; 
    event_type = {'B-AOn-Exc','B-AOn-Inh','B-AOff-Exc','B-AOff-Inh','B-Arb-Exc','B-Arb-Inh','NB-AOn-Exc','NB-AOn-Inh','NB-AOff-Exc','NB-AOff-Inh','NB-Arb-Exc','NB-Arb-Inh'};
    sic = {[Lb];[Lbs];[Ab_On];[Abs_On];[Ab_Off];[Abs_Off];[Ab_Offc];[Abs_Offc]};
    clear all_gFR all_exc all_inh all_gV all_exc_inh
    prop_names = {'resp','N_Resp_Trials','zMI','zMINaN','HaFD','HiFD','cells_pooled'};
    cell_sel = {'good_FR','vals','exc','inh'};
    varName = {'all_gFR','all_gV','all_exc','all_inh'};
    for cii = 1:length(cell_sel);
        cmdTxt = sprintf('clear %s',varName{cii});eval(cmdTxt);
    end
    pni = 7;
    all_exc_inh = [];
    for ii = 1:length(sic)
        sit = sic{ii};
        tRs = o.Rs(:,sit);
        props1 = get_props_Rs(tRs,ntrials);
        for cii = 1:length(cell_sel)
            cmdTxt = sprintf('gFR = props1.%s;',cell_sel{cii});eval(cmdTxt);
            if pni == 1
                rf = find_percent(gFR);
                cmdTxt = sprintf('%s(:,ii) = mean(rf,2);',varName{cii}); eval(cmdTxt)
            else
                if pni == 7
                    cmdTxt = sprintf('%s(:,ii) = cell_list_op(gFR,[],''or'',1);',varName{cii});eval(cmdTxt);
                else
                    if pni == 3 || pni == 4
                        cmdtxt = sprintf('rf = props1.%s;',prop_names{3}); eval(cmdtxt);
                        if pni == 3
                            cmdTxt = sprintf('%s(:,ii) = mean(exec_fun_on_cell_mat(rf,''nanmean'',gFR),2);',varName{cii});eval(cmdTxt);
                        else
                            for rrr = 1:size(rf,1)
                                for ccc = 1:size(rf,2)
                                    temp = isnan(rf{rrr,ccc});
                                    rf1(rrr,ccc) = 100*sum(temp(gFR{rrr,ccc}))/size(rf{rrr,ccc},1);
                                end
                            end
                            cmdTxt = sprintf('%s(:,ii) = mean(rf1,2);',varName{cii}); eval(cmdTxt)
                        end
                    else
                        cmdtxt = sprintf('rf = props1.%s;',prop_names{pni}); eval(cmdtxt);
                        cmdTxt = sprintf('%s(:,ii) = mean(exec_fun_on_cell_mat(rf,''mean'',gFR),2);',varName{cii});eval(cmdTxt);
                    end
                end
            end
        end
        all_exc_inh = [all_exc_inh all_exc(:,ii) all_inh(:,ii)];
    end
%%
jj= 1;
for ii = 1:2:8
    conjness(:,jj) = find_percent(cell_list_op(all_gV(:,ii:(ii+1)),[],'and',1));
    jj = jj + 1;
end
    
avar = conjness;
[within,dvn,xlabels,withinD] = make_within_table({'ET'},[4]); withinD3 = withinD;
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{'hsd'}});
ra.ranova
print_for_manuscript(ra)

%%
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[5 5 1.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.01 0.25],'widthHeightAdjustment',[10 -420]);
MY = 5; ysp = 1; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'Percent of Cells'};

stp = 0.3;magfac; widths = ([1.25 ])*magfac; gap = 0.01*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'ET','hsd'},[1.5 1 1]);
xdata = make_xdata([4],[1 1.5]);   
tcolors = mData.dcolors(4:end);
axes(ff.h_axes(1,1))
[hbs,maxY]  = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
% make_bars_hollow(hbs(13:end));
set_axes_limits(gca,[0.25 xdata(end)+.75],[mY MY]); format_axes(gca); xticks = xdata; 
set(gca,'xtick',xticks,'xticklabels',{'Light','AOn','AOff','Arb'});xtickangle(45)
put_axes_labels(gca,{[],[0 0 0]},{ylabeltxt,[0 0 0]});
ht = set_axes_top_text_no_line(gcf,gca,'Conjunction',[0 -0.051 0 0]);
save_pdf(ff.hf,mData.pdf_folder,sprintf('light_vs_air.pdf'),600);


