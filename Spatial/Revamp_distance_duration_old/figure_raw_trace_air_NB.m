function figure_raw_trace_air_NB

udata = evalin('base','udata1');
an = 1;
data_an = udata{an};

field_names = fieldnames(data_an);
for ii = 1:length(field_names)
    varname = field_names{ii};
    cmdTxt = sprintf('%s = data_an.%s;',varname,varname);eval(cmdTxt);
end

n = 0;
%%
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1:10),1,2);
% figure(300);clf; ha = gca;
ytxt = 49; fs = 6;
ff = makeFigureRowsCols(2020,[10 4 5.9 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.05 0.3],'widthHeightAdjustment',[-70 -350]);
plot(ts,speed,'color','b','LineWidth',0.5); hold on;
plot(ts,C3*43,'color',colors{1}); re = find_rising_edge(C3,0.5,0.05); text(ts(re),ytxt,'C3','FontSize',fs);
plot(ts,C4*43,'color',colors{2}); re = find_rising_edge(C4,0.5,0.05); text(ts(re),ytxt,'C4','FontSize',fs);
plot(ts,C5*43,'color',colors{3}); re = find_rising_edge(C5,0.5,0.05); text(ts(re),ytxt,'C5','FontSize',fs);
xlim([250,1180]); ylim([0 50]);
format_axes(gca);
box off; ylabel('Speed (cm/s)'); xlabel('Time (s)');

onsets = find_rising_edge(air,0.5,0.05);
offsets = find_falling_edge(air,-0.5,0.05);

ylims = ylim; ylims(2) = 43
[TLx TLy] = ds2nfu(ts(onsets(1)),ylims(2)-0);
[BLx BLy] = ds2nfu(ts(onsets(1)),ylims(1));
aH = (TLy - BLy);
% for ii = 1:length(onsets)
for ii = [11:40]%length(onsets)
    [BRx BRy] = ds2nfu(ts(offsets(ii)),ylims(1));
    [BLx BLy] = ds2nfu(ts(onsets(ii)),ylims(1));
    aW = (BRx-BLx);
    annotation('rectangle',[BLx BLy aW aH],'facealpha',0.2,'linestyle','none','facecolor','k');
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('fig_raw_trace_speeds.pdf'),600);




% MY = 235.7; ysp = 3.25; mY = 0; ystf = 3.12; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
% [hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),ra,{'AP:BT','hsd',(0.05/2)},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
% format_axes(gca);
% set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
%     'XTick',xdata,'XTickLabel',{'T01','T02','T03','T04','T05','T06','T07','T08','T09','T10'});xtickangle(20);
% ylabel('Time (s)')
% % set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,10,{'Time-Bin','Dist-Bin'},{[0 0]});
% % ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
% save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);