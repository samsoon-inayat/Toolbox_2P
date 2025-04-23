%%
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1:3),1,2);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[10 4 1.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.3 0.35],'widthHeightAdjustment',[-550 -380]);
MY = 55; ysp = 0.5; mY = -1; ystf = 0.52; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),ra,{'AP','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
% make_bars_hollow(hbs(2))
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'AOn','AOff'});xtickangle(20);
ylabel({'Ske. Speed'});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0 0]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);
%%
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1:3),1,2);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[10 4 1.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.3 0.35],'widthHeightAdjustment',[-550 -380]);
MY = 1.5; ysp = 1; mY = -0.5; ystf = 0.41; ysigf = 0.01;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raR{1},{'PT','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
% make_bars_hollow(hbs(2))
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Time-Speed','Dist-Speed'});xtickangle(20);
ylabel({'Ske. Speed'});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0 0]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);

%%
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1:3),1,2);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[10 4 1.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.3 0.35],'widthHeightAdjustment',[-550 -380]);
MY = 1.5; ysp = 1; mY = -0.5; ystf = 0.41; ysigf = 0.01;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raR{2},{'DT','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
% make_bars_hollow(hbs(2))
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Time','Dist'});xtickangle(20);
ylabel({'Ske. Speed'});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0 0]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);


%%
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1:3),1,2);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[10 4 1.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.3 0.35],'widthHeightAdjustment',[-550 -380]);
MY = 0.1; ysp = 0.01; mY = -0.5; ystf = 0.01; ysigf = 0.01;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raR{2},{'PT','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
% make_bars_hollow(hbs(2))
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'FR-Time','FR-Dist','FR-Speed'});xtickangle(20);
ylabel({'Ske. Speed'});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Pooled'},{[0 0]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);


%%
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1:3),1,2);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[10 4 1.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.3 0.35],'widthHeightAdjustment',[-550 -380]);
MY = 1; ysp = 0.1; mY = 0; ystf = 0.01; ysigf = 0.01;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raR{1},{'PT','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
% make_bars_hollow(hbs(2))
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'FR-Time','FR-Dist','FR-Speed'});xtickangle(20);
ylabel({'Ske. Speed'});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Pooled'},{[0 0]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);


%%
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1:10),1,2);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[10 4 3.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.07 0.35],'widthHeightAdjustment',[-100 -380]);
MY = 1.7; ysp = 3.25; mY = 0; ystf = 3.12; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
% MY = 0.07; ysp = 3.25; mY = -0.07; ystf = 3.12; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raRR{1},{'TN','hsd',(0.05/2)},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'T01','T02','T03','T04','T05','T06','T07','T08','T09','T10'});xtickangle(20);
ylabel('Time (s)')
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,10,{'Time-Bin','Dist-Bin'},{[0 0]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);

%%
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1:10),1,3);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[4 4 6.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.3 0.35],'widthHeightAdjustment',[-550 -380]);
MY = 0.51; ysp = 0.025; mY = 0; ystf = 0.12; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raR{1},{'AP:TN','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
% make_bars_hollow(hbs(2))
format_axes(gca);
xlbl = {'T01','T02','T03','T04','T05','T06','T07','T08','T09','T10'};
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',xlbl);xtickangle(20);
% ylabel({'Mutual','Information'});
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0 0]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);


%%
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1:10),1,2);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[4 4 6.5 3],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.07 0.15],'widthHeightAdjustment',[-100 -380]);
MY = 0.7; ysp = 0.05; mY = 0; ystf = 0.02; ysigf = 0.015;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),ra,{'AP:TT','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',tuningTypes([2 3 5]));xtickangle(20);
ylabel('MI')
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,10,{'Time-Bin','Dist-Bin'},{[0 0]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);

%% Figure 2B Air-On PC
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.colors(1:2),1,2);
% figure(300);clf; ha = gca
ff = makeFigureRowsCols(2020,[5 4 1.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.2 0.35],'widthHeightAdjustment',[-250 -550]);
MY = 0.71; ysp = 0.155; mY = 0; ystf = 0.073; ysigf = 0.025;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raR{1},{'BT:PT','hsd',0.05},[1 1.5],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Ti-Sp','Di-Sp'});xtickangle(20);
shift_ticklabels(gca,-3,0)
ylabel('PC')
axes_title_shifts_line = [0 0 0 0]; axes_title_shifts_text = [0.1 0.6155 0 0];
ht = axes_title(ff,{1},{'Air-On (pooled)'},axes_title_shifts_line,axes_title_shifts_text,'no');
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Time-Bin','Dist-Bin'},{[0 0.00531]});
% set_sub_graph_text(ff,1,{'Pooled'},[0.05 -0.375 0.3 0],[0.25 -0.041 0 0]);
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Pooled'},{[0 -0.075]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('Pooled'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);
%% Figure 2D Air-Off PC
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.colors(1:2),1,2);
ff = makeFigureRowsCols(2020,[5 4 1.5 1],'RowsCols',[1 2],'spaceRowsCols',[0.07 0.071],'rightUpShifts',[0.2 0.35],'widthHeightAdjustment',[-170 -550]);
axes(ff.h_axes(1,1))
MY = 0.1; ysp = 0.05; mY = -0.55; ystf = 0.23; ysigf = 0.025;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raR{2},{'PT','hsd',0.05},[1 1.5],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Ti-Sp','Di-Sp'});xtickangle(30);
ylabel('PC')
shift_ticklabels(gca,-3,0);
axes_title_shifts_line = [-0.1 0 0.5 0]; axes_title_shifts_text = [0.1 0.6155 0 0];
ht = axes_title(ff,{1},{'Air-Off (pooled)'},axes_title_shifts_line,axes_title_shifts_text,'no');

tcolors = repmat(mData.colors(3:4),1,2);
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,2),raR{2},{'BT','hsd',0.05},[1 1.5],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Time-Bin','Dist-Bin'},'YTickLabel',[]);xtickangle(30);
% ylabel('PC')

% axes_title_shifts_line = [0 0 0 0]; axes_title_shifts_text = [0.1 0.45 0 0];ht = axes_title(ff,{1:2},{'Air-Off (pooled)'},axes_title_shifts_line,axes_title_shifts_text,'yes');
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);


%% Figure 2D Air-ON MI
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.colors(1:2),1,2);
ff = makeFigureRowsCols(2020,[5 4 1.5 1],'RowsCols',[1 2],'spaceRowsCols',[0.07 0.071],'rightUpShifts',[0.2 0.35],'widthHeightAdjustment',[-170 -550]);
axes(ff.h_axes(1,1))
MY = 1; ysp = 0.01; mY = 0; ystf = 0.13; ysigf = 0.025;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raR{1},{'PT','hsd',0.05},[1 1.5],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Ti-Sp','Di-Sp'});xtickangle(30);
ylabel('MI')
shift_ticklabels(gca,-3,0);
axes_title_shifts_line = [-0.1 0 0.5 0]; axes_title_shifts_text = [0.1 0.6155 0 0];
ht = axes_title(ff,{1},{'Air-On (pooled)'},axes_title_shifts_line,axes_title_shifts_text,'no');


tcolors = repmat(mData.colors(3:4),1,2);
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,2),raR{1},{'BT','hsd',0.05},[1 1.5],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Time-Bin','Dist-Bin'},'YTickLabel',[]);xtickangle(30);
% ylabel('PC')

% axes_title_shifts_line = [0 0 0 0]; axes_title_shifts_text = [0.1 0.45 0 0];ht = axes_title(ff,{1:2},{'Air-Off (pooled)'},axes_title_shifts_line,axes_title_shifts_text,'yes');
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);

winopen(fullfile(mData.pdf_folder,'bar_graphs.pdf'))

%% Figure 2D Air-OFF MI
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.colors(1:2),1,2);
ff = makeFigureRowsCols(2020,[5 4 1.5 1],'RowsCols',[1 2],'spaceRowsCols',[0.07 0.071],'rightUpShifts',[0.2 0.35],'widthHeightAdjustment',[-170 -550]);
axes(ff.h_axes(1,1))
MY = 1; ysp = 0.01; mY = 0; ystf = 0.13; ysigf = 0.025;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raR{2},{'PT','hsd',0.05},[1 1.5],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Ti-Sp','Di-Sp'});xtickangle(30);
ylabel('MI')
shift_ticklabels(gca,-3,0);
axes_title_shifts_line = [-0.1 0 0.5 0]; axes_title_shifts_text = [0.1 0.6155 0 0];
ht = axes_title(ff,{1},{'Air-Off (pooled)'},axes_title_shifts_line,axes_title_shifts_text,'no');


tcolors = repmat(mData.colors(3:4),1,2);
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,2),raR{2},{'BT','hsd',0.05},[1 1.5],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Time-Bin','Dist-Bin'},'YTickLabel',[]);xtickangle(30);
% ylabel('PC')

% axes_title_shifts_line = [0 0 0 0]; axes_title_shifts_text = [0.1 0.45 0 0];ht = axes_title(ff,{1:2},{'Air-Off (pooled)'},axes_title_shifts_line,axes_title_shifts_text,'yes');
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);

winopen(fullfile(mData.pdf_folder,'bar_graphs.pdf'))

%% Figure 3 Active Cells
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.colors(1:3),1,2);
ff = makeFigureRowsCols(2020,[5 4 1.5 1],'RowsCols',[1 2],'spaceRowsCols',[0.07 0.071],'rightUpShifts',[0.2 0.35],'widthHeightAdjustment',[-170 -550]);
axes(ff.h_axes(1,1))
MY = 115.25; ysp = 15; mY = 0; ystf = 15; ysigf = 0.025;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),ra,{'CN','hsd',0.05},[1 1.5],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'C3','C4','C5'});xtickangle(30);
ylabel('Cells (%)');
shift_ticklabels(gca,-3,0);
axes_title_shifts_line = [-0.1 0 0.5 0]; axes_title_shifts_text = [0.1 0.6155 0 0];
ht = axes_title(ff,{1},{'Active Cells (Pooled)'},axes_title_shifts_line,axes_title_shifts_text,'no');

tcolors = repmat(mData.colors(5:6),1,2);
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,2),ra,{'AP','hsd',0.05},[1 1.5],tcolors,[mY MY ysp ystf ysigf],mData);
set(ff.h_axes(1,2),'Position',(get(ff.h_axes(1,2),'Position')+[0 0 -0.05 0]))
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Air-On','Air-Off'},'YTickLabel',[]);xtickangle(30);
% ylabel('PC')

% axes_title_shifts_line = [0 0 0 0]; axes_title_shifts_text = [0.1 0.45 0 0];ht = axes_title(ff,{1:2},{'Air-Off (pooled)'},axes_title_shifts_line,axes_title_shifts_text,'yes');
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);

% winopen(fullfile(mData.pdf_folder,'bar_graphs.pdf'))
%% Figure 3 Active Cells distributions comparison
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.colors(1:11),1,3);
ff = makeFigureRowsCols(2020,[5 4 3.5 1],'RowsCols',[1 2],'spaceRowsCols',[0.07 0.071],'rightUpShifts',[0.1 0.35],'widthHeightAdjustment',[-120 -550]);
axes(ff.h_axes(1,1))
MY = 40; ysp = 1; mY = 0; ystf = 1; ysigf = 0.025;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raR{1},{'CN:NT','hsd',0.05},[1 1.5],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'C3','C4','C5'});xtickangle(30);
ylabel('Cells (%)');
shift_ticklabels(gca,-3,0);
axes_title_shifts_line = [-0.1 0 0.5 0]; axes_title_shifts_text = [0.1 0.6155 0 0];
ht = axes_title(ff,{1},{'Active Cells (Pooled)'},axes_title_shifts_line,axes_title_shifts_text,'no');

tcolors = repmat(mData.colors(1:11),1,2);
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,2),raR{2},{'NT','hsd',0.05},[1 1.5],tcolors,[mY MY ysp ystf ysigf],mData);
set(ff.h_axes(1,2),'Position',(get(ff.h_axes(1,2),'Position')+[0 0 -0.05 0]))
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Air-On','Air-Off'},'YTickLabel',[]);xtickangle(30);
% ylabel('PC')

% axes_title_shifts_line = [0 0 0 0]; axes_title_shifts_text = [0.1 0.45 0 0];ht = axes_title(ff,{1:2},{'Air-Off (pooled)'},axes_title_shifts_line,axes_title_shifts_text,'yes');
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);

% winopen(fullfile(mData.pdf_folder,'bar_graphs.pdf'))

%% Figure 3 conjunctive complementary1 complementary2 between air-on and air-off
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.colors(8:10),1,2);
ff = makeFigureRowsCols(2020,[5 4 1.5 1],'RowsCols',[1 2],'spaceRowsCols',[0.07 0.071],'rightUpShifts',[0.2 0.35],'widthHeightAdjustment',[-170 -550]);
axes(ff.h_axes(1,1))
MY = 115.25; ysp = 15; mY = 0; ystf = 5; ysigf = 0.025;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),ra,{'CT','hsd',0.05},[1 1.5],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Conj','Comp1','Comp2'});xtickangle(30);
ylabel('Cells (%)');
shift_ticklabels(gca,-2,0);
axes_title_shifts_line = [-0.1 0 0.5 0]; axes_title_shifts_text = [0.1 0.6155 0 0];
ht = axes_title(ff,{1},{'Population Overlap (Pooled)'},axes_title_shifts_line,axes_title_shifts_text,'no');

tcolors = repmat(mData.colors(1:3),1,2);
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,2),ra,{'CN','hsd',0.05},[1 1.5],tcolors,[mY MY ysp ystf ysigf],mData);
set(ff.h_axes(1,2),'Position',(get(ff.h_axes(1,2),'Position')+[0 0 -0.05 0]))
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'C3','C4','C5'},'YTickLabel',[]);xtickangle(30);
% ylabel('PC')

% axes_title_shifts_line = [0 0 0 0]; axes_title_shifts_text = [0.1 0.45 0 0];ht = axes_title(ff,{1:2},{'Air-Off (pooled)'},axes_title_shifts_line,axes_title_shifts_text,'yes');
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);

% winopen(fullfile(mData.pdf_folder,'bar_graphs.pdf'))

%% Figure 4 Tuning type time, distance, speed cells percentages
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.colors(1:3),1,3);
ff = makeFigureRowsCols(2020,[5 4 6.5 5],'RowsCols',[1 1],'spaceRowsCols',[0.07 0.091],'rightUpShifts',[0.0512 0.35],'widthHeightAdjustment',[-100 -250]);
MY = 0.25; ysp = 0.015; mY = -0.1; ystf = 0.015; ysigf = 0.025;titletxt = ''; ylabeltxt = {'PC'}; % for all cells (vals) MY = 80
MY = 25; ysp = 0.51; mY = 0; ystf = 0.05; ysigf = 0.025;titletxt = ''; ylabeltxt = {'MI'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),ra_AP{1},{'BT:TT','hsd',0.05},[1 1.75],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'T','D','S'});xtickangle(30);
% ylabel('Cells (%)');
shift_ticklabels(gca,-3,0);
axes_title_shifts_line = [0.1 0 0.5 0]; axes_title_shifts_text = [0.1 0.76155 0 0];
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'PC','MI'},{[0.1 0.01]});
ht = axes_title(ff,{1},{'Air-On (Pooled)'},axes_title_shifts_line,axes_title_shifts_text,'no');

% [hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,2),ra_AP1_BT{1},{'CT','hsd',0.05},[1 1.5],tcolors,[mY MY ysp ystf ysigf],mData);
% format_axes(gca);
% set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
%     'XTick',xdata,'XTickLabel',{'T','D','S','TD','TS','DS','TDS'});xtickangle(30);
% % ylabel('Cells (%)');
% shift_ticklabels(gca,-3,0);
% axes_title_shifts_line = [0.1 0 0.25 0]; axes_title_shifts_text = [0.1 0.26155 0 0];
% ht = axes_title(ff,{2},{'Air-Off (Pooled)'},axes_title_shifts_line,axes_title_shifts_text,'no');

%% Figure 4 Tuning type time, distance, speed cells percentages
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(8:10),1,3);
ff = makeFigureRowsCols(2020,[5 4 1.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0.091],'rightUpShifts',[0.175 0.35],'widthHeightAdjustment',[-450 -550]);
MY = 0.15; ysp = 0.015; mY = -0.051; ystf = 0.015; ysigf = 0.0025;titletxt = ''; ylabeltxt = {'PC'}; % for all cells (vals) MY = 80
MY = 20; ysp = 5; mY = 0; ystf = 5; ysigf = 0.025;titletxt = ''; ylabeltxt = {'Cells (%)'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),ra_AP{1},{'TT','hsd',0.05},[1 1.5],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Time','Dist','Speed'});xtickangle(30);
ylabel(ylabeltxt);
shift_ticklabels(gca,-3,0);
axes_title_shifts_line = [0 0 0.5 0]; axes_title_shifts_text = [0 0.6155 0 0];
ht = axes_title(ff,{1},{'Air-On (Pooled)'},axes_title_shifts_line,axes_title_shifts_text,'no');

save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);

% winopen(fullfile(mData.pdf_folder,'bar_graphs.pdf'))


%% Response Fidelity
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.colors(1:3),1,2);
ff = makeFigureRowsCols(2020,[5 4 1.5 1],'RowsCols',[1 2],'spaceRowsCols',[0.07 0.071],'rightUpShifts',[0.2 0.35],'widthHeightAdjustment',[-170 -550]);
axes(ff.h_axes(1,1))
MY = 80; ysp = 10.05; mY = 0; ystf = 10.23; ysigf = 0.025;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raR{2},{'CN','hsd',0.05},[1 1.5],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'C3','C4','C5'});xtickangle(30);
ylabel('Trials (%)')
shift_ticklabels(gca,-3,0);
axes_title_shifts_line = [-0.1 0 0.5 0]; axes_title_shifts_text = [0.1 0.6155 0 0];
ht = axes_title(ff,{1},{'Air-Off (pooled)'},axes_title_shifts_line,axes_title_shifts_text,'no');
delete(ff.h_axes(1,2))
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);

%% Low, Medium, High Response Fidelity
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.colors(1:3),1,3);
ff = makeFigureRowsCols(2020,[5 4 1.5 1],'RowsCols',[1 2],'spaceRowsCols',[0.07 0.071],'rightUpShifts',[0.2 0.35],'widthHeightAdjustment',[-170 -550]);
axes(ff.h_axes(1,1))
MY = 80; ysp = 10.05; mY = 0; ystf = 10.23; ysigf = 0.025;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raR{2},{'RF','hsd',0.05},[1 1.5],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'L','M','H'});xtickangle(30);
ylabel('Trials (%)')
shift_ticklabels(gca,-3,0);
axes_title_shifts_line = [-0.1 0 0.5 0]; axes_title_shifts_text = [0.1 0.6155 0 0];
ht = axes_title(ff,{1},{'Air-Off (pooled)'},axes_title_shifts_line,axes_title_shifts_text,'no');
delete(ff.h_axes(1,2))
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);
%%
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1:3),1,10);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[2 4 6.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.07 0.35],'widthHeightAdjustment',[-100 -380]);
MY = 1; ysp = 0.25; mY = -1; ystf = 0.12; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raR{2},{'BT:TT','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Ti-Di','Ti-Sp','Di-Sp'});xtickangle(20);
ylabel('MI')
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,10,{'Time-Bin','Dist-Bin'},{[0 0]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);


%%
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1:10),1,3);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[10 4 6.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.07 0.35],'widthHeightAdjustment',[-100 -380]);
MY = 1.5; ysp = 0.25; mY = -1; ystf = 0.12; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raR{1},{'PT:TN','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Ti-Di','Ti-Sp','Di-Sp'});xtickangle(20);
ylabel('MI')
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,10,{'Time-Bin','Dist-Bin'},{[0 0]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);

%%
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1:10),1,2);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[10 4 3.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.07 0.35],'widthHeightAdjustment',[-100 -380]);
MY = 0.5; ysp = 0.5; mY = 0; ystf = 0.012; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raR{1},{'CN:PT','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Time-Dist','Time-Speed','Dist-Speed'});xtickangle(20);
ylabel('MI')
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,10,{'Time-Bin','Dist-Bin'},{[0 0]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);

%%  raw traces Figure 2A
% plotting the concatenated trials
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1:3),1,2);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[1 4 6.9 2],'RowsCols',[2 4],'spaceRowsCols',[0.12 0.021],'rightUpShifts',[0.05 0.15],'widthHeightAdjustment',[-40 -180]);
an = 1; cn = 1; aps = [1 1 2 2]; pos_air = [0.051 0.0155 0.3 0]; pos_PCMI = [0.01 -0.035 0.1 0]; pos_an_heading = [-0.03 0.075 0.3 0];
for gn = 1:4
    ap = aps(gn);
    if mod(gn,2) == 0
        out = outD;
    else
        out = outT;
    end
    timecc = out.atimecc{an,cn,ap}; distcc = out.adistcc{an,cn,ap}; speedcc = out.aspeedcc{an,cn,ap}; speedcc(speedcc < 0) = 0;
    trialsc = out.atrialcc{an,cn,ap}; tr_idx = find(trialsc == 5,1,'last');
    xlims = [1 tr_idx];
    axes(ff.h_axes(1,gn))
    xs = 1:length(timecc);
    [ax, h1, h2] = plotyy(xs,timecc,xs,speedcc); box off;
    set(ax(1),'xlim',xlims);
    if gn == 1
        ylabel(ax(1),'Time (s)');
        set(ax(2),'YTickLabel',[]);
    else
        if gn < 4
            set(ax(2),'YTickLabel',[]);
            set(ax(1),'YTickLabel',[]);
        else
            set(ax(1),'YTickLabel',[]);
            set(ax(2),'YTick',[0 20 40 60]);
        end
    end
    format_axes(ax(1));
    set(ax(1),'YColor','b','ylim',[0 20]); set(h1, 'Color', 'b');
    set(ax(2),'xlim',xlims,'ylim',[0 60]); set(h2, 'Color', 'r');
    if gn == 4
        ylabel(ax(2),'Speed (cm/s)');
    end
    format_axes(ax(2)); set(ax(2),'YColor','r');
    PC = corr(timecc,speedcc); MI = calc_metric_MI(timecc,speedcc,10,0);
    % set_axes_top_text_no_line(ff.hf,ax(1),'Animal - 01',[0.031 0.051 0.3 0])
    if gn == 1
        set_axes_top_text_no_line(ff.hf,ax(1),'Animal # 1 - C3 - First 5 Concatenated Trials',pos_an_heading)
        set_axes_top_text_no_line(ff.hf,ax(1),'Air-On',pos_air)
    else
        if gn >= 3
            set_axes_top_text_no_line(ff.hf,ax(1),'Air-Off',pos_air)
        else
            set_axes_top_text_no_line(ff.hf,ax(1),'Air-On',pos_air)
        end
    end

    ht = set_axes_top_text_no_line(ff.hf,ax(1),sprintf('PC = %.4f, MI = %.4f',PC,MI),pos_PCMI); set(ht,'color','m')
    
    axes(ff.h_axes(2,gn))
    [ax, h1, h2] = plotyy(xs,distcc,xs,speedcc);  box off;
    set(ax(1),'xlim',xlims); 
    if gn == 1
        ylabel(ax(1),'Distance (cm)'); 
        set(ax(2),'YTickLabel',[]);
    else
        if gn < 4
            set(ax(2),'YTickLabel',[]);
            set(ax(1),'YTickLabel',[]);
        else
            set(ax(1),'YTickLabel',[]);
            set(ax(2),'YTick',[0 20 40 60]);
        end
    end
        format_axes(ax(1));
    set(ax(1),'YColor','k','ylim',[0 200]); set(h1, 'Color', 'k');
    set(ax(2),'xlim',xlims,'ylim',[0 60]); set(h2, 'Color', 'r');
    if gn == 4
        ylabel(ax(2),'Speed (cm/s)');
    end
    format_axes(ax(2)); set(ax(2),'YColor','r');
    if mod(gn,2) == 0
        xlabel('Distance Bins');
    else
        xlabel('Time Bins');
    end
    PC = corr(distcc,speedcc); MI = calc_metric_MI(distcc,speedcc,10,0);
    % set_axes_top_text_no_line(ff.hf,ax(1),'Animal-01, C3, Air-On',[0.1 0 0.1 0])
    ht = set_axes_top_text_no_line(ff.hf,ax(1),sprintf('PC = %.4f, MI = %.4f',PC,MI),pos_PCMI); set(ht,'color','m')
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('raw_plots.pdf'),600);

