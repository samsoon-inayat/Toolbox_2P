function behaviorPlot(aei)
if ~exist('aei','var')
    aei = evalin('base','ei');
end

if ~iscell(aei)
    b = aei;
    clear aei;
    aei{1}.b = b;
end

for ee = 1:length(aei)
    ei = aei{ee};
    try
        plotMarkers(ei.b,ei.b.air_puff_r,ei.b.stim_r,1001);
    catch
        plotMarkers(ei.b,ei.b.air_puff_r,[],1001);
    end

    ylim([-0.1 1.4]);
    % xlim([300 550]);
    set(gca,'FontSize',12,'FontWeight','Bold','linewidth',1.5);
    ylabel('A.U');

    legendText = {'Air Puff','Belt Marker','Normalized Speed','Stim (e.g. LED or Sound)'};
    thisCols = {'b','r','m','g'};
    x1 = 350; x2 = x1+15; y1 = (1.1:0.05:100); y1 = y1(1:4); y2 = y1;
    legendFontSize = 11;
    legs = legendText;
    for ii = 1:length(legs)
        plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols{ii},'linewidth',2);
        text(x2+0.5,y1(ii),sprintf('%s',legs{ii}),'Color',thisCols{ii},'FontSize',legendFontSize);
    end
    box off;
    title(ei.recordingFolder);
    if length(aei) > 1
        pause;
    end
end
% save2pdf('behavior.pdf',gcf,600);