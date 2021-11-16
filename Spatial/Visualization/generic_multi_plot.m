function generic_multi_plot(figNum,NRC,fhandle,data)

numberOfRows = NRC(1);
numberOfCols = NRC(2);
graphsOnOneFigure = numberOfRows * numberOfCols;
numberOfData = NRC(3);
numberOfGroups = ceil(numberOfData/graphsOnOneFigure);
graph_nums = 1:numberOfData;
indices = (reshape(1:graphsOnOneFigure,numberOfCols,numberOfRows))';
groupIndices = reshape(NaN(graphsOnOneFigure*numberOfGroups,1),numberOfRows,numberOfCols,numberOfGroups);
ind = 1;
for gg = 1:numberOfGroups
    for rr = 1:numberOfRows
        for cc = 1:numberOfCols
            if ind <= length(graph_nums)
                groupIndices(rr,cc,gg) = graph_nums(ind);
                ind = ind + 1;
            end
        end
    end
end

ff = makeFigureRowsCols(figNum,[5 5 7 5],'RowsCols',[numberOfRows numberOfCols],...
    'spaceRowsCols',[0.1 0.05],'rightUpShifts',[0.05 0.07],'widthHeightAdjustment',...
    [-70 -115]);

gg = 1;
while 1
    for rr = 1:numberOfRows
        for cc = 1:numberOfCols
            cni = groupIndices(rr,cc,gg);
            if isnan(cni)
                cla(ff.h_axes(rr,cc));
                set(ff.h_axes(rr,cc),'visible','off');
                continue;
            end
            set(ff.h_axes(rr,cc),'visible','on');
            cmdTxt = sprintf('%s(ff.h_axes(rr,cc),data,cni);',fhandle);
            eval(cmdTxt);
            n = 0;
        end
    end
    display('Press key');
    gg = keyboardInput(gg,[1 numberOfGroups],[1 6],'');
    if gg < 0
        break;
    end
end
if gg < 0
    return;
end

function plotSpeedTuning(ah,d,ii)
axes(ah);
cla;
rs = [d.cl(ii) d.cs(ii) d.cg(ii)];
[~,mind] = max(rs);
plot(d.bcs,d.FR(ii,:),'c');hold on;
pl(1) = plot(d.bcs,d.fFRl(ii,:),'k','linewidth',1);
pl(2) = plot(d.bcs,d.fFRs(ii,:),'m','linewidth',1);
pl(3) = plot(d.bcs,d.fFRg(ii,:),'b','linewidth',1);
set(pl(mind),'linewidth',2);
title(sprintf('Cell %d',ii));
% legend('FR','Lin','Sig','Gauss');
% xlabel('Speed (cm/sec)');
% ylabel('FR (AU)');
