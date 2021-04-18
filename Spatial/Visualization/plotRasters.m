function plotRasters(Ai,ccs)
% ccs = 1:size(A.rasters,3);
if ~iscell(Ai)
    A = Ai;
    numberOfRows = 4;
    numberOfCols = 4;
    graphsOnOneFigure = numberOfRows * numberOfCols;
    numberOfData = length(ccs);
    numberOfGroups = ceil(numberOfData/graphsOnOneFigure);
    numberOfFullGroups = floor(numberOfData/graphsOnOneFigure);
    indices = NaN(1,(numberOfGroups*graphsOnOneFigure));
    indices(1:numberOfData) = 1:numberOfData;
    groupIndices = reshape(indices,numberOfRows,numberOfCols,numberOfGroups);
    for gg = 1:numberOfGroups
        groupIndices(:,:,gg) = groupIndices(:,:,gg)';
    end

    ff = makeFigureRowsCols(1005,[NaN 0.1 7 5],'RowsCols',[numberOfRows numberOfCols],...
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
                cn = ccs(cni);
                thisRaster = A.rasters(:,:,cn);
                mSig = nanmean(thisRaster);
                axes(ff.h_axes(rr,cc));
                ft = fittype(A.formula);
                xs = 1:size(A.rasters,2);
                coeff = A.coeff(:,cn);
                fitplot = ft(coeff(1),coeff(2),coeff(3),xs);
                imagesc(thisRaster,[0 max(thisRaster(:))]);hold on;
                plot(size(thisRaster,1)*fitplot/max(fitplot),'linewidth',1.5,'color','r');
                plot(size(thisRaster,1)*mSig/max(mSig),'linewidth',1.5,'color','m');
                title(sprintf('%d-SI(%.2f)-Center(%.1f)-MaxFR(%.1f)',cn,A.SI(cn),A.centers(cn),max(thisRaster(:))));
                box off;
                set(gca,'FontSize',10,'FontWeight','Normal','linewidth',0.75,'Ydir','normal','TickDir','out');
                xticks = [1:10:size(thisRaster,2)];
                set(gca,'XTick',xticks,'XTickLabels',A.xs(xticks));
            end
        end
    %     showCells(1000,ei,ccs(cng));
    %     save2pdf('temp.pdf',gcf,600);
        display('Press key');
        gg = keyboardInput(gg,[1 numberOfGroups],[1 6],'');
        if gg < 0
            break;
        end
    end
else
    numberOfRows = 4;
    numberOfCols = size(Ai,2);
    graphsOnOneFigure = numberOfRows;
    numberOfData = length(ccs);
    numberOfGroups = ceil(numberOfData/graphsOnOneFigure);
    indices = NaN(1,(numberOfGroups*graphsOnOneFigure));
    indices(1:numberOfData) = 1:numberOfData;
    groupIndices = reshape(indices,numberOfRows,1,numberOfGroups);
%     for gg = 1:numberOfGroups
%         groupIndices(:,:,gg) = groupIndicesI(:,:,gg)';
%     end

    ff = makeFigureRowsCols(1005,[NaN 0.1 7 5],'RowsCols',[numberOfRows numberOfCols],...
        'spaceRowsCols',[0.1 0.05],'rightUpShifts',[0.05 0.07],'widthHeightAdjustment',...
        [-70 -115]);

    gg = 1;
    while 1
        for rr = 1:numberOfRows
            for cc = 1:numberOfCols
                cni = groupIndices(rr,1,gg);
                if isnan(cni)
                    cla(ff.h_axes(rr,cc));
                    set(ff.h_axes(rr,cc),'visible','off');
                    continue;
                end
                set(ff.h_axes(rr,cc),'visible','on');
                cn = ccs(cni);
                A = Ai{cc};
                thisRaster = A.rasters(:,:,cn);
                mSig = applyGaussFilt(nanmean(thisRaster),5);
                axes(ff.h_axes(rr,cc));
%                 ft = fittype(A.formula);
                xs = 1:size(A.rasters,2);
                uxs = linspace(1,length(xs),length(xs)*10);
%                 coeff = A.coeff(:,cn);
                imagesc(thisRaster,[0 max(thisRaster(:))]);colorbar;hold on;
                try
                    fitplot = gauss_fit(xs,A.gauss_fit_on_mean.coefficients_Rs_mean(cn,1:3),A.gauss_fit_on_mean.gauss1Formula);
                    plot(xs,size(thisRaster,1)*fitplot/max(fitplot),'linewidth',1.5,'color','r');
                catch
                end
                plot(size(thisRaster,1)*mSig/max(mSig),'linewidth',1.5,'color','m');
                try
                title(sprintf('%d-SI(%.2f)-Center(%.1f)-MaxFR(%.1f)-Rsq(%.2f)-HaFD(%.3f)-HiFD(%.3f)',cn,A.SI(cn),A.centers(cn),...
                    max(mSig),A.rs(cn),A.fractal_dim.HaFD(cn),A.fractal_dim.HiFD(cn)));
                catch
                    title(sprintf('%d',cn));
                end
                box off;
                set(gca,'FontSize',8,'FontWeight','Normal','linewidth',0.75,'Ydir','normal','TickDir','out');
                xticks = [1:10:size(thisRaster,2)];
                set(gca,'XTick',xticks,'XTickLabels',A.xs(xticks));
            end
        end
    %     showCells(1000,ei,ccs(cng));
    %     save2pdf('temp.pdf',gcf,600);
        display('Press key');
        gg = keyboardInput(gg,[1 numberOfGroups],[1 6],'');
        if gg < 0
            break;
        end
    end
end