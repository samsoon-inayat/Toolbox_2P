function see_imaging_quality

protocol_C = '10_C';
protocol_A = '10_A';
ei_C = evalin('base','ei10_C1');
ei_A = evalin('base','ei10_A');
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ET_C = evalin('base',sprintf('ET10_CD'));
ET_A = evalin('base',sprintf('ET10_CC'));
selAnimals_C = 1:length(ei_C)
selAnimals_A = 1:length(ei_A)

figure(100);clf;
plotNums = (reshape(1:10,2,5))';
for ii = 1:10
    subplot(2,5,ii)
    if ii <=5
        try
        cell_list = ei_C{ii}.plane{1}.tP.iscell(:,1); cell_list(isnan(cell_list)) = 0;
        catch
            continue;
        end
        showCells(gca,ei_C{ii},1,[]);
%         showCells(gca,ei_C{ii},1,1:length(cell_list));
        title(ii);
    else
        cell_list = ei_A{ii-5}.plane{1}.tP.iscell(:,1); cell_list(isnan(cell_list)) = 0;
        showCells(gca,ei_A{ii-5},1,[]);
%         showCells(gca,ei_A{ii-5},1,1:length(cell_list));
        title(ii);
    end
end