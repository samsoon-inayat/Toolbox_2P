function [odata,mData] = getRasterData(iei,context,stimMarker)
mData = [];


allContexts = contextDefinitions;


for aa = 1:length(iei)
    ei = iei{aa};
    belt_length = ei.b.belt_length;
    datap = [];
    for pp = 1:length(ei.plane)
        
        if strcmp(stimMarker,'airOnsets27') || strcmp(stimMarker,'airOffsets27') || strcmp(stimMarker,'airOnsets11') || strcmp(stimMarker,'airOffsets11')...
                || strcmp(stimMarker,'light') || strcmp(stimMarker,'airOnsets010') || strcmp(stimMarker,'airOffsets010') ...
                || strcmp(stimMarker,'airOnsets01') || strcmp(stimMarker,'airOffsets01')
            [~,data] = getDataContexts(ei.plane{pp},context,stimMarker);
        else
            [data,~] = getDataContexts(ei.plane{pp},context,stimMarker);
        end
        if pp == 1 & (strcmp(stimMarker,'air') || strcmp(stimMarker,'belt'))
            data = getSpeedRastersContexts(data,ei.plane{pp},context);
        end
        for ii = 1:length(data)
            [aa pp ii]
            if isequal([aa pp ii],[4 1 1])
                n = 0;
            end
            mrfs = data{ii}.gauss_fit_on_mean;
            [rs,coeff] = getMRFS_vals(mrfs);
            as = coeff(:,1);
            bs = coeff(:,2);
            cs = coeff(:,3);
            PWs = 2.36*cs./sqrt(2)*(belt_length/50);
            data{ii}.rs = rs;
            data{ii}.coeff = coeff';
            data{ii}.pws = PWs';
            data{ii}.centers = (bs * belt_length/50)';
            data{ii}.peaks = as';
            data{ii}.sel = (rs > 0.5) & (bs > 1 & bs < 50) & (PWs > 5 & PWs < 120)';% & (data{ii}.SI > 4));
            data{ii}.placeCells5 = logical((rs > 0.5) & (bs > 1 & bs < 50)' & (PWs > 5 & PWs < 120)' & (data{ii}.SI > 5));
            data{ii}.placeCells4 = logical((rs > 0.5) & (bs > 1 & bs < 50)' & (PWs > 5 & PWs < 120)' & (data{ii}.SI > 4));
            data{ii}.placeCells3 = logical((rs > 0.5) & (bs > 1 & bs < 50)' & (PWs > 5 & PWs < 120)' & (data{ii}.SI > 3));
            data{ii}.percentPlaceCells5 = sum(data{ii}.placeCells5)/length(data{ii}.placeCells5);
            data{ii}.percentPlaceCells4 = sum(data{ii}.placeCells4)/length(data{ii}.placeCells4);
            data{ii}.percentPlaceCells3 = sum(data{ii}.placeCells3)/length(data{ii}.placeCells3);
            data{ii}.formula = mrfs.gauss1Formula;
            data{ii}.belt_length = belt_length;
        end
        datap{pp} = data;
    end
    odata{aa} = datap;
end



function  data = getSpeedRastersContexts(data,ei,context)
for ii = 1:length(context)
    data{ii}.speedRasters = ei.contexts(context(ii)).rastersSpeed;
end
