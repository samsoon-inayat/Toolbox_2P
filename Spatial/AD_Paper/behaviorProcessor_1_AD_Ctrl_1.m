function behaviorProcessor_1_AD_Ctrl
mData = evalin('base','mData');
training_data = evalin('base','training_data_C');
% {'173706';'183761';'183745';'183628';'183762'};
selRows = [4 9 8 6 10]; selCols = [1:3];
data_C = get_training_data(training_data,selRows);

training_data = evalin('base','training_data_A');
% {'183227';'183228';'183329';'001567';'001569'};
selRows = [4 5 6 1 2]; selCols = [1:3];
data_A = get_training_data(training_data,selRows);
n = 0;

%%
while 1
    HiFD_C = find_HiFD_here(data_C,10);
    HiFD_A = find_HiFD_here(data_A,10);
    %%
    ms_C = find_mean_speed(data_C,10);
    ms_A = find_mean_speed(data_A,10);
    %%
    break
end
%%
while 1
    var_C = []; var_A = [];
    for dd = 1:3
        var_C = [var_C squeeze(ms_C(:,dd,2:end))];
        var_A = [var_A squeeze(ms_A(:,dd,2:end))];
    end
    
   
    [within,dvn,xlabels,withinD] = make_within_table({'Days','Trials'},[3,9]); awithinD = withinD;
    dataT = make_between_table({var_C;var_A},dvn);
    ra = RMA(dataT,within,{0.05,{}});
    ra.ranova
    %%
    [within,dvn,xlabels,withinD] = make_within_table({'Trials'},[3]); awithinD = withinD;
    dataT = make_between_table({var_C;var_A},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd'}});
    ra.ranova
    %%
    break;
end



function speeds = get_training_data(temp,selRows,selCols)
numberOfTrials = findNumberOfTrials(temp);
[rr,cc] = find(numberOfTrials < 20);
ei1 = temp.bs;
i = 1;
for ii = 1:length(selRows)
    for jj = 1:3
        ei(i,jj) = ei1(selRows(ii),jj);
    end
    i = i + 1;
end
for rr = 1:size(ei,1)
    for cc = 1:size(ei,2)
        b = ei{rr,cc};
        if isempty(b)
            continue;
        end
        for ii = 1:length(b.air_puff_r)
            speeds{rr,cc,ii} = b.fSpeed(b.air_puff_r(ii):b.air_puff_f(ii));
        end
    end
end


function HiFD = find_HiFD_here(data,NTrials)
HiFD = NaN(size(data,1),size(data,2),NTrials);
for rr = 1:size(data,1)
    for cc = 1:size(data,2)
        parfor tt = 1:NTrials
            HiFD(rr,cc,tt) = Higuchi_FD(data{rr,cc,tt},20);
        end
    end
end

function HiFD = find_mean_speed(data,NTrials)
HiFD = NaN(size(data,1),size(data,2),NTrials);
for rr = 1:size(data,1)
    for cc = 1:size(data,2)
        for tt = 1:NTrials
            HiFD(rr,cc,tt) = nanmean(data{rr,cc,tt});
        end
    end
end




function numberOfTrials = findNumberOfTrials(training_data)
bs = training_data.bs;
numberOfTrials = NaN(size(bs));
for ii = 1:size(bs,1)
    for jj = 1:size(bs,2)
        thisb = bs{ii,jj};
        if isempty(thisb)
            continue;
        end
        numberOfTrials(ii,jj) = length(thisb.air_puff_r);
    end
end
