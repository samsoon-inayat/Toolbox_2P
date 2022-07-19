function behaviorProcessor_1_AD_Ctrl
mData = evalin('base','mData');
training_data_C = evalin('base','training_data_C');
% {'173706';'183761';'183745';'183628';'183762'};
selRows_C = [4 9 8 6 10]; selCols = [1:3];
data_C = get_training_data(training_data_C,selRows_C);

training_data_A = evalin('base','training_data_A');
% {'183227';'183228';'183329';'001567';'001569'};
selRows_A = [4 5 6 1 2]; selCols = [1:3];
data_A = get_training_data(training_data_A,selRows_A);
n = 0;
%%

figure(100);clf;
training_data = training_data_C; sel_rows = selRows_C;
training_data = training_data_A; sel_rows = selRows_A;
trials = [1 30]; timePre = 1; timePost = 1; T = [0 20*60];
for ani = 1:length(sel_rows)
    an = sel_rows(ani);
    for dn = 1:3
        subplot(3,1,dn);
        b = training_data.bs{an,dn}; %t1 = b.air_puff_r(trials(1)) - round(1e6 * timePre/b.si); t2 = b.air_puff_f(trials(2)) - round(1e6 * timePost/b.si);
        t1 = find(b.ts-T(1)>0,1,'first')-1; t2 = find(b.ts-T(2)>0,1,'first')-1;
        t1 = b.air_puff_r(1); t2 = b.air_puff_f(end);
        xs = b.ts(t1:t2); ys1 = b.fSpeed(t1:t2); ys2 = b.air_puff_raw(t1:t2);
        plot(xs,ys1);%hold on; plot(xs,ys2*max(ys1)); axis off;
        format_axes_b(gca);
    end
    title(ani)
pause;
end

%% Raw speeds 
while 1
    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[10 3 4.6 2],'RowsCols',[6 1],'spaceRowsCols',[0.1 -0.02],'rightUpShifts',[0.07 0.1],'widthHeightAdjustment',[10 -110]);
    MY = 50; ysp = 3; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'Percent of Cells'};
    stp = 0.25*magfac; widths = ([3.7 1.3 1.3 1.3 1.3 0.5 0.5 0.5]+0.25)*magfac; gap = 0.16*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{''});
    
    training_data = training_data_C; an = 9; trials = [1 30]; timePre = 1; timePost = 1; T = [0 15*60];
    for dn = 1:3
        axes(ff.h_axes(dn,1));
        b = training_data.bs{an,dn}; %t1 = b.air_puff_r(trials(1)) - round(1e6 * timePre/b.si); t2 = b.air_puff_f(trials(2)) - round(1e6 * timePost/b.si);
        t1 = find(b.ts-T(1)>0,1,'first')-1; t2 = find(b.ts-T(2)>0,1,'first')-1;
        xs = b.ts(t1:t2); ys1 = b.fSpeed(t1:t2); ys2 = b.air_puff_raw(t1:t2);
        plot(xs,ys1);hold on; plot(xs,ys2*max(ys1)); axis off;
        format_axes_b(gca);xlim([0 340])
    end
    
    
    %%
    break;
end


%%
while 1
    HiFD_C = find_HiFD_here(data_C,20);
    HiFD_A = find_HiFD_here(data_A,20);
    %%
    ms_C = find_mean_speed(data_C,20);
    ms_A = find_mean_speed(data_A,20);
    %%
    break
end
%%
while 1
    var_C = []; var_A = [];
    for dd = 1:3
%         var_C = [var_C squeeze(ms_C(:,dd,:))];
%         var_A = [var_A squeeze(ms_A(:,dd,:))];
        var_C = [var_C squeeze(HiFD_C(:,dd,2:end))];
        var_A = [var_A squeeze(HiFD_A(:,dd,2:end))];
    end
    
   
    [within,dvn,xlabels,withinD] = make_within_table({'Days','Trials'},[3,19]); awithinD = withinD;
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

%%
while 1
    figure(100);clf
    plot(1:57,mean(var_C)); hold on;
    plot(1:57,mean(var_A)); 
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
