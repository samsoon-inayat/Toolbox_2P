function rastersO = find_resp_mdata(rasters,Dur,timeBefore,effective_sampling_rate,titlestr,ei,tRs,isCell)

binWidths = evalin('base','binWidths');

if isnan(timeBefore)
    rastersO.mdata = [];
    rastersO.resp = findResponsiveRasters(rasters,Dur,[],tRs,isCell);
    rastersO.rasters = rasters;
    rastersO.title = titlestr;
    number_of_columns = size(rasters,2);
    totalTime = number_of_columns * binWidths(2);
    rastersO.xs = linspace(0,totalTime,number_of_columns);
    return;
end
    

number_of_columns = size(rasters,2);
column_index = round(timeBefore * effective_sampling_rate);
cis = [];
cis(1,:) = [1 column_index (number_of_columns-column_index+1)];
column_index = cis - 1; column_index(1) = []; column_index(3) = number_of_columns;
cis(2,:) = column_index;
resp = findResponsiveRasters(rasters,Dur,cis);
number_of_columns = size(rasters,2);
totalTime = number_of_columns/effective_sampling_rate;
rastersO.xs = linspace(0,totalTime,number_of_columns);
% % mdata.xs = round(mdata.xs - mdata.xs(column_index),1);
mdata.cis = cis(1,2:3);
rastersO.mdata = mdata;
rastersO.resp = resp;
rastersO.rasters = rasters;
rastersO.title = titlestr;