function [within,dataT_var_names,xlabels] = make_within_table(var_names,nc,nwf)
nfac = length(var_names);
if ~exist('nwf','var')
    nwf = 1;
end
wvar1 = NaN(nc*nwf,1);
wvar2 = wvar1;
ind = 1;
for ii = 1:length(var_names)
    FL(ii) = var_names{ii}(1);
end
for ii = 1:nc
    for jj = 1:nwf
        wvar1(ind) = ii;
        wvar2(ind) = jj;
        if nfac == 1
            dataT_var_names{ind} = sprintf('%s%d',FL(1),ii);
            xlabels{ind} = sprintf('%s%d',FL(1),ii);
        end
        if nfac == 2
            dataT_var_names{ind} = sprintf('%s%d_%s%d',FL(1),ii,FL(2),jj);
            xlabels{ind} = sprintf('%s%d-%s%d',FL(1),ii,FL(2),jj);
        end
        ind = ind + 1;
    end
end

if nfac == 2
    within = array2table([wvar1 wvar2]);
end
if nfac == 1
    within = array2table([wvar1]);
end
within.Properties.VariableNames = var_names;
for ii = 1:length(var_names)
    cmdTxt = sprintf('within.%s = categorical(within.%s);',var_names{ii},var_names{ii});
    eval(cmdTxt);
end







