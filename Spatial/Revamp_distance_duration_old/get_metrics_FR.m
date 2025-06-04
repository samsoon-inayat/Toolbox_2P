function met_vals = get_metrics_FR(out,var_name,trialconcat,nshuffles)

noofbMI = 10; %nshuffles = 5;
MI_fun = @(x,y,noofbMI,nshuffles) calc_metric_MI(x,y,noofbMI,nshuffles);
PC_fun = @(x,y,nshuffles) calc_metric_PC(x,y,nshuffles);

field_names = fieldnames(out);
for ii = 1:length(field_names)
    varname = field_names{ii};
    cmdTxt = sprintf('clear %s;',varname);eval(cmdTxt);
    cmdTxt = sprintf('%s = out.%s;',varname,varname);eval(cmdTxt);
end
met_vals = {};
for an = 1:5
    for cn = 1:3
        for ap = 1:2
            timecc = atimecc{an,cn,ap}; distcc = adistcc{an,cn,ap}; speedcc = aspeedcc{an,cn,ap}; FRcc = aFRcc{an,cn,ap}; trialcc = atrialcc{an,cn,ap};
            idx_us = strfind(var_name,'_');
            var1 = var_name(1:(idx_us-1)); var2 = var_name((idx_us+1):end);
            cmdTxt = sprintf('var1v = %scc;',var1);eval(cmdTxt); cmdTxt = sprintf('var2v = %scc;',var2);eval(cmdTxt);
            num_neurons = size(var1v,2);
            % grouped_var1v = accumarray(trialcc, var1v, [], @(x) {x}); grouped_var2v = accumarray(trialcc, var2v, [], @(x) {x});
            if strcmp(trialconcat,'trials')
                if strcmp(mfun,'MI')
                    thisvar = arrayfun(@(i) MI_fun(grouped_var1v{i}, grouped_var2v {i},noofbMI,nshuffles), 1:length(grouped_var1v));
                end
                if strcmp(mfun,'PC')
                    thisvar = arrayfun(@(i) PC_fun(grouped_var1v{i}, grouped_var2v {i},nshuffles), 1:length(grouped_var1v));
                end
            end
            if strcmp(trialconcat,'concatenate')
                clear metric_vals
                if nshuffles == 0
                    metric_valsM = NaN(num_neurons,1);
                    metric_valsP = NaN(num_neurons,1);
                else
                    metric_valsM = NaN(num_neurons,3);
                    metric_valsP = NaN(num_neurons,3);
                end
                thisFR = var1v; thissig = var2v;
                parfor nn = 1:num_neurons
                    metric_valsM(nn,:) = MI_fun(thisFR(:,nn), thissig,noofbMI,nshuffles);
                    metric_valsP(nn,:) = PC_fun(thisFR(:,nn), thissig,nshuffles);
                end
                met_vals{an,cn,ap}.MI = metric_valsM;
                met_vals{an,cn,ap}.PC = metric_valsP;
            end
            end
        end
    end
end