function avar = get_metrics_one_phase(out,var_name,mfun,trialconcat,phase)

noofbMI = 10; nshuffles = 0;
MI_fun = @(x,y,noofbMI,nshuffles) calc_metric_MI(x,y,noofbMI,nshuffles);
PC_fun = @(x,y,nshuffles) calc_metric_PC(x,y,nshuffles);

field_names = fieldnames(out);
for ii = 1:length(field_names)
    varname = field_names{ii};
    cmdTxt = sprintf('clear %s;',varname);eval(cmdTxt);
    cmdTxt = sprintf('%s = out.%s;',varname,varname);eval(cmdTxt);
end
avar = [];
for an = 1:5
    anvar = [];
    for cn = 1:3
        for ap = phase
            timecc = atimecc{an,cn,ap}; distcc = adistcc{an,cn,ap}; speedcc = aspeedcc{an,cn,ap}; FRcc = aFRcc{an,cn,ap}; trialcc = atrialcc{an,cn,ap};
            idx_us = strfind(var_name,'_');
            var1 = var_name(1:(idx_us-1)); var2 = var_name((idx_us+1):end);
            cmdTxt = sprintf('var1v = %scc;',var1);eval(cmdTxt); cmdTxt = sprintf('var2v = %scc;',var2);eval(cmdTxt);
            grouped_var1v = accumarray(trialcc, var1v, [], @(x) {x}); grouped_var2v = accumarray(trialcc, var2v, [], @(x) {x});
            if strcmp(trialconcat,'trials')
                if strcmp(mfun,'MI')
                    thisvar = arrayfun(@(i) MI_fun(grouped_var1v{i}, grouped_var2v {i},noofbMI,nshuffles), 1:length(grouped_var1v));
                end
                if strcmp(mfun,'PC')
                    thisvar = arrayfun(@(i) PC_fun(grouped_var1v{i}, grouped_var2v {i},nshuffles), 1:length(grouped_var1v));
                end
            end
            if strcmp(trialconcat,'concatenate')
                if strcmp(mfun,'MI')
                    thisvar = MI_fun(var1v, var2v,noofbMI,nshuffles);
                end
                if strcmp(mfun,'PC')
                    thisvar = PC_fun(var1v, var2v,nshuffles);
                end
            end
            anvar = [anvar thisvar];
        end
    end
    avar = [avar;anvar];
end