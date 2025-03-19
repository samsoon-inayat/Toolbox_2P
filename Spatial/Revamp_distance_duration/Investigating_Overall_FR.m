function Investigating_Overall
%%
%Listing all the variables that I have related to the animal I have may be 
% . Time, Distance, Speed, and neural
%firing rates.
% experimental configurations include C3, C4, and C5 as well as presence or
% absence of stimuli include the air on and off phases and also the light
% stimulus in C4. Now first just to quantify motion, I can look at the time
% distance and speed and how they are changing with respect to air and
% light stimuli. This would be to quantify motion. The repeated measures
% ANOVA at the animal level will provide me with the information whether
% there was a certain effect of categorical variables which may be
% consistent across animals if I find significant differences.

udataT = evalin('base','udataT'); %Time Bins
udataD = evalin('base','udataD'); %Dist Bins
ei = evalin('base','ei');
configurations = {'C3','C4','C5'};
air_phases = {'ON','OFF'};
bin_types = {'time_bin','dist_bin'};
variable_combs = {'time_dist','time_speed','dist_speed','FR_time','FR_dist','FR_speed'};
variable_combs = {'time_dist','time_speed','dist_speed'};
variable_combs = {'FR_time','FR_dist','FR_speed'};
azMI_vals = []; azPC_vals = [];
ow = [0 0] ;
for an = 1:5
    zMI_vals = []; zPC_vals = [];
    for bn = 1:length(bin_types)
        if bn == 1
            data_an = udataT{an};
        else
            data_an = udataD{an};
        end
        for cn = 1:length(configurations)
            for ap = 1:length(air_phases)
                out = get_continuous_variables(data_an,air_phases{ap},configurations{cn}); % concatenate the variables, time, distance, speed, and neural firing rates for trials in air phase (on or off) and configuration (C3, C4, or C5)
                % find the pairwise metrics including mutual information and
                % correlation
                for vn = 1:length(variable_combs)
                    tvar = variable_combs{vn};
                    spos = strfind(tvar,'_'); var1 = tvar(1:(spos-1)); var2 = tvar((spos+1):end);
                    cmdTxt = sprintf('var1v = out.%s;',var1); eval(cmdTxt); cmdTxt = sprintf('var2v = out.%s;',var2); eval(cmdTxt);
                    params = {'no_of_bins_for_MI',10,'no_of_shuffles_for_norm',500,'animal_info',ei{an},'overwrite_processing',ow,'air_phase',air_phases{ap},...
                        'configuration',configurations{cn},'variables',variable_combs{vn},'bin_type',bin_types{bn},'trial_type','concatenated'};
                    met = myMetrics(var1v,var2v,params); MIs = met.MI; PCs = met.PC;
                    valsMI = sum(isnan(MIs(:,2)))/size(MIs,1); valsPC = sum(isnan(PCs(:,2)))/size(MIs,1);
                    zMI_vals = [zMI_vals valsMI]; zPC_vals = [zPC_vals valsPC];
                end
            end
        end
    end
    azMI_vals(an,:) = zMI_vals; azPC_vals(an,:) = zPC_vals;
end

n = 0;
% (Intercept):Bin_Type:Air_Phase:Corr_Type [F(2,8) = 28.96, p < 0.001, η2 = .14] <--

%%
[within,dvn,xlabels,awithinD] = make_within_table({'Bin_Type','Conf_Num','Air_Phase','Corr_Type'},[2,3,2,3]);
data_matrix = azPC_vals;
dataT = make_between_table({data_matrix},dvn);
ra = RMA(dataT,within,{0.05,{''}});
ra.ranova
print_for_manuscript(ra)





