function out = get_the_binned_data(udata,bintype,binwidth)
dist_bins = 0:binwidth:1000; % set a large number of bins
time_bins = 0:binwidth:1000;
configurations = {'C3','C4','C5'};
air_phases = {'ON','OFF'};
atimecc = {}; adistcc = {}; aspeedcc = {}; aFRcc = {}; atrialcc = {};  o_atimecc = {}; o_adistcc = {}; o_aspeedcc = {}; o_atrialcc = {}; 
for an = 1:5
    data_an = udata{an};
    field_names = fieldnames(data_an);
    for ii = 1:length(field_names)
        varname = field_names{ii};
        cmdTxt = sprintf('clear %s;',varname);eval(cmdTxt);
        cmdTxt = sprintf('%s = data_an.%s;',varname,varname);eval(cmdTxt);
    end
    firing_rate = data_an.firing_rate;
    % firing_rate = data_an.ca_signal;
    for cn = 1:length(configurations)
        for ap = 1:2
            timecc = []; distcc = []; speedcc = []; FRcc = []; trialcc = []; active_cells = [];
            otimecc = []; odistcc = []; ospeedcc = []; otrialcc = []; oFRcc = [];
            for tn = 1:10
                [an cn ap tn]
                % get time to complete trial
                if ap == 1
                    cmdTxt = sprintf('air_trials = air_trials_on .* %s;',configurations{cn}); eval(cmdTxt);
                else
                    cmdTxt = sprintf('air_trials = air_trials_off .* %s;',configurations{cn}); eval(cmdTxt);
                end
                idx = find(air_trials == tn & frf);
                tts = ts(idx)-ts(idx(1)); tds = ds(idx)-ds(idx(1)); tds(tds<0) = 0; 
                tsp = speed(idx);

                otimecc = [otimecc;tts']; odistcc = [odistcc;tds']; ospeedcc = [ospeedcc;tsp']; otrialcc = [otrialcc;tn*ones(size(tts'))];

                idx_fr = frf_n(idx);
                FR = firing_rate(:,idx_fr);
                if strcmp(bintype,'time')
                    bin_indices = discretize(tts,time_bins);
                end
                if strcmp(bintype,'dist')
                    bin_indices = discretize(tds,dist_bins);
                end
                time_binned = do_the_binning(bin_indices,tts); time_binned = time_binned - time_binned(1);
                dist_binned = do_the_binning(bin_indices,tds); dist_binned = dist_binned- dist_binned(1);
                speed_binned = do_the_binning(bin_indices,tsp);
                timecc = [timecc;time_binned]; distcc = [distcc;dist_binned]; speedcc = [speedcc;speed_binned];
                trialcc = [trialcc;tn*ones(size(time_binned))];
                FR_binned = [];
                parfor neuron_idx = 1:size(firing_rate,1)
                    % Calculate the binned firing rate for each neuron using the frame-based bin indices
                    % FR_binned(neuron_idx, :) = accumarray(bin_indices', FR(neuron_idx, :),[], @mean, NaN)';
                    FR_binned(neuron_idx, :) = do_the_binning(bin_indices,FR(neuron_idx, :))
                end
                FR_binned = FR_binned';
                FRcc = [FRcc;FR_binned];
                oFRcc = [oFRcc FR];
                active_cells = [active_cells sum(FR,2)~=0];
            end
            % figure(100);clf;
            % subplot 311;plot(timecc); subplot 312; plot(distcc); subplot 313; plot(speedcc);
            % pause;
            atimecc{an,cn,ap} = timecc; adistcc{an,cn,ap} = distcc; aspeedcc{an,cn,ap} = speedcc; aFRcc{an,cn,ap} = FRcc; atrialcc{an,cn,ap} = trialcc;
            o_atimecc{an,cn,ap} = otimecc; o_adistcc{an,cn,ap} = odistcc; o_aspeedcc{an,cn,ap} = ospeedcc; o_atrialcc{an,cn,ap} = otrialcc;
            [rasters{an,cn,ap},rx{an,cn,ap},ry{an,cn,ap}] = build_rasters(trialcc,FRcc');
            aactive_cells{an,cn,ap} = sum(oFRcc,2)~=0;
            active_cells_trials{an,cn,ap} = active_cells;
            n = 0;
        end
    end
end

out.atrialcc = atrialcc;
out.atimecc = atimecc;
out.adistcc = adistcc;
out.aspeedcc = aspeedcc;
out.aFRcc = aFRcc;
out.arasters = rasters;
out.arx = rx;
out.ary = ry;
out.active_cells = aactive_cells;
out.active_cells_trials = active_cells_trials;



out.otrial = o_atrialcc;
out.otime = o_atimecc;
out.odist = o_adistcc;
out.ospeed = o_aspeedcc;

function binned = do_the_binning(bin_indices,signal)
maxindex = max(bin_indices);
binned = accumarray(bin_indices',signal,[maxindex,1],@mean,NaN);
binned = fillmissing(binned, 'linear');  % Linear interpolation
                % time_binned = time_binned - time_binned(1);