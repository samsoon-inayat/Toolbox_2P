function get_durations_context(ei)

if iscell(ei)
    ei = ei{1};
end
b = ei.b;

if isfield(b,'air_puff_f')
    for ii = 1:length(b.air_puff_r)
        air_puff_durations(1,ii) = b.ts(b.air_puff_f(ii)) - b.ts(b.air_puff_r(ii));
    end
    air_puff_durations
end

if isfield(b,'stim_f')
    for ii = 1:length(b.stim_r)
        stim_durations(1,ii) = b.ts(b.stim_f(ii)) - b.ts(b.stim_r(ii));
    end
    stim_durations
end