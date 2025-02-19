function mbl = get_mean_belt_length(ei_C,protocol)
an = 1;
if ~isempty(strfind(protocol,'10'))
    if iscell(ei_C)
        for an = 1:length(ei_C)
            for cn = 1:4
                onsets = ei_C{an}.plane{1}.contexts(cn).markers.air_onsets;
                offsets = ei_C{an}.plane{1}.contexts(cn).markers.air_offsets;
                b = ei_C{an}.b;
                belt_lengths(cn,:) = b.dist(offsets)-b.dist(onsets);
            end
            mbl{an} = (mean(belt_lengths,2))';
        end
    else
    for cn = 1:4
        onsets = ei_C{an}.plane{1}.contexts(cn).markers.air_onsets;
        offsets = ei_C{an}.plane{1}.contexts(cn).markers.air_offsets;
        b = ei_C{an}.b;
        belt_lengths(cn,:) = b.dist(offsets)-b.dist(onsets);
    end
    mbl = (mean(belt_lengths,2))';
    end
    return 
end

if ~isempty(strfind(protocol,'15'))
    if iscell(ei_C)
        for an = 1:length(ei_C)
            for cn = 1:3
                try
                    onsets = ei_C{an}.plane{1}.contexts(cn).markers.air_onsets;
                    offsets = ei_C{an}.plane{1}.contexts(cn).markers.air_offsets;
                catch
                    onsets = ei_C{an}.plane{1}.contexts(cn+2).markers.air_onsets;
                    offsets = ei_C{an}.plane{1}.contexts(cn+2).markers.air_offsets;
                end
                b = ei_C{an}.b;
                belt_lengths(cn,:) = b.dist(offsets)-b.dist(onsets);
            end
            mbl{an} = (mean(belt_lengths,2))';
        end
    else
    for cn = 1:3
        onsets = ei_C{an}.plane{1}.contexts(cn+2).markers.air_onsets;
        offsets = ei_C{an}.plane{1}.contexts(cn+2).markers.air_offsets;
        b = ei_C{an}.b;
        belt_lengths(cn,:) = b.dist(offsets)-b.dist(onsets);
    end
    mbl = (mean(belt_lengths,2))';
    end
    return 
end