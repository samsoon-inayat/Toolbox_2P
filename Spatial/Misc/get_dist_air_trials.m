function [BL1] = get_belt_length(ei)

b = ei.b;

d = b.dist;

onsets = b.air_puff_r;
offsets = b.air_puff_f;

if length(onsets) == length(offsets)
    BL1 = d(offsets)-d(onsets);
else
    if length(onsets) > length(offsets)
        onsets(end) = [];
        BL1 = d(offsets)-d(onsets);
    end
    n = 0;;
end

