function [BL1,onsetsO,offsetsO,inds] = get_belt_length(ei)

b = ei.b;
d = b.dist;

onsets = b.photo_sensor_f(1:(length(b.photo_sensor_f)-1));
offsets = b.photo_sensor_f(2:end);

BL = d(offsets)-d(onsets);

mBL = max(BL);

inds = find(BL > (0.75*mBL));

onsetsO = onsets(inds);
offsetsO = offsets(inds);

BL1 = d(offsetsO) - d(onsetsO);

