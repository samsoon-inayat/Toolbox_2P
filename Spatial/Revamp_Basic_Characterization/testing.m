function testing(ei)

if ~exist('ei','var')
    ei = evalin('base','data(1)');
end

contexts = getContexts(ei);
b = ei{1}.b;
onsets = contexts(1).markers.motionOnsetsOffsets_onsets;
offsets = contexts(1).markers.motionOnsetsOffsets_offsets;
plotMarkers(b,onsets,offsets,100,'motion')
n = 0;