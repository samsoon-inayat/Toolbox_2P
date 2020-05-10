function ei = placeCellProperties(ei,stimMarker,ow)

if ~exist('ow','var')
    ow = 0;
end

for ii = 1:length(ei)
    aa = ii;
    % extract rasters from ei
    for pp = 1:length(ei{aa}.plane)
        thispFolder = ei{aa}.plane{pp}.folder;
        disp(sprintf('Processing context %s',thispFolder));
        contexts = ei{ii}.plane{pp}.contexts;
        for jj = 1:length(contexts)
            thisContext = contexts(jj);
            Rsm = thisContext.rasters;
            if isempty(Rsm)
                continue;
            end
            try
                cmdTxt = sprintf('Rs = Rsm.%sD;',stimMarker);
                eval(cmdTxt);
            catch
                continue;
            end
    %         cmdTxt = sprintf('Rs = Rsm.%sD;',stimMarker);
    %         eval(cmdTxt);
            fileName = fullfile(thispFolder,sprintf('placeProperties_%d_%s.mat',jj,stimMarker));
            if exist(fileName,'file') && ow == 0;
                temp = load(fileName);
                placeProps = temp.placeProps;
            else
                [pcw,pcc] = getProps(Rs,ei{ii}.b.belt_length);
                placeProps.pcw = pcw;
                placeProps.pcc = pcc;
                save(fileName,'placeProps');
            end
            Rs.placeProps = placeProps;
            cmdTxt = sprintf('Rsm.%sD = Rs;',stimMarker);
            eval(cmdTxt);
            thisContext.rasters = Rsm;
            contexts(jj) = thisContext;
        end
        ei{ii}.plane{pp}.contexts = contexts;
    end
end
display('Done!!!');


function [pcws,pccs] = getProps(Rs,BL)

trials = 3:size(Rs.rasters,1);
for ii = 1:size(Rs.rasters,3)
    raster = Rs.rasters(trials,:,ii);
    mSig = nanmean(raster);
    props = place_cell_properties(mSig,raster,'cm_per_bin',BL/size(Rs.rasters,2));
    pcws(ii) = props.width;
    pccs(ii) = props.center;
end
