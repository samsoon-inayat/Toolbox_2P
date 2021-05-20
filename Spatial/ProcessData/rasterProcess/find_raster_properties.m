function ei = find_raster_properties(ei)
allContexts = contextDefinitions;
for aa = 1:length(ei)
    ptei = ei{aa};
    tei = ptei;
    for pp = 1:length(ei{aa}.plane)
        thispFolder = tei.plane{pp}.folder;
        disp(thispFolder);
        tei.folders.thispFolder = thispFolder;
%         tei.areCells = find(tei.plane{pp}.tP.areCells);
%         tei.deconv = tei.plane{pp}.tP.deconv;
%         tei.ops1 = tei.plane{pp}.tP.ops;
%         tei.tP = tei.plane{pp}.tP;
        tplane = tei.plane{pp};
        contexts = tplane.contexts;
        for contextNumber = 1:length(contexts)
            thisContext = contexts(contextNumber);
            rasters = thisContext.rasters;
            fields = fieldnames(rasters);
            for ff = 1:length(fields)
                cmdTxt = sprintf('thisSubset = rasters.%s;',fields{ff});
                eval(cmdTxt);
                Rs = thisSubset.sp_rasters1;
                Dur = thisSubset.duration1;
                MIs = NaN(size(Rs,3),1);
                parfor rr = 1:size(Rs,3)
                    MIs(rr) = info_metrics_S_onlyMI_2(Rs(:,:,rr),[],4,Dur,0);
                end
                thisSubset.MIs = MIs;
                cmdTxt = sprintf('rasters.%s = thisSubset;',fields{ff});
                eval(cmdTxt);
            end
            thisContext.rasters = rasters;
            contexts(contextNumber) = thisContext;
        end
        ptei.plane{pp}.contexts = contexts;
%         ptei.b.belt_length = temp.belt_length;
    end
    ei{aa} = ptei;
end

