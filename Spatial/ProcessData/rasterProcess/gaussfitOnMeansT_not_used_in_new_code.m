function ei = gaussfitOnMeansT(ei,stimMarker,ow)

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
                cmdTxt = sprintf('Rs = Rsm.%sT;',stimMarker);
                eval(cmdTxt);
            catch
                continue;
            end
            fileName = fullfile(thispFolder,sprintf('gaussfitomT_%d_%s.mat',jj,stimMarker));
            if exist(fileName,'file') && ow == 0
                temp = load(fileName);
                mean_raster_fits = temp.mean_raster_fits;
            else
                mean_raster_fits = getFits(Rs,1:size(Rs.rasters,3),3:size(Rs.rasters,1));
                save(fileName,'mean_raster_fits');
            end
            Rs.gauss_fit_on_mean = mean_raster_fits;
            cmdTxt = sprintf('Rsm.%sT = Rs;',stimMarker);
            eval(cmdTxt);

    %         cmdTxt = sprintf('Rs = Rsm.%sT;',stimMarker);
    %         eval(cmdTxt);
    %         fileName = fullfile(ei{ii}.folders.thispFolder,sprintf('gaussfitomT_%d_%s.mat',jj,stimMarker));
    %         if exist(fileName,'file')
    %             temp = load(fileName);
    %             mean_raster_fits = temp.mean_raster_fits;
    %         else
    %             mean_raster_fits = getFits(Rs,1:size(Rs.rasters,3),3:size(Rs.rasters,1));
    %             save(fileName,'mean_raster_fits');
    %         end
    %         Rs.gauss_fit_on_mean = mean_raster_fits;
    %         cmdTxt = sprintf('Rsm.%sT = Rs;',stimMarker);
    %         eval(cmdTxt);

            thisContext.rasters = Rsm;
            contexts(jj) = thisContext;
        end
        ei{ii}.plane{pp}.contexts = contexts;
    end
end
display('Gaussian Fitting Done');