function frames = load_tif_files(fns,varargin)

h = waitbar(1/length(fns),'Loading image frames');


if iscell(fns)
    for ii = 1:length(fns)
        waitbar(ii/length(fns),h,sprintf('Loading image frames %d/%d',ii,length(fns)));
        frames(:,:,ii) = double(imread(fns{ii}));
    end
else
    imginfo = imfinfo(fns);
    for ii = 1:length(imginfo)
        waitbar(ii/length(fns),h,sprintf('Loading image frames %d/%d',ii,length(fns)));
        frames(:,:,ii) = double(imread(fns,ii));
    end
end
close(h);