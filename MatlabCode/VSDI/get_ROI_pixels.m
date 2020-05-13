function pixels = get_ROI_pixels (center,numberOfPixels)

if isstruct(center)
    names = fieldnames(center);
    for ii = 2:length(names)
        cmdText = sprintf('thisCenter = center.%s;',names{ii});
        eval(cmdText);
        thispixels(:,1) = (thisCenter(1)-numberOfPixels(1)):(thisCenter(1)+numberOfPixels(1));
        thispixels(:,2) = (thisCenter(2)-numberOfPixels(2)):(thisCenter(2)+numberOfPixels(2));
        pixels{ii-1} = thispixels;
    end
else
    thispixels(:,1) = (center(1)-numberOfPixels(1)):(center(1)+numberOfPixels(1));
    thispixels(:,2) = (center(2)-numberOfPixels(2)):(center(2)+numberOfPixels(2));
    pixels{1} = thispixels;
end
pixels = pixels';