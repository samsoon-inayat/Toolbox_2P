function data = applyMask (data,mask)
for ii = 1:size(data,3)
    data(:,:,ii) = data(:,:,ii) .* mask;
end
mdata = min(data(:));
imask = ~mask * mdata;
for ii = 1:size(data,3)
    data(:,:,ii) = data(:,:,ii) + imask;
end
