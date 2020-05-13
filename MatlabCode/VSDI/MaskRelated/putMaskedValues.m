function data = putMaskedValues (data,mask,values)
maskI = find(mask);
for ii = 1:size(data,3)
    temp = data(:,:,ii);
    temp(maskI) = values(:,ii);
    data(:,:,ii) = temp;
end

