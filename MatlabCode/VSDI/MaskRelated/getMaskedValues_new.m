function values = getMaskedValues (data,mask)
maskI = find(mask);
values = zeros(length(maskI),size(data,3));
for ii = 1:size(data,3)
    temp = data(:,:,ii);
    values(:,ii) = temp(maskI);
end
