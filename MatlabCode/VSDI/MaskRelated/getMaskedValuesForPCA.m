function varargout = getMaskedValuesForPCA (data,mask)
[r c] = find(mask==1);
maskI = sub2ind(size(mask),r,c);
values = [];
for ii = 1:size(data,3)
    temp = data(:,:,ii);
    values = [values temp(maskI)];
    if ~isreal(data)
        maxvalues(ii) = max(abs(temp(maskI)));
    else
        maxvalues(ii) = max(temp(maskI));
    end
end

if nargout == 1
    varargout{1} = values';
end

if nargout == 2
    varargout{1} = values';
    varargout{2} = maxvalues;
end
