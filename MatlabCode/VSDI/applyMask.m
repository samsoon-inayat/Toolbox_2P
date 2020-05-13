function data = applyMask (data,mask,varargin)
for ii = 1:size(data,3)
    data(:,:,ii) = data(:,:,ii) .* mask;
end
if nargin == 3
    mdata = varargin{1};
else
    mdata = min(data(:));
end
imask = ~mask * mdata;
for ii = 1:size(data,3)
    data(:,:,ii) = data(:,:,ii) + imask;
end
