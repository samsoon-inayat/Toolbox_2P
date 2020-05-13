function mask = makeMask(sizeImg,pixels)

mask = zeros(sizeImg);
[mesX mesY] = meshgrid(pixels(:,1),pixels(:,2));
mesX = reshape(mesX,1,numel(mesX));
mesY = reshape(mesY,1,numel(mesY));
% maskI = sub2ind(sizeImg,pixels(:,2),pixels(:,1));
maskI = sub2ind(sizeImg,mesY,mesX);
mask(maskI) = ones(size(maskI));
